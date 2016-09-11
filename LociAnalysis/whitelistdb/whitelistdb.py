
import re
import os
import os.path as op
import logging
import subprocess
import time
import copy

from collections import defaultdict

from pbcore.io import openDataSet, DataSet

import LociAnalysis.refdb as refdb

NPROC = 1

def CallDataSetCreate( inputBam ):
    outputPath = op.dirname( inputBam )

    outputXml = op.join(outputPath,
                        re.sub("bam$", "subreadset.xml", op.basename( inputBam )))

    datasetCmd = ['dataset', 'create', outputXml, op.abspath(inputBam)]

    logging.trace("Calling dataset-create with command line '%s'", ' '.join(datasetCmd))
    proc = subprocess.Popen(datasetCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    logging.trace("Finished running dataset-create")

    if proc.returncode != 0:
        logging.error("dataset-create failed. Stderr was %s", stderr)
        raise RuntimeError(" exited with returncode {e}"
                           .format(e=proc.returncode))

    return outputXml

def CallBlasr( query, refFn, refSa=None, outputPath=None, name=None, nproc=None ):

    locus = op.basename( refFn ).split('.')[0] if name is None else name

    if query.endswith(".subreadset.xml"):
        outputM1 = op.join(outputPath,
                           re.sub("subreadset.xml$", "{0}.m1".format(locus), op.basename( query )))
    elif query.endswith(".bam"):
        outputM1 = op.join(outputPath,
                           re.sub("bam$", "{0}.m1".format(locus), op.basename( query )))
    else:
        outputM1 = op.join(outputPath, op.basename(query) + ".{0}.m1".format(locus))

    if op.isfile(outputM1) and op.getsize(outputM1) > 0:
        logging.debug("Existing output file for for locus '{0}', skipping alignment".format(locus))
        return outputM1

    blasrCmd = ['blasr', query, refFn, '--bestn', '1', '--out', outputM1, '--fastSDP',
                             '--minSubreadLength', '1000', '--minAlnLength', '1000']

    if nproc is not None:
        blasrCmd.extend(['--nproc', str(nproc)])
    if refSa is not None:
        blasrCmd.extend(["--sa", refSa])

    logging.trace("Calling Blasr with command line '%s'", ' '.join(blasrCmd))
    proc = subprocess.Popen(blasrCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    logging.trace("Finished running Blasr")

    if proc.returncode != 0:
        logging.error("Blasr alignment failed. Stderr was %s", stderr)
        raise RuntimeError(" exited with returncode {e}"
                           .format(e=proc.returncode))

    return outputM1


class WhitelistDb(object):

    _files      = []    # Track temporary files to be deleted
    _scores     = defaultdict(int)
    _reads      = defaultdict(list)
    _loci       = defaultdict(list)
    _whitelists = {}

    def __init__( self, refDb, query, dataset=None, alnDir=None, combined=None, nproc=NPROC ):
        logging.info("Building whitelist database for '{0}'".format(query))
        tStart = time.time()

        self._refDb    = self._getRefDb( refDb )
        self._queryFn  = self._getQuery( query )
        self._queryDs  = self._getDataSet( dataset )
        self._alnDir   = self._getAlnDir( alnDir )
        self._combined = self._getCombinations( combined )
        self._nproc    = nproc

        for locus, (refFn, refSa) in self._refDb.iteritems():
            m1 = CallBlasr(self._queryFn, refFn, refSa, self._alnDir, locus, self._nproc)
            self._updateMapping(m1, locus)
        self._createLociReference()
        self._combineLoci()
        self._writeWhitelists()

        tEnd = time.time()
        logging.info("Finished building whitelist database in {0}s".format(round(tEnd - tStart, 3)))

    def __enter__( self ):
        return self

    def __exit__( self, exception_type, exception_value, traceback ):
        return True

    def _updateMapping(self, m1, locus):
        with open( m1 ) as handle:
            for line in handle:
                parts = line.strip().split()
                query = '/'.join(parts[0].split('/')[:3])
                try:
                    score = abs(int(parts[4]))
                except:
                    continue
                refScore = self._scores[query]
                if score > refScore:
                    self._scores[query] = score
                    self._reads[query] = [locus]
                elif score == refScore:
                    self._reads[query].append( locus )

    def _combineLoci( self ):
        if self._combined is None:
            return
        for name, loci in self._combined.iteritems():
            pool = []
            if name in self._loci.keys():
                pool += self._loci[name]
            for locus in loci:
                pool += self._loci[locus]
            self._loci[name] = pool

    def _getRefDb( self, refDb ):
        if isinstance(refDb, refdb.RefDb):
            return refDb
        else:
            return refdb.RefDb( refDb )

    def _getQuery( self, query ):
        if not op.exists( query ) or not op.isfile( query ):
            msg = "Query must be a valid file!"
            logging.error( msg )
            raise RuntimeError( msg )
        if not query.endswith(".xml") and not query.endswith(".bam"):
            msg = "Query must be either a BAM or DataSet XML!"
            logging.error( msg )
            raise RuntimeError( msg )
        return query

    def _getDataSet( self, ds ):
        if isinstance( ds, DataSet ):
            return ds
        elif ds is not None:
            msg = "Invalid object supplied as DataSet of type: {0}!".format(type(ds))
            logging.error( msg )
            raise RuntimeError( msg )
        else:
            try:
                ds = openDataSet( self._queryFn )
            except:
                msg = "Could not open input file as DataSet: {0}!".format(self._queryFn)
                logging.error( msg )
                raise RuntimeError( msg )
        return ds

    def _getAlnDir( self, alnDir ):
        if alnDir is None:
            alnDir = self._queryFn + "_aln"
        if not op.isdir( alnDir ):
            try:
                os.mkdir( alnDir )
            except:
                msg = "Could not create output directory: {0}".format(alnDir)
                logging.error( msg )
                raise RuntimeError( msg )
        return alnDir

    def _getCombinations( self, combinations ):
        if combinations is None:
            return None
        for name, loci in combinations.iteritems():
            for locus in loci:
                if locus not in self._refDb.keys():
                    msg = "Locus not found ({0})! Can only combine reads from known loci!".format( locus )
                    logging.error( msg )
                    raise RuntimeError( msg )
        return combinations

    def _createLociReference( self ):
        count = 0
        for read, vals in self._reads.iteritems():
            count += 1
            for locus in vals:
                self._loci[locus].append( read )
        logging.debug("Found {0} subreads with at least one good alignment".format(count))

    def _writeWhitelist( self, locus ):
        outputTxt = op.join( self._alnDir, locus + ".subreads.txt" )
        if op.isfile(outputTxt) and op.getsize(outputTxt) > 0:
            logging.debug("Existing locus-specific whitelist file found for '{0}', skipping filtering".format(locus))
            return outputTxt

        with open( outputTxt, 'w' ) as handle:
            for qname in self._loci[locus]:
                handle.write( qname + "\n" )

        return outputTxt

    def _writeWhitelistDataset( self, locus ):
        subreads = self._loci[locus]

        outputXml = op.join( self._alnDir, locus + ".subreadset.xml" )
        if op.isfile(outputXml) and op.getsize(outputXml) > 0:
            logging.debug("Existing locus-specific dataset file found for '{0}', skipping filtering".format(locus))
            return outputXml

        sset = copy.deepcopy( self._queryDs )
        sset.filters.addRequirement(qname=[('=', q) for q in self._loci[locus]])
        sset.write( outputXml )
        logging.debug("Wrote a SubreadSet with {0} whitelisted subreads for locus '{1}'".format(len(subreads), locus))

        return outputXml

    def _writeWhitelists( self ):
        for locus in self._loci.keys():
            self._whitelists[locus] = self._writeWhitelist( locus )
        logging.debug("Found an existing whitelist for {0} loci".format(len(self._whitelists.keys())))

    def _writeWhitelistDatasets( self ):
        for locus in self._loci.keys():
            self._whitelists[locus] = self._writeWhitelistDataset( locus )
        logging.debug("Found a whitelisted SubreadSet for {0} loci".format(len(self._whitelists.keys())))

    def keys(self):
        return sorted(self._whitelists.keys())

    def iteritems(self):
        for locus in self.keys():
            yield (locus, self._whitelists[locus])

    def items(self):
        return [(k,v) for k,v in self.iteritems()]

    def __iter__(self):
        for locus in sorted(self.keys()):
            yield locus

    def __getitem__(self, name):
        return self._whitelists[name]


if __name__ == "__main__":
    import sys

    ref   = sys.argv[1]
    query = sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else NPROC

    logging.basicConfig(level=logging.DEBUG)
    with WhitelistDb( ref, query, nproc=nproc ) as db:
        for locus in db.keys():
            ds = db[locus]
            print locus, ds
