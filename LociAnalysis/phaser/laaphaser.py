
from __future__ import print_function

import csv
import logging
import os
import os.path

from collections import namedtuple
from subprocess import Popen, PIPE
from shutil import rmtree
from tempfile import mkdtemp

from pbcore.io import FastqReader

from LociAnalysis.results import PhasingResult

ILLEGAL_OPTS = set(["--doBc", "--resultFile", "--reportsFile", "--subreadsReportPrefix", "--noChimeraFilter"])

def which(exe, env=os.environ):
    def isExe(p):
        return os.path.isfile(p) and os.access(p, os.X_OK)
    if os.path.split(exe)[0] and isExe(exe):
        return exe
    for d in env.get("PATH", "").split(os.pathsep):
        exePath = os.path.join(d.strip("\""), exe)
        if isExe(exePath):
            return exePath
    return None

class LaaPhaser(object):

    def __init__(self, barcode, dataset, locus=None, nproc=1, **kwargs):
        self._barcode = barcode
        self._dataset = dataset
        self._locus   = locus
        self._nproc   = nproc
        self._tmpdir  = None
        self._records = None
        self._kwargs  = self._validateKwargs(kwargs)

    def _validateKwargs( self, kwargs ):
        invalidOpts = set(kwargs.keys()) & ILLEGAL_OPTS
        if invalidOpts:
            raise RuntimeError("invalid options for laa: '{0}'".format(" ".join(invalidOpts)))
        return kwargs

    def _parseSummaryCsv( self ):
        recData = {}
        with open(os.path.join(self._tmpdir, "amplicon_analysis_summary.csv")) as summ:
            rdr = csv.reader(summ)
            hdr = next(rdr)
            cluster      = hdr.index("CoarseCluster")
            phase        = hdr.index("Phase")
            coverage     = hdr.index("TotalCoverage")
            readQuality  = hdr.index("PredictedAccuracy")
            didConverge  = hdr.index("ConsensusConverged")
            isNoise      = hdr.index("NoiseSequence")
            isDup        = hdr.index("IsDuplicate")
            dupOf        = hdr.index("DuplicateOf")
            isChimera    = hdr.index("IsChimera")
            chimeraScore = hdr.index("ChimeraScore")
            parentA      = hdr.index("ParentSequenceA")
            parentB      = hdr.index("ParentSequenceB")
            crossover    = hdr.index("CrossoverPosition")
            for row in rdr:
                recId = row[1]
                data  = {"cluster":      row[cluster],
                         "phase":        row[phase],
                         "coverage":     row[coverage],
                         "readQuality":  row[readQuality],
                         "didConverge":  row[didConverge],
                         "isNoise":      row[isNoise],
                         "isDup":        row[isDup],
                         "dupOf":        row[dupOf],
                         "isChimera":    row[isChimera],
                         "chimeraScore": row[chimeraScore],
                         "parentA":      row[parentA],
                         "parentB":      row[parentB],
                         "crossover":    row[crossover]}
                recData[recId] = data
        return recData

    def _parseSubreadCsv( self ):
        """
        Convert a CSV matrix into a dictionary-of-dictionaries,
        the outer dictionary indexed by result name, the inner
        dictonary indexed by the subread id, and the value the
        weight associated with that Result-Subread pair if it
        was non-zero
        """
        subreadData = None

        if self._barcode is not None:
            subreadCsv = os.path.join(self._tmpdir, "amplicon_analysis_subreads.{0}.csv".format(self._barcode))
        else:
            subreadCsv = os.path.join(self._tmpdir, "amplicon_analysis_subreads.csv")

        try:
            with open(subreadCsv) as handle:
                reader = csv.reader(handle)
                header = next(reader)
                subreadData = dict((name, dict()) for name in header[1:])
                for row in reader:
                    for i in range(1, len(row)):
                        weight = float(row[i])
                        if weight > 0.0:
                            subreadData[header[i]][row[0]] = weight
        except:
            # If there were parsing errors, return an empty dict
            return {}
        return subreadData

    def _parseSequences( self ):
        """
        Read the two expected output FASTQ files, and combine them into
        a single list of tuples containing the record and which file it
        originated from
        """
        records = []
        for fname, isJunk in (("amplicon_analysis.fastq", False), ("amplicon_analysis_chimeras_noise.fastq", True)):
            for record in FastqReader(os.path.join(self._tmpdir, fname)):
                records.append( (record, isJunk) )
        return records

    def _formatResults( self, sequences, summaryData, subreadData ):
        records = []
        for seqRecord, isJunk in sequences:
            recId    = seqRecord.id
            summary  = summaryData[recId]
            subreads = subreadData[recId]
            records.append(PhasingResult(self._barcode, self._locus, seqRecord, summary, subreads, isJunk))
        return records

    def __enter__(self):
        self._tmpdir = mkdtemp()
        laa = which("laa")
        if not laa:
            raise RuntimeError("laa not on PATH")
        cmd = [laa, "-n", str(self._nproc)]
        if self._barcode is not None:
            cmd.extend([ "--doBc", self._barcode ])
        for key, value in self._kwargs.iteritems():
            cmd.extend([ "--{0}".format(key), str(value) ])
        cmd.append(self._dataset)
        logging.trace("running `{0}` in '{1}'".format(" ".join(cmd), self._tmpdir))
        proc = Popen(cmd, cwd=self._tmpdir, stderr=PIPE, close_fds=True)
        proc.wait()
        if proc.returncode != 0:
            raise RuntimeError("`{0}` failed with exit code {1}:\n{2}".format(' '.join(cmd), proc.returncode, proc.stderr.read()))

        # Parse the various expected output files, then combine them into PhasingResults
        sequences     = self._parseSequences()
        summaryData   = self._parseSummaryCsv()
        subreadData   = self._parseSubreadCsv()
        self._results = self._formatResults( sequences, summaryData, subreadData)

        return self

    def __exit__(self, typ, val, traceback):
        rmtree(self._tmpdir)
        self._tmpdir = None
        self._records = None

    def __iter__(self):
        if self._results is None:
            raise RuntimeError("LaaPhaser is a context object! Use it as such to generate records")
        return iter(self._results)

# for testing purposes
if __name__ == "__main__":
    import sys
    bc, ds = sys.argv[1:]
    logging.basicConfig(level=logging.INFO)
    with LaaPhaser(bc, os.path.abspath(ds)) as lp:
        for r in lp:
            print(repr(r))
