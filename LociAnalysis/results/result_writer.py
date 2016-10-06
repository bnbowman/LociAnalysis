
import csv
import logging
import os
import os.path as op

from collections import defaultdict

from pbcore.io import FastqWriter, FastqRecord

SUMMARY_HEADER = ["BarcodeName", "FastaName", "CoarseCluster", "Phase", "TotalCoverage", "SequenceLength",
                  "PredictedAccuracy", "ConsensusConverged", "NoiseSequence", "IsDuplicate", "DuplicateOf",
                  "IsChimera", "ChimeraScore", "ParentSequenceA", "ParentSequenceB", "CrossoverPosition"]

NumReadsLambda = lambda x: int(x.split('NumReads')[1])

class ResultWriter(object):

    # Primary class-variables
    _directory   = None
    _goodFastq   = None
    _junkFastq   = None
    _summaryCsv  = None
    _subreadRoot = None

    # Secondary class-variables for writing out subread matrices
    _currBarcode = None
    _subreadCols = []
    _subreadData = defaultdict(lambda : defaultdict(int))

    def __init__(self, directory):
        self._directory   = self._validateDirectory( directory )
        self._goodFastq   = self._openFastqWriter( "loci_analysis.fastq" )
        self._junkFastq   = self._openFastqWriter( "loci_analysis_chimeras_noise.fastq" )
        self._summaryCsv  = self._openCsvWriter( "loci_analysis_summary.csv" )
        self._subreadRoot = op.join( self._directory, "loci_analysis_subreads." )

        self._writerSummaryCsvHeader()

    def _validateDirectory( self, directory ):
        if not op.exists( directory ):
            try:
                os.mkdirs( directory )
            except:
                msg = "Could not create result directory: {0}".format( directory )
                logging.error( msg )
                raise RuntimeError( msg )
        if not os.access( directory, os.W_OK):
            msg = "Result directory is not writable: {0}".format( directory )
            logging.error( msg )
            raise RuntimeError( msg )
        return directory

    def _openFastqWriter( self, filename ):
        filepath = op.join( self._directory, filename )
        try:
            writer = FastqWriter( filepath )
        except:
            msg = "Could not open FASTQ output for writing: {0}".format( filename )
            logging.error( msg )
            raise RuntimeError( msg )
        return writer

    def _openCsvWriter( self, filename ):
        filepath = op.join( self._directory, filename )
        try:
            handle = open( filepath, 'wb' )
            writer = csv.writer( handle )
        except:
            msg = "Could not open CSV output for writing: {0}".format( filename )
            logging.error( msg )
            raise RuntimeError( msg )
        return writer

    def _writerSummaryCsvHeader( self ):
        self._summaryCsv.writerow( SUMMARY_HEADER )

    def _checkBarcode( self, barcode ):
        # Check that the barcode is either unset, or matches
        if self._currBarcode is not None:
            if barcode != self._currBarcode:
                msg = "Barcode mismatch ({0} and {1})! Call finalizeSubreadCsv before changing samples".format( barcode, self._currBarcode )
                logging.error( msg )
                raise RuntimeError( msg )
        else:
            self._currBarcode = barcode

    def _writeSummary( self, result ):
        row = [result.barcode, result.id, result.summary["cluster"], result.summary["phase"],
               result.summary["coverage"], str(len(result.sequence)), result.summary["readQuality"],
               result.summary["didConverge"], result.summary["isNoise"], result.summary["isDup"],
               result.summary["dupOf"], result.summary["isChimera"], result.summary["chimeraScore"],
               result.summary["parentA"], result.summary["parentB"], result.summary["crossover"]]
        self._summaryCsv.writerow( row )

    def _addSubreadData( self, result ):
        # All columns (ResultIds) must be unique for us to store data correctly
        if result.id in self._subreadCols:
            msg = "Duplicate Result Id: {0}".format(result.id)
            logging.error( msg )
            raise RuntimeError( msg )
        else:
            self._subreadCols.append( result.id )

        # We store the weight for a subread-result pair with the subread first
        #  for fast row-wise access during writing
        for subread, weight in result.subreads.iteritems():
            self._subreadData[subread][result.id] = weight

    def finalizeSubreadCsv( self ):
        # If we have a barcode, we have subread data that needs to be written out
        if self._currBarcode:
            csv = self._openCsvWriter( self._subreadRoot + self._currBarcode + ".csv" )

            # Header is "SubreadId" followed by result ids in descending order
            cols = ["SubreadId"] + sorted(self._subreadCols, key=NumReadsLambda, reverse=True)
            csv.writerow( cols )

            # Each row is the data for one subread
            for subread, colData in self._subreadData.iteritems():
                row = [subread] + [colData[rid] for rid in cols[1:]]
                csv.writerow( row )

        # Finally, reset subread-related class variables for the next sample
        self._currBarcode = None
        self._subreadCols = []
        self._subreadData = defaultdict(lambda : defaultdict(int))

    def writeResult( self, result ):
        # First check that the barcode for this result is sensible
        self._checkBarcode( result.barcode )

        # If so, write the FASTQ and Summary data to our current handles
        if not result.isJunk:
            self._goodFastq.writeRecord( result.record )
        else:
            self._junkFastq.writeRecord( result.record )
        self._writeSummary( result )

        # Finally, add the subread data to the current store
        self._addSubreadData( result )
