#! /usr/bin/env python

# Author: Brett Bowman

from __future__ import absolute_import

import logging
import itertools
import time

from pbcore.io import openDataSet

from LociAnalysis.options import (options,
                                  parseOptions)
from LociAnalysis.barcodes import getBarcodes
from LociAnalysis.refdb import RefDb
from LociAnalysis.whitelistdb import WhitelistDb
from LociAnalysis.phaser import LaaPhaser
from LociAnalysis.results import ResultWriter

import LociAnalysis.logger  # Enable TRACE-level logging

class LociAnalysis(object):
    """
    The main driver class for locus-specific Amplicon Analysis
    """
    def __init__(self):
        self._inputFn     = None
        self._inputDs     = None
        self._barcodes    = None
        self._refDb       = None
        self._whitelistDb = None

    def _setupLogging(self):
        if options.quiet:
            logLevel = logging.ERROR
        elif options.verbosity >= 3:
            logLevel = logging.TRACE
        elif options.verbosity == 2:
            logLevel = logging.DEBUG
        elif options.verbosity == 1:
            logLevel = logging.INFO
        else:
            logLevel = logging.WARNING
        logFormat = ">|> %(asctime)s -|- %(levelname)s -|- %(lineno)d -|--|- %(message)s"
        logging.Formatter.converter = time.gmtime
        logging.basicConfig(level=logLevel, format=logFormat)

    def main(self):
        parseOptions()
        self._setupLogging()

        self._inputFn      = options.inputFilename
        self._inputDs      = self._openDataSet(options.inputFilename)
        self._resultWriter = ResultWriter(options.outputDirectory)
        self._barcodes     = getBarcodes(self._inputDs, options.doBc)
        self._refDb        = RefDb(options.referenceDirectory)
        self._whitelistDb  = WhitelistDb(self._refDb, self._inputFn,
                                         dataset=self._inputDs,
                                         combined=options.combineLoci,
                                         nproc=options.nproc)

        if self._barcodes:
            for barcode in self._barcodes:
                self._phaseSample( barcode )
        else:
            self._phaseSample( None )

        self._resultWriter.finalizeSubreadCsv()

    def _phaseSample(self, barcode):
        if barcode is not None:
            logging.info("Processing loci for barcode '{0}'".format(barcode))
        else:
            logging.info("Processing loci for full dataset")

        for locus, locusWl in self._whitelistDb.iteritems():
            if options.doLoci is not None and locus not in options.doLoci:
                logging.debug("Locus '{0}' not specified by the user, skipping".format(locus))
                continue
            if options.ignoreLoci is not None and locus in options.ignoreLoci:
                logging.debug("User elected to ignore Locus '{0}', skipping".format(locus))
                continue

            if barcode is not None:
                logging.info("Phasing locus '{0}' for barcode '{1}'".format(locus, barcode))
            else:
                logging.info("Phasing locus '{0}'".format(locus))

            rngSeed            = self._getOption( locus, "rngSeed" )
            minBarcodeScore    = self._getOption( locus, "minBarcodeScore" )
            minLength          = self._getOption( locus, "minLength" )
            maxLength          = self._getOption( locus, "maxLength" )
            minReadScore       = self._getOption( locus, "minReadScore" )
            minSnr             = self._getOption( locus, "minSnr" )
            maxReads           = self._getOption( locus, "maxReads" )
            maxClusteringReads = self._getOption( locus, "maxClusteringReads" )
            skipRate           = self._getOption( locus, "skipRate" )

            with LaaPhaser(barcode, self._inputFn, locus, nproc=options.nproc,
                           whitelist=locusWl, rngSeed=rngSeed, minLength=minLength,
                           maxLength=maxLength, minReadScore=minReadScore,
                           minSnr=minSnr, maxReads=maxReads, skipRate=skipRate,
                           maxClusteringReads=maxClusteringReads,
                           minBarcodeScore=minBarcodeScore) as phaser:
                for result in phaser:
                    self._resultWriter.writeResult( result )


    def _openDataSet( self, fn ):
        try:
            ds = openDataSet( fn )
        except:
            msg = "Input filename must be a valid BAM or SubreadSet XML"
            logging.error( msg )
            raise RuntimeError( msg )
        return ds

    def _getOption( self, locus, opt ):
        optByLocus = "{0}ByLocus".format(opt)

        optDict = vars(options)
        if optByLocus in optDict and optDict[optByLocus] is not None:
            if locus in optDict[optByLocus]:
                return optDict[optByLocus][locus]
        return optDict[opt]

    @property
    def barcodes(self):
        return self._barcodes


def main():
    LociAnalysis().main()

if __name__ == "__main__":
    main()
