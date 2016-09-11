
import logging

from pbcore.io import FastqRecord

class PhasingResult(object):

    _record   = None
    _locus    = None
    _summary  = None
    _subreads = None

    def __init__(self, barcode, locus, record, summary, subreads, isJunk ):
        self._barcode  = "0" if barcode is None else barcode
        self._locus    = locus
        self._record   = self._formatRecord( record )
        self._locus    = locus
        self._summary  = summary
        self._subreads = subreads
        self._isJunk   = isJunk

        self._validateRecord()
        self._validateIsJunk()

    def _formatRecord( self, record ):
        """
        Direct LAA Outputs can clash across loci, so we splice
        in the name of the locus as a workaround
        """
        idParts = record.id.split('_', 1)
        newId = "{0}_Locus{1}_{2}".format(idParts[0], self._locus, idParts[1])
        return FastqRecord(newId, record.sequence, record.quality)

    def _validateRecord( self ):
        if not isinstance( self._record, FastqRecord ):
            raise RuntimeError("Sequence record is not a valid FASTQ!")

    def _validateIsJunk( self ):
        if not isinstance( self._isJunk, bool ):
            raise RuntimeError("IsJunk argument must be boolean!")

    @property
    def record(self):
        return self._record

    @property
    def id(self):
        return self._record.id

    @property
    def sequence(self):
        return self._record.sequence

    @property
    def quality(self):
        return self._record.quality

    @property
    def qualityString(self):
        return self._record.qualityString

    @property
    def barcode(self):
        return self._barcode

    @property
    def locus(self):
        return self._locus

    @property
    def summary(self):
        return self._summary

    @property
    def subreads(self):
        return self._subreads

    @property
    def isJunk(self):
        return self._isJunk
