#################################################################################
# Copyright (c) 2015-2016, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#################################################################################

# Author: Brett Bowman

import logging
import itertools

def getDataSetBarcodes( dataset ):
    bcs = set()

    # If the index has no barcode information, the dataset isn't barcoded
    if not dataset.isBarcoded:
        return []

    # Otherwise return every unique pair of Fwd/Rev barcodes
    for bc in itertools.izip(dataset.index.bcForward,
                             dataset.index.bcReverse):
        bcs.add(bc)
    bcs = ["{0}--{1}".format(f, r) for f, r in sorted(list(bcs)) if f >= 0 if r >= 0]
    return bcs

def getDoBcBarcodes( bcStr ):
    if bcStr is None:
        return []
    pairs = []
    for bc in bcStr.split(','):
        pair = bc.split('--')
        assert len(pair) == 2, "Invalid barcode string: {0}".format(bc)
        assert pair[0].isdigit(), "Invalid barcode index: {0}".format(pair[0])
        assert pair[1].isdigit(), "Invalid barcode index: {0}".format(pair[0])
        pairs.append(bc)
    return pairs

def barcodeIntersection( dsBarcodes, optBarcodes ):
    intersect = []
    for bc in optBarcodes:
        if bc not in dsBarcodes:
            msg = "Invalid Barcode: '{0}' not found in DataSet!".format(bc)
            logging.error( msg )
            raise RuntimeError( msg )
        intersect.append( bc )
    return sorted(intersect)

def getBarcodes( dataset, bcOpt ):
    logging.debug("Scanning the input data for unique barcodes...")
    dsBarcodes  = getDataSetBarcodes( dataset )
    optBarcodes = getDoBcBarcodes( bcOpt )

    if not dsBarcodes:
        logging.warn("No barcode data detected: has this DataSet been barcoded?")

        # If the dataset wasn't barcoded we shouldn't see a doBc option
        if optBarcodes:
            msg = "Invalid Barcode: Dataset is not barcoded!"
            logging.error( msg )
            raise RuntimeError( msg )

        intersection = dsBarcodes
    elif not optBarcodes:
        # If there were no reads specified via doBc, return all barcodes
        intersection = dsBarcodes
    else:
        # Otherwise find the intersection
        intersection = barcodeIntersection( dsBarcodes, optBarcodes )

    if not intersection:
        logging.debug("Input data not barcoded - no barcode-pairs to analyze")
    else:
        logging.debug("Found {0} valid barcode-pair(s) to analyze".format(len(intersection)))

    return intersection
