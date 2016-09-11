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

from __future__ import absolute_import
import argparse
import os
import os.path as op
import sys

PRESETS = ["classI", "fiveLoci", "gendx"]

options = argparse.Namespace()

def parseOptions():
    """
    Parse and sanity-check the options
    """
    desc = "Run Long Amplicon Analysis v2 indepedently on different loci and combine the results"
    parser = argparse.ArgumentParser(description=desc, add_help=True)

    def canonicalizedFilePath(path):
        return op.abspath(op.expanduser(path))

    def checkInputDirectory(path):
        if not op.isdir(path):
            parser.error("Could not find input directory: {0}".format( path ))

    def checkInputFile(path):
        if not op.isfile(path):
            parser.error("Input file not found: {0}".format( path ))

    def checkOutputDirectory(path):
        if not op.isdir(path):
            try:
                os.makedirs(path)
            except:
                parser.error("Could not create output directory: {0}".format( path ))
        if not os.access( path, os.W_OK ):
            parser.error("Output directory is not writable: {0}".format( path ))

    def parseList( value ):
        if value is None:
            return None
        return value.split(',')

    def parseDict( value ):
        if value is None:
            return None
        return {p.split(':')[0]:p.split(':')[1] for p in value.split(',')}

    def parseDictOfLists( value ):
        if value is None:
            return None
        return {p.split(':')[0]:p.split(':')[1:] for p in value.split(',')}

    basics = parser.add_argument_group("Basic required options")
    basics.add_argument(
        "referenceDirectory",
        type=canonicalizedFilePath,
        help="The input directory of per-locus reference sequences")
    basics.add_argument(
        "inputFilename",
        type=canonicalizedFilePath,
        help="The filename of the query BAM or DataSet")
    basics.add_argument(
        "-o", "--outputDirectory",
        default="loci_analysis",
        metavar="STRING",
        type=canonicalizedFilePath,
        help="The output folder for combined results")
    basics.add_argument(
        "--verbose", "-v",
        dest="verbosity",
        action="count",
        help="Set the verbosity level.")
    basics.add_argument(
        "--quiet",
        dest="quiet",
        action="store_true",
        help="Turn off all logging, including warnings")
    basics.add_argument(
        "--rngSeed",
        type=int,
        metavar="INT",
        default=42,
        help="RNG seed, modulates which reads are chosen when they exceed the number needed. Default = 42")
    basics.add_argument(
        "-n", "--nproc",
        type=int,
        metavar="INT",
        default=1,
        help="The number of processors to be used")

    barcoding = parser.add_argument_group("Barcode Options")
    barcoding.add_argument(
        "--doBc",
        metavar="STRING",
        help="Comma-separated list of barcode-index pairs to analyze, in the form '0--0'. Default = All")
    barcoding.add_argument(
        "--minBarcodeScore",
        type=int,
        metavar="INT",
        default=0,
        help="Minimum average barcode score to require of subreads. Default = 0")

    filtering = parser.add_argument_group("Data Filtering Options")
    filtering.add_argument(
        "-l", "--minLength",
        type=int,
        metavar="INT",
        default=3000,
        help="Minimum length of input reads. Default = 3000")
    filtering.add_argument(
        "-L", "--maxLength",
        type=int,
        metavar="INT",
        default=0,
        help="Minimum length of input reads, set <1 to disable. Default = 0")
    filtering.add_argument(
        "-s", "--minReadScore",
        type=float,
        metavar="FLOAT",
        default=0.75,
        help="Minimum read score of input reads. Default = 0.75")
    filtering.add_argument(
        "--minSnr",
        type=float,
        metavar="FLOAT",
        default=3.75,
        help="Minimum SNR of input reads. Default = 3.75")

    filtering = parser.add_argument_group("Coarse Clustering Options")
    filtering.add_argument(
        "-r", "--maxReads",
        type=int,
        metavar="INT",
        default=1000,
        help="Minimum length of input reads. Default = 1000")
    filtering.add_argument(
        "-c", "--maxClusteringReads",
        type=int,
        metavar="INT",
        default=250,
        help="Minimum length of input reads. Default = 250")
    filtering.add_argument(
        "--skipRate",
        type=float,
        metavar="FLOAT",
        default=0.0,
        help="Skip some high-scoring alignments to disperse the cluster more. Default = 250")

    locus = parser.add_argument_group("Locus Options")
    locus.add_argument(
        "--doLoci",
        metavar="STRING",
        type=parseList,
        help="Comma-separated list of loci to analyze. Default = All")
    locus.add_argument(
        "--ignoreLoci",
        metavar="STRING",
        type=parseList,
        help="Comma-separated list of loci to ignore. Default = None")
    locus.add_argument(
        "--combineLoci",
        metavar="STRING",
        type=parseDictOfLists,
        help="Comma-separated list of loci to combine, in form 'NewName:LocusA:LocusB'. "
             "Useful for capturing good reads associated with the wrong loci. Default = None")

    per_locus = parser.add_argument_group("Per-Locus Options",
        "Locus-level options over-ride global values set by other options. "
        "For locus-level options that require values, they are specified in "
        "the form 'Locus:Value', e.g. A:3000.  Multiple values can be specified "
        "for different alleles separated by commas, e.g. A:3000,B:3200.")
    per_locus.add_argument(
        "--minLengthByLocus",
        metavar="STRING",
        type=parseDict,
        help="Per-locus minimum read length. Default = 3000")
    per_locus.add_argument(
        "--maxLengthByLocus",
        metavar="STRING",
        type=parseDict,
        help="Per-locus maximum read length. Default = 0")
    per_locus.add_argument(
        "--maxReadsByLocus",
        metavar="STRING",
        type=parseDict,
        help="Per-locus maximum number of reads. Default = 1000")
    per_locus.add_argument(
        "--maxClusteringReadsByLocus",
        metavar="STRING",
        type=parseDict,
        help="Per-locus maximum number of clustering reads. Default = 250")
    per_locus.add_argument(
        "--minReadScoreByLocus",
        metavar="STRING",
        type=parseDict,
        help="Per-locus minimum read score of input reads. Default = 0.75")
    per_locus.add_argument(
        "--minSnrByLocus",
        metavar="STRING",
        type=parseDict,
        help="Per-locus minimum SNR of input reads. Default = 3.75")

    presets = parser.add_argument_group("Preset Design Options",
        "Though LociAnalysis is designed to support any arbitrary "
        "combination of amplicon targets and sizes, a couple of presets "
        "are available for common HLA applications.")
    presets.add_argument(
        "--classI",
        dest="classI",
        action="store_true",
        help="Set defaults for full-length HLA A,B,C")
    presets.add_argument(
        "--fiveLoci",
        dest="fiveLoci",
        action="store_true",
        help="Set defaults for full-length HLA A,B,C and the active sites of DQB1,DRB1")
    presets.add_argument(
        "--gendx",
        dest="gendx",
        action="store_true",
        help="Set defaults for the GenDx NGSgo kit, containing A,B,C,DQA,DQB,DPA,DPB,DRB1 and DRB345")

    parser.parse_args(namespace=options)

    # Check that we don't have multiple competing presets
    optDict = vars(options)
    for i in range(len(PRESETS)-1):
        fst = PRESETS[i]
        for snd in PRESETS[i+1:]:
            if optDict[fst] and optDict[snd]:
                parser.error("Contradictory Options: {0} and {1} cannot both be True".format(fst, snd))

    # Validate expected inputs and output directory
    checkInputDirectory(options.referenceDirectory)
    checkInputFile(options.inputFilename)
    checkOutputDirectory(options.outputDirectory)
