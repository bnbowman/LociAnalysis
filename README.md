# LociAnalysis

LociAnalysis is a tool for reference-aware, but otherwise de novo analysis
of long amplicon data from SMRT Sequencing.  It uses BLASR to "sort" sequence
data based on which references it is most similar too, and Long Amplicon
Analysis to analyze the data in each bin to generate polished and phased
de novo consensus sequences.

LociAnalysis can be thought of as a wrapper around Long Amplicon Analysis
to allow it to handle larger or more complex samples than it was designed for,
or in effect a form of "Meta Long Amplicon Analysis".  Effort has been made
to perserve option names and functions between Long Amplicon Analysis
and LociAnalysis to facilitate moving work between the two tools.

### Background

Alleles of clinical interest from the genes in the Major Histocompatibility 
Complex (MHC) range from ~3kb in length for genes in Class I to over 16kb
for some alleles from Class II.  Even carefully designed PCR assays that target
only the active sites of the target genes and break longer loci into multiple amplicons,
such as the GenDx NGSgo® HLA Sequencing kit, can still contain 13 primer pairs for 
amplicons from 11 different HLA genes, ranging in size from 400bp to ~6kb.

Long Amplicon Analysis, on the other hand was designed to analyze multiplexed 
amplicon data from a single locus or from a single closely-related gene family, with 
amplicons and alleles within ~1kb in length of each other.  In addition, the initial 
clustering step was not developed with samples containing more than 5 amplicons
in mind.

To handle samples of this complexity and to work around these limitations in the design 
of LAA, LociAnalysis was concieved.  BLASR is used to sort data by sequence
similarity into different bins, and Long Amplicon Analysis is run indepdently on each
bin, optionally with parameters tailored to the contents of that "bin".  This allows 
LociAnalysis scale to handle both more amplicons, and amplicons more divergent
sizes and prevalence rates than LAA can handle directly,

### Input

LociAnalysis requires only two inputs:

1. A SMRT Analysis v3 compliant data
2. A directory of reference FASTA files, one per "bin"

Sequence "bins" for data analysis are taken from the filenames themselves,
such that any sequences most similar to references in "A.fasta" or "A_gen.fasta" 
(the preferred IMGT filename) are assigned to bin "A" for data from the HLA-A
locus.

### Output

LociAnalysis was designed to mimic the output of Long Amplicon Analysis as
closely as possible, only with the prefix "amplicon\_analysis" replaced with
"loci\_analysis"

- loci\_analysis.fastq
- loci\_analysis\_chimeras_noise.fastq
- loci\_analysis\_summary.csv
- loci\_analysis\_subreads.csv

### Requirements

- SMRT Analysis >= v3.1
- virtualenv

LociAnalysis requires a local installation of SMRT Analysis containing LAA v2, 
and a version of BLASR that is POSIX-compliant, which means SMRT Analysis 
v3.1 or greater.

In addition, since LociAnalysis is an experimental tool not officially supported
by PacBio, Python's "virtualenv" is required so that LociAnalysis doesn't
contaminate the SMRT Analysis installation.

### Installation

(1) Activate the local SMRT Analysis Environment
	
	$ module load smrtanlysis/3.1
(2) Create a local virtualenv for LociAnalysis

	$ virtualenv LociAnalysisEnv
(3) Activate the new environment

	$ source LociAnalysisEnv/bin/activate
(4) Install LociAnalysis from GitHub

	$ pip install git+https://github.com/bnbowman/LociAnalysis
### Contact 

Please direct all inquiries to the following e-mail address:
Brett Bowman (bbowman@pacificbiosciences.com)