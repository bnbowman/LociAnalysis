
import logging
import subprocess
import os.path
import time

from glob import glob


def CallSaWriter( inputFasta ):
    saWriterCmd = ['sawriter', inputFasta]

    logging.debug("Calling sawriter with command line '%s'", ' '.join(saWriterCmd))
    proc = subprocess.Popen(saWriterCmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    logging.debug("Finished running sawriter")

    if proc.returncode != 0:
        logging.error("sawriter failed. Stderr was %s", stderr)
        raise RuntimeError(" exited with returncode {e}"
                           .format(e=proc.returncode))

    return True


class RefDb(object):

    def __init__(self, dbPath, writeSuffixArrays=True):
        logging.info("Building reference database from path '{0}'".format(dbPath))
        tStart = time.time()

        fastas = []
        for suffix in ("fa", "fna", "fasta"):
            fastas.extend(glob(os.path.join(dbPath, "*.{0}".format(suffix))))

        logging.info("Found {0} reference fasta files".format(len(fastas)))

        refs = dict()
        for fasta in fastas:
            bn = os.path.basename(fasta)
            suffixArray = fasta + ".sa"

            loci, _ = os.path.splitext(bn)
            loci = loci.split('.')[0]
            if loci.endswith("_gen"):  # Catch and clean-up the common IMGT genomic suffix
                loci = loci[:-4]

            if not os.path.exists(suffixArray):
                if writeSuffixArrays:
                    CallSaWriter( fasta )
                else:
                    logging.warn("missing suffix array for : '{0}'".format(fasta))
                    suffixArray = None

            if loci in refs.keys():
                msg = "duplicate references for locus '{0}' found".format(loci)
                logging.error(msg)
                raise RuntimeError(msg)

            refs[loci] = (fasta, suffixArray)

        self._refs = refs
        logging.debug("Found references for the following loci : {0}".format(", ".join(sorted(self._refs.keys()))))

        tEnd = time.time()
        logging.info("Finished building reference database in {0}s".format(round(tEnd - tStart, 3)))

    def __iter__(self):
        for locus in sorted(self._refs.keys()):
            yield locus

    def iteritems(self):
        for locus in self:
            yield (locus, self._refs[locus])

    def keys(self):
        for locus in sorted(self._refs.keys()):
            yield locus

    def __getitem__(self, name):
        return self._refs[name]

if __name__ == "__main__":
    import sys

    dbPath = sys.argv[1]

    logging.basicConfig(level=logging.DEBUG)
    db = RefDb( dbPath )
