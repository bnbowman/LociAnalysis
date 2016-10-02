
import logging
from subprocess import Popen, PIPE

from LociAnalysis.which import which

def LongAmpliconAnalysisRawString():
    laa = which('laa')
    if not laa:
        raise RuntimeError("laa not on PATH")
    cmd = [laa, "--version"]
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, close_fds=True)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError("`{0}` failed with exit code {1}:\n{2}".format(' '.join(cmd), proc.returncode, proc.stderr.read()))
    return proc.stdout.read()

def LongAmpliconAnalysisVersion( raw ):
    return raw.splitlines()[0].split('|')[0].strip()[4:]

def SmrtAnalysisVersion( laa_version ):
    if laa_version == "2.0.0 (commit 1edd539)":
        return "3.1.0"
    if laa_version == "2.0.0 (commit 0171f76)":
        return "3.1.1"
    elif laa_version.startswith("2.0.0"):
        return "3.1.Roche"
    else:
        return "3.2"

_LONG_AMPLICON_RAW_STR = LongAmpliconAnalysisRawString()
LONG_AMPLICON_VERSION = LongAmpliconAnalysisVersion(_LONG_AMPLICON_RAW_STR)
SMRT_ANALYSIS_VERSION = SmrtAnalysisVersion(LONG_AMPLICON_VERSION)
