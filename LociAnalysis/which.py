
import os
import os.path

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
