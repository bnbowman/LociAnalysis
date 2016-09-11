
import logging

TRACE_LEVELV_NUM = 5

logging.addLevelName(TRACE_LEVELV_NUM, "TRACE")

def traceFn(msg, *args, **kwargs):
    if len(logging.root.handlers) == 0:
        logging.basicConfig()
    logging.root.log(TRACE_LEVELV_NUM, msg, *args, **kwargs)

#def traceFn(self, message, *args, **kws):
    #if self.isEnabledFor(TRACE_LEVELV_NUM):
    #    self._log(TRACE_LEVELV_NUM, message, args, **kws)

logging.trace = traceFn
logging.TRACE = TRACE_LEVELV_NUM
