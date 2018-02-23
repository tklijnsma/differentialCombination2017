import logging

# Add level
logging.TRACE = 9
logging.addLevelName(logging.TRACE, 'TRACE')
def trace(self, message, *args, **kws):
    # Yes, logger takes its '*args' as 'args'.
    if self.isEnabledFor(logging.TRACE):
        self._log(logging.TRACE, message, args, **kws) 
logging.Logger.trace = trace
logging.trace = logging.getLogger().trace

FORMAT = '[%(levelname)-7s %(filename)s:%(lineno)s:%(funcName)s] %(message)s'
logging.basicConfig(format=FORMAT, level=logging.INFO)

def set_level_debug():
    logging.getLogger().setLevel(logging.DEBUG)

def set_level_trace():
    logging.getLogger().setLevel(logging.TRACE)

def enable_testmode():
    fmt = logging.Formatter('[testmode]' + FORMAT)
    for h in logging.getLogger().handlers:
        h.setFormatter(fmt)