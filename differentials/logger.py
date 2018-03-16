import logging

standard_format = '[%(levelname)-7s %(filename)-12.12s:%(lineno)-4s:%(funcName)-14.14s] %(message)s'

def set_basic_format():
    # Add level
    logging.TRACE = 9
    logging.addLevelName(logging.TRACE, 'TRACE')
    def trace(self, message, *args, **kws):
        # Yes, logger takes its '*args' as 'args'.
        if self.isEnabledFor(logging.TRACE):
            self._log(logging.TRACE, message, args, **kws) 
    logging.Logger.trace = trace
    logging.trace = logging.getLogger().trace

    logging.basicConfig(format=standard_format, level=logging.INFO)
    set_formatter(standard_format)
    set_level_info()

def set_formatter(fmtstr):
    fmt = logging.Formatter(fmtstr)
    for h in logging.getLogger().handlers:
        h.setFormatter(fmt)    

def set_level_info():
    logging.getLogger().setLevel(logging.INFO)

def set_level_debug():
    logging.getLogger().setLevel(logging.DEBUG)

def set_level_trace():
    logging.getLogger().setLevel(logging.TRACE)

def enable_testmode():
    test_format = '[testmode]' + standard_format
    set_formatter(test_format)