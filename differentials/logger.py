import logging

logging.basicConfig(level=logging.INFO)

def info(text):
    logging.info(text)

def debug(text):
    logging.debug(text)

def warning(text):
    logging.warning(text)

def set_level_debug():
    logging.getLogger().setLevel(logging.DEBUG)
