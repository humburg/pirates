"""Utility functions used in other modules"""

import logging
import gzip

def console_handler():
    """Create a handler for logging to the console.

    Returns:
        :obj:`logging.Handler`: The requested handler.
    """
    handler = logging.StreamHandler()
    formatter = logging.Formatter('[%(levelname)s] %(name)s - %(message)s')
    handler.setFormatter(formatter)
    return handler

def get_logger(name, level='NOTSET', handlers=None):
    """Create logger with consistent formatting.

    Args:
        name (:obj:`str`): Name to use for logger.
        level (:obj:`str`): Logging level.
        handlers (:obj:`list`): List of logging handlers to add to this logger.

    Returns:
        :obj:`logging.Logger`: The requested logger.
    """
    level = getattr(logging, level, logging.NOTSET)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    if handlers is not None:
        for handler in handlers:
            logger.addHandler(handler)
    return logger

def smart_open(file_name):
    """Choose function to open a file based on its extension

    Args:
        file_name (:obj:`str`): Name of the file to be opened.

    Returns:
        :obj:`function` appropriate for opening a file with the given name.
    """
    access_fun = open
    if file_name.endswith('.gz'):
        access_fun = gzip.open
    return access_fun
