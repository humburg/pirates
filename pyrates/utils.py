"""Utility functions used in other modules"""

import logging

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
