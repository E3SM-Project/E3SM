"""Logger module for setting up a custom logger."""
import logging
import logging.handlers
import shutil

LOG_FILENAME = "e3sm_diags_run.log"


def custom_logger(name: str, propagate: bool = True) -> logging.Logger:
    """Sets up a custom logger.

    Parameters
    ----------
    name : str
        Name of the file where this function is called.
    propagate : bool, optional
        Whether to propagate logger messages or not, by default True.

    Returns
    -------
    logging.Logger
        The logger.

    Examples
    ---------
    Detailed information, typically of interest only when diagnosing problems:

    >>> logger.debug("")

    Confirmation that things are working as expected:

    >>> logger.info("")

    An indication that something unexpected happened, or indicative of some
    problem in the near future:

    >>> logger.warning("")

    The software has not been able to perform some function due to a more
    serious problem:

    >>> logger.error("")

    Similar to ``logger.error()``, but also outputs stack trace:

    >>> logger.exception("", exc_info=True)

    A serious error, indicating that the program itself may be unable to
    continue running:

    >>> logger.critical("")
    """
    log_format = "%(asctime)s [%(levelname)s]: %(filename)s(%(funcName)s:%(lineno)s) >> %(message)s"
    log_filemode = "w"

    # Setup
    logging.basicConfig(
        format=log_format,
        filename=LOG_FILENAME,
        filemode=log_filemode,
        level=logging.INFO,
    )
    logging.captureWarnings(True)

    logger = logging.getLogger(name)
    logger.propagate = propagate

    # Console output
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(log_format))
    logger.addHandler(console_handler)

    return logger


def move_log_to_prov_dir(results_dir: str):
    """Moves the e3sm diags log file to the provenance directory.

    This function should be called at the end of the diagnostic run to capture
    all console outputs.

    Parameters
    ----------
    results_dir : str
        The results directory for the run.
    """
    provenance_dir = f"{results_dir}/prov/{LOG_FILENAME}"

    logger = custom_logger(__name__)
    logger.info(f"Log file saved in {provenance_dir}")

    shutil.move(LOG_FILENAME, provenance_dir)
