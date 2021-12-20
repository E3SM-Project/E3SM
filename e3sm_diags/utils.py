import subprocess
from typing import List


def run_command_and_get_stderr(command: str) -> List[str]:
    """Runs the test command and captures the stderr for further processing.

    The Python logger uses to stderr for streaming, not stdout. To capture
    file paths, stderr must be piped using ``capture_output=True``. This means
    logger messages won't output in the console during the test run until it
    completes (or fails), which is done via a print statement.

    Parameters
    ----------
    command : str
        The test command.

    Returns
    -------
    List[str]
        List of strings from stderr, decoded with "utf-8".

    Notes
    -----
    If capture_output is true, stdout and stderr will be captured. When used,
    the internal Popen object is automatically created with stdout=PIPE and
    stderr=PIPE. The stdout and stderr arguments may not be supplied at the
    same time as capture_output. If you wish to capture and combine both streams
    into one, use stdout=PIPE and stderr=STDOUT instead of capture_output.

    References
    ----------
    https://docs.python.org/3/library/subprocess.html
    """
    print("\nRunning tests, please wait for log output.")
    proc: subprocess.CompletedProcess = subprocess.run(
        command.split(), capture_output=True
    )
    stderr = proc.stderr.decode("utf-8").splitlines()

    print(*stderr, sep="\n")
    return stderr
