import subprocess
from typing import List


def run_cmd_and_pipe_stderr(command: str) -> List[str]:
    """Runs the test command and pipes the stderr for further processing.

    E3SM diags uses the Python logging module for logging runs. The Python
    logger uses stderr for streaming, rather than stdout (e.g., print
    statements). To capture information such as file paths (e.g.,
    ``reference_dir``) from the logger, stderr must be piped using
    ``capture_output=True``. Be aware, piping stderr results in logger messages
    not outputting to the console when running tests. The workaround is to
    perform a normal print of the entire piped stderr outputs once testing
    complete.

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
