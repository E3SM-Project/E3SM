import ast
import configparser
import os
import shutil
import subprocess
from typing import List

from PIL import Image, ImageChops, ImageDraw

from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter import SET_TO_PARAMETERS
from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.meridional_mean_2d_parameter import MeridionalMean2dParameter
from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
from tests.integration.config import TEST_ROOT_PATH

logger = custom_logger(__name__)


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


def _get_test_params() -> List[CoreParameter]:
    param = CoreParameter()
    ts_param = AreaMeanTimeSeriesParameter()

    m2d_param = MeridionalMean2dParameter()
    m2d_param.plevs = [
        200.0,
        500.0,
    ]
    z2d_param = ZonalMean2dParameter()
    z2d_param.plevs = [
        200.0,
        300.0,
    ]

    enso_param = EnsoDiagsParameter()
    enso_param.test_name = "e3sm_v1"

    params = [param, ts_param, m2d_param, z2d_param, enso_param]

    return params


def _convert_cfg_to_param_objs(cfg_path: str) -> List[CoreParameter]:
    """Convert diagnostic cfg entries to parameter objects.

    NOTE: ast.literal_eval is not considered "safe" on untrusted data.
    The reason why it is used is because `configparser.ConfigParser`
    doesn't work well with parsing Python types from strings in
    `.cfg` files, resulting in things such as nested strings or string
    representation of lists. Since we are only calling literal_eval on
    `.cfg` files hosted in this repo, there is minimal risk here.

    Returns
    -------
    List[CoreParameter]
        A list of CoreParameter objects, one for each diagnotic set.
    Notes
    -----
    This function seems to be a duplicate of `CoreParser._get_cfg_paramters()`'.
    """
    config = configparser.ConfigParser()
    config.read(cfg_path)
    params = []

    for set_name in config.sections():
        param = SET_TO_PARAMETERS[set_name]()

        for option in config.options(set_name):
            val = config.get(set_name, option)
            val = ast.literal_eval(val)

            setattr(param, option, val)

        params.append(param)

    return params


def _count_images(directory: str):
    """Count the number of images of type file_type in directory"""
    count = 0

    for _, __, files in os.walk(directory):
        for f in files:
            if f.endswith("png"):
                count += 1

    return count


def _compare_images(
    mismatched_images: List[str],
    image_name: str,
    path_to_actual_png: str,
    path_to_expected_png: str,
) -> List[str]:
    # https://stackoverflow.com/questions/35176639/compare-images-python-pil

    actual_png = Image.open(path_to_actual_png).convert("RGB")
    expected_png = Image.open(path_to_expected_png).convert("RGB")
    diff = ImageChops.difference(actual_png, expected_png)

    diff_dir = f"{TEST_ROOT_PATH}image_check_failures"
    if not os.path.isdir(diff_dir):
        os.mkdir(diff_dir)

    bbox = diff.getbbox()
    # If `diff.getbbox()` is None, then the images are in theory equal
    if bbox is None:
        pass
    else:
        # Sometimes, a few pixels will differ, but the two images appear identical.
        # https://codereview.stackexchange.com/questions/55902/fastest-way-to-count-non-zero-pixels-using-python-and-pillow
        nonzero_pixels = (
            diff.crop(bbox)
            .point(lambda x: 255 if x else 0)
            .convert("L")
            .point(bool)
            .getdata()
        )
        num_nonzero_pixels = sum(nonzero_pixels)
        logger.info("\npath_to_actual_png={}".format(path_to_actual_png))
        logger.info("path_to_expected_png={}".format(path_to_expected_png))
        logger.info("diff has {} nonzero pixels.".format(num_nonzero_pixels))
        width, height = expected_png.size
        num_pixels = width * height
        logger.info("total number of pixels={}".format(num_pixels))
        fraction = num_nonzero_pixels / num_pixels
        logger.info("num_nonzero_pixels/num_pixels fraction={}".format(fraction))

        # Fraction of mismatched pixels should be less than 0.02%
        if fraction >= 0.0002:
            mismatched_images.append(image_name)

            simple_image_name = image_name.split("/")[-1].split(".")[0]
            shutil.copy(
                path_to_actual_png,
                os.path.join(diff_dir, "{}_actual.png".format(simple_image_name)),
            )
            shutil.copy(
                path_to_expected_png,
                os.path.join(diff_dir, "{}_expected.png".format(simple_image_name)),
            )
            # https://stackoverflow.com/questions/41405632/draw-a-rectangle-and-a-text-in-it-using-pil
            draw = ImageDraw.Draw(diff)
            (left, upper, right, lower) = diff.getbbox()
            draw.rectangle(((left, upper), (right, lower)), outline="red")
            diff.save(
                os.path.join(diff_dir, "{}_diff.png".format(simple_image_name)),
                "PNG",
            )

    return mismatched_images
