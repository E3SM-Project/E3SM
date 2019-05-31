"""
Utilities that are used by many different viewers.
"""

import os
import shutil
import datetime
import acme_diags
from bs4 import BeautifulSoup

CURRENT_TIMESTAMP = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def _copy_acme_logo(root_dir):
    """
    Copy over e3sm_logo.png to root_dir/viewer.
    """
    src_pth = os.path.join(acme_diags.INSTALL_PATH, "e3sm_logo.png")
    dest_path = os.path.join(root_dir, "viewer")
    if not os.path.isfile(os.path.join(dest_path, "e3sm_logo.png")):
        shutil.copy(src_pth, dest_path)


def _get_acme_logo_path(root_dir, html_path):
    """
    Based of the root dir of the viewer and the current
    dir of the html, get the relative path of the E3SM logo.
    """
    _copy_acme_logo(root_dir)

    # When current_dir = myresults-07-11/viewer/index.html, the image is in
    #   myresults-07-11/viewer/viewer/e3sm_logo.png
    # So there's no need to move some number of directories up.
    # That's why we have - 3
    relative_dir = html_path.replace(root_dir + '/', '')
    dirs_to_go_up = len(relative_dir.split('/')) - 1
    pth = os.path.join('.')

    for _ in range(0, dirs_to_go_up):
        pth = os.path.join(pth, '..')
    
    return os.path.join(pth, 'viewer', 'e3sm_logo.png')


def add_header(root_dir, path, parameters):
    """
    Add the header to the html located at the path.
    """
    # We're inserting the following in the body under navbar navbar-default
    # <div id="e3sm-header" style="background-color:#dbe6c5; float:left; width:45%">
    # 	<p style="margin-left:5em">
    # 		<b>E3SM Diagnostics Package [VERSION]</b><br>
    # 		Test: [SOMETHING]<br>
    # 		Reference: [SOMETHING]<br>
    # 		Created: [DATE]<br>
    # 	</p>
    # </div>
    # <div id="e3sm-header2" style="background-color:#dbe6c5; float:right; width:55%">
    # 	<img src="e3sm_logo.png" alt="logo" style="width:201px; height:91px; background-color:#dbe6c5">
    # </div>

    test_name = parameters[0].short_test_name if parameters[0].short_test_name else parameters[0].test_name
    if parameters[0].run_type == 'model_vs_obs':
        ref_name = 'Observation and Reanalysis'
    else:
        ref_name = parameters[0].short_ref_name if parameters[0].short_ref_name else parameters[0].ref_name            

    soup = BeautifulSoup(open(path), "lxml")
    old_header = soup.find_all("nav", "navbar navbar-default")
    if len(old_header) is not 0:
        old_header[0].decompose()

    header_div = soup.new_tag(
        "div", id="e3sm-header", style="background-color:#dbe6c5; float:left; width:45%")
    p = soup.new_tag("p", style="margin-left:5em")

    bolded_title = soup.new_tag("b")
    bolded_title.append("E3SM Diagnostics Package {}".format(acme_diags.__version__))
    bolded_title.append(soup.new_tag("br"))
    p.append(bolded_title)

    p.append("Test: {}".format(test_name))
    p.append(soup.new_tag("br"))

    p.append("Reference: {}".format(ref_name))
    p.append(soup.new_tag("br"))

    p.append("Created: {}".format(CURRENT_TIMESTAMP))

    header_div.append(p)
    soup.body.insert(0, header_div)

    logo_path = _get_acme_logo_path(root_dir, path)
    img_div = soup.new_tag("div", id="e3sm-header2",
                           style="background-color:#dbe6c5; float:right; width:55%")
    img = soup.new_tag("img", src=logo_path, alt="logo",
                       style="width:201px; height:91px; background-color:#dbe6c5")
    img_div.append(img)
    soup.body.insert(1, img_div)

    html = soup.prettify("utf-8")
    with open(path, "wb") as f:
        f.write(html)


def h1_to_h3(path):
    """
    Change any <h1> to <h3> because h1 is just too big.
    """
    soup = BeautifulSoup(open(path), "lxml")
    h1 = soup.find('h1')
    if h1 is None:
        return
    h1.name = 'h3'

    html = soup.prettify("utf-8")
    with open(path, "wb") as f:
        f.write(html)
