from __future__ import print_function

import subprocess
import argparse
from six.moves.urllib.parse import urlparse
import string
import os.path

#-------------------------------------------------------------------------------

def download_file(url, destination, proxy):

    import subprocess

    print(url, proxy)

    if (proxy == "none"):
        process = subprocess.Popen(["wget", "-O", "%s" %(destination), "%s" %(url)], stdout=subprocess.PIPE)
    else:
        process = subprocess.Popen(["wget", "-O", "%s" %(destination), "-e", "use_proxy=yes", "-e", proxy, "%s" %(url)], stdout=subprocess.PIPE)

    while process.poll() is None:
        line = process.stdout.readline() # This blocks until it receives a newline.
        print(line)
    print(process.stdout.read())

#-------------------------------------------------------------------------------

parser = argparse.ArgumentParser(description='Download data needed for DEMSI test cases.')

parser.add_argument('-p', '--proxytype', dest='proxytype', help='Proxy type', choices=['none','lanl'], default="none")

args = parser.parse_args()

# proxies
proxies = {"none": {"http":  "none",
                    "https": "none"},
           "lanl": {"http":  "http_proxy=http://proxyout.lanl.gov:8080",
                    "https": "https_proxy=http://proxyout.lanl.gov:8080"}}

url = "https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/testdata/MPAS-Seaice_test_dataset_V1/"
proxy = proxies[args.proxytype][urlparse(url).scheme]



manifest = url + "manifest"
download_file(manifest, "manifest", proxy)

manifestFile = open("manifest","r")

filenames = manifestFile.readlines()
for filename in filenames:
    if (not os.path.exists(os.path.dirname(string.strip(filename)))):
        os.makedirs(os.path.dirname(string.strip(filename)))
    download_file(url + string.strip(filename), string.strip(filename), proxy)

manifestFile.close()
    

# environment variable to set
print("Set the MPAS-Seaice data environment variable:")
print("export MPAS_SEAICE_DOMAINS_DIR=%s" %(os.getcwd()))
