import urllib.request
import os
import sys
import subprocess


def download_tree(locDir):
    url = "http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/"\
          "UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"
    treePath = os.path.join(locDir, "public-latest.all.masked.pb.gz")
    urllib.request.urlretrieve(url, treePath)
    return treePath


def convert_tree(locDir):
    print(locDir)
    treePath = os.path.join(locDir, "public-latest.all.masked.pb.gz")
    varCmd = f"matUtils extract -i {treePath} -C lineagePaths.txt"
    sys.stdout.flush()  # force python to flush
    completed = subprocess.run(varCmd, shell=True, executable="/bin/bash",
                               stdout=subprocess.DEVNULL)
    return completed


def get_curated_lineage_data(locDir):
    url2 = "https://raw.githubusercontent.com/outbreak-info/outbreak.info/"\
           "master/web/src/assets/genomics/curated_lineages.json"
    urllib.request.urlretrieve(url2,
                               os.path.join(locDir,
                                            "curated_lineages.json"))


if __name__ == '__main__':
    download_tree()
    get_curated_lineage_data()
    convert_tree()
