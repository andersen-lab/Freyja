import urllib.request
import os
import sys
import subprocess
import requests


def download_tree(locDir):
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/"\
          "UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"
    treePath = os.path.join(locDir, "public-latest.all.masked.pb.gz")
    urllib.request.urlretrieve(url, treePath)
    return treePath


def convert_tree(loc_dir):
    print(loc_dir)
    tree_path = os.path.join(loc_dir, "public-latest.all.masked.pb.gz")
    var_cmd = f"matUtils extract -i {tree_path} -C lineagePaths.txt"
    sys.stdout.flush()  # force python to flush
    # note we could add check=True here in order to immediately raise if this command fails
    # instead, how we have it, we just log the response code
    return_code = subprocess.run(var_cmd, shell=True, executable="/bin/bash",
                               stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return return_code


def get_curated_lineage_data(locDir):
    url2 = "https://raw.githubusercontent.com/outbreak-info/outbreak.info/"\
           "master/web/src/assets/genomics/curated_lineages.json"
    urllib.request.urlretrieve(url2,
                               os.path.join(locDir,
                                            "curated_lineages.json"))


def get_cl_lineages(locDir):
    # for now, use lineages metadata created using patch
    r = requests.get('https://raw.githubusercontent.com/outbreak-info/' +
                     'outbreak.info/master/curated_reports_prep/lineages.yml')
    # r = requests.get('https://raw.githubusercontent.com/cov-lineages' +
    #                  '/lineages-website/master/data/lineages.yml')
    if r.status_code == 200:
        with open(os.path.join(locDir, 'lineages.yml'), 'w+') as f:
            f.write(r.text)


if __name__ == '__main__':
    download_tree()
    get_curated_lineage_data()
    convert_tree()
