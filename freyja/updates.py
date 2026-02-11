import urllib.request
import os
import sys
import subprocess
import requests
import yaml
from datetime import datetime
import json
import re


def get_pathogen_config(locDir):
    config_path = os.path.join(locDir, 'pathogen_config.yml')

    if not os.path.exists(config_path):  # Check if file exists
        return None  # Return None

    with open(config_path, 'r') as f:
        try:
            pathogen_config = yaml.safe_load(f)
        except yaml.YAMLError:
            raise ValueError('Error loading pathogen update config')

    return pathogen_config


def download_config(locDir):
    url = "https://raw.githubusercontent.com/andersen-lab/"\
          "Freyja/main/freyja/data/pathogen_config.yml"
    lpath = os.path.join(locDir, 'pathogen_config.yml')
    urllib.request.urlretrieve(url, lpath)
    return lpath


def download_tree(locDir):
    url = "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/"\
          "UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"
    treePath = os.path.join(locDir, "public-latest.all.masked.pb.gz")
    urllib.request.urlretrieve(url, treePath)
    return treePath


def get_latest_git_folder_date(org, repo, branch, pathogen_dir):
    api_url = (
        f"https://api.github.com/repos/{org}/{repo}/contents/"
        f"{pathogen_dir}"
        f"?ref={branch}"
    )
    try:
        with urllib.request.urlopen(api_url) as r:
            items = json.loads(r.read().decode())
    except Exception as e:
        print(f"Error fetching GitHub directory listing: {e}")
        return None

    # Match YYYY-MM-DD
    date_pattern = re.compile(r"^\d{4}-\d{2}-\d{2}$")

    dates = [
        item["name"] for item in items
        if date_pattern.match(item["name"])]

    if not dates:
        return None

    # Lexical sort works for YYYY-MM-DD format
    return sorted(dates)[-1]


def download_barcodes(locDir, pathogen='SARS-CoV-2'):
    if pathogen == 'SARS-CoV-2':
        base = (
            "https://raw.githubusercontent.com/"
            "andersen-lab/Freyja/main/freyja/data"
            )

        urllib.request.urlretrieve(f"{base}/usher_barcodes.feather",
                                   os.path.join(locDir,
                                                "usher_barcodes.feather"))
        urllib.request.urlretrieve(f"{base}/last_barcode_update.txt",
                                   os.path.join(locDir,
                                                "last_barcode_update.txt"))
        urllib.request.urlretrieve(f"{base}/lineage_mutations.json",
                                   os.path.join(locDir,
                                                "lineage_mutations.json"))
    else:
        pathogen_config = get_pathogen_config(locDir)
        entry = pathogen_config[pathogen][0]
        barcode_url = entry['barcodes']
        pathogen_name = entry['name']
        parts = barcode_url.split("/")
        org = parts[3]
        repo = parts[4]
        branch = parts[7]
        pathogen_dir = parts[8]
        output_format = parts[10].split(".")[1]
        date_from_github = get_latest_git_folder_date(
            org, repo, branch, pathogen_dir)
        if date_from_github:
            timestamp = date_from_github  # YYYY-MM-DD
        else:
            timestamp = datetime.now().strftime("%Y-%m-%d")  # fallback
        out_csv = os.path.join(locDir,
                               f"{pathogen_name}_barcodes.{output_format}")
        urllib.request.urlretrieve(barcode_url, out_csv)
        print(barcode_url)
        update_file = os.path.join(locDir, "last_barcode_update_other.txt")
        update_line = f"{pathogen_name}:{timestamp}\n"

        lines = []
        if os.path.exists(update_file):
            with open(update_file, "r") as f:
                lines = f.readlines()

        updated = False
        for i, line in enumerate(lines):
            if line.startswith(f"{pathogen_name}:"):
                lines[i] = update_line
                updated = True
                break

        if not updated:
            lines.append(update_line)

        with open(update_file, "w") as f:
            f.writelines(lines)


def convert_tree(loc_dir):
    print(f"Writing updated files to: {loc_dir}")
    tree_path = os.path.join(loc_dir, "public-latest.all.masked.pb.gz")
    var_cmd = f"matUtils extract -i {tree_path} -C lineagePaths.txt"
    sys.stdout.flush()  # force python to flush
    return_code = subprocess.run(var_cmd, shell=True, executable="/bin/bash",
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.PIPE)
    return return_code


def convert_tree_custom(tree_path):
    print(f"Reading custom tree at: {tree_path}")
    var_cmd = f"matUtils extract -i {tree_path} -C lineagePaths.txt"
    sys.stdout.flush()  # force python to flush
    return_code = subprocess.run(var_cmd, shell=True, executable="/bin/bash",
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.PIPE)
    return return_code


def get_curated_lineage_data(locDir, pathogen):
    if pathogen == 'SARS-CoV-2':
        url2 = "https://raw.githubusercontent.com/"\
               "outbreak-info/outbreak.info/"\
               "master/web/src/assets/genomics/curated_lineages.json"
        urllib.request.urlretrieve(url2,
                                   os.path.join(locDir,
                                                "curated_lineages.json"))


def get_cl_lineages(locDir, pathogen='SARS-CoV-2'):
    if pathogen == 'SARS-CoV-2':
        # for now, use lineages metadata created using patch
        r = requests.get('https://raw.githubusercontent.com/outbreak-info/' +
                         'outbreak.info/master/curated_reports_prep/' +
                         'lineages.yml')
        if r.status_code == 200:
            with open(os.path.join(locDir, 'lineages.yml'), 'w+') as f:
                f.write(r.text)
    else:
        pathogen_config = get_pathogen_config(locDir)
        if pathogen_config and pathogen in pathogen_config:
            if isinstance(pathogen_config[pathogen],
                          list) and pathogen_config[pathogen]:
                if 'lineageyml' in pathogen_config[pathogen][0]:
                    url = pathogen_config[pathogen][0]['lineageyml']
                    r = requests.get(url)
                    if r.status_code == 200:
                        file_path = os.path.join(
                            locDir, f"{pathogen_config[pathogen][0]['name']}"
                            "_lineages.yml")
                        with open(file_path, 'w+') as f:
                            f.write(r.text)
        else:
            print(f'Lineage hierarchy not yet set up for {pathogen}')


if __name__ == '__main__':
    download_tree()
    get_curated_lineage_data()
    convert_tree()
