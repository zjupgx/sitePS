import argparse
import pandas as pd
import subprocess
import os
import shutil
from .preprocess2csv import preprocess
from .sites_add_PS import site_add_PS
import concurrent.futures
from typing import Optional, List


def download(sp_id):
    print(f"Downloading {sp_id} Proteomes")
    baseurl = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/"
    UP_df = pd.read_csv(path + "/proteomes_all_2023_12_01.tsv")
    try:
        UP_line = UP_df[UP_df['Tax_ID'] == sp_id].iloc[0]
        UP_name = UP_line.iloc[0]
        UP_genus = UP_line.iloc[2].title()
    except Exception:
        print(f"==Error in {sp_id} Proteomes, maybe {sp_id} not in our uniprot name database or we can not find the url"
              f" of {sp_id}. Please manually download this proteomes to Basic_101species_fa==")
        return None
    gz_name = UP_name + '_' + str(sp_id) + '.fasta.gz'
    url = baseurl + UP_genus + '/' + UP_name + '/' + gz_name
    try:
        subprocess.run(f"wget -c {url}", shell=True, stderr=subprocess.PIPE, check=True)
    except Exception:
        print(f"==Error in {sp_id} Proteomes, maybe {sp_id} not in our uniprot name database or we can not find"
              f" the url of {sp_id}. Please manually download this proteomes to Basic_101species_fa==")
        return None
    subprocess.run(f"gzip -d {gz_name}", shell=True)
    subprocess.run(f"sed -i '/>/s/>.*|\\(.*\\)|.*/>\\1/' {gz_name[:-3]}", shell=True)
    new_name = str(sp_id) + '.fa'
    subprocess.run(f"mv {gz_name[:-3]} {new_name}", shell=True)


def blastdb(species):
    sp_name = species[:-3]
    species = path + "/Basic_101species_fa/" + species
    # print(f"Make {sp_name} balstdb")
    subprocess.run(f"mkdir {path + '/Basic_101species_db/' + sp_name}", shell=True)
    db_path = path + "/Basic_101species_db/" + sp_name + "/" + sp_name
    try:
        subprocess.run(f"makeblastdb -in {species} -dbtype prot -out {db_path}", shell=True, check=True)
    except Exception:
        print(f"==Making blastdb failed in {sp_name}==")


def blastp(species, evalue, query):
    sp_name = species[:-3]
    db_path = path + "/Basic_101species_db/" + sp_name + "/" + sp_name
    out_file = path + "/Blastpout/" + sp_name + '.out'
    print(f"Starting blastp {sp_name}")
    try:
        subprocess.run(f"blastp -evalue {evalue} -query {query} -db {db_path} -out {out_file} -num_threads 1 -outfmt 6", shell=True, check=True)
    except Exception:
        print(f"==Making blastp alignment failed in {sp_name}==")


def sitePS(
    args_query: str,
    args_output: str,
    args_evalue: str = '1e-4',
    args_threads: str = '64',
    args_add: str = None,
    args_delete: Optional[List[int]] = None
):
    dir_fas = "Basic_101species_fa"
    Git_url = "https://github.com/Dabaixian/sitePS/archive/refs/heads/main.zip"
    global path
    path = os.getcwd()

    if args_add:
        with open(args_add) as add:
            file_add = add.readlines()
        Tax_id_adds = [int(item.split(',')[0]) for item in file_add]
        subprocess.run(f"cat {args_add} >> Basic_101species_name.list", shell=True, check=True)
    if args_delete:
        Tax_id_deletes = [int(item) for item in args_delete]

    if os.path.exists("Blastpout"):
        shutil.rmtree("Blastpout")
        os.mkdir("Blastpout")
    else:
        os.mkdir("Blastpout")
    if os.path.exists("Basic_101species_db"):
        if args_add or args_delete:
            shutil.rmtree("Basic_101species_db")
            os.mkdir("Basic_101species_db")
    else:
        os.mkdir("Basic_101species_db")

    with open(args_query) as q:
        query = q.readlines()
    query_name = query[0][1:].rstrip('\n')
    if '9606' not in query_name:
        subprocess.run(f"sed -i 's/>/>9606_/' {args_query}", shell=True, check=True)

    if not os.path.exists("Basic_101species_fa"):
        print("First use, downloading default 101 species name list")
        subprocess.run(f"wget -c {Git_url}", shell=True, check=True)
        subprocess.run(f"unzip main.zip", shell=True, check=True)
        subprocess.run(f"mv sitePS-main/* .", shell=True, check=True)
        subprocess.run(f"rm -r main.zip sitePS-main", shell=True, check=True)
        print("Downloading default 101 species name")
        subprocess.run("mkdir Basic_101species_fa", shell=True, check=True)
        os.chdir("Basic_101species_fa")
        with open(path + '/Basic_101species_name.list') as f:
            file = f.readlines()
        file.pop(0)
        default_id = [int(item.split(',')[0]) for item in file]
        if args_add:
            default_id.extend(Tax_id_adds)
            subprocess.run(f"cat {args_add} >> Basic_101species_name.list", shell=True, check=True)
        if args_delete:
            default_id = list(set(default_id) - set(Tax_id_deletes))
        if 9606 in default_id:
            default_id.remove(9606)

        with concurrent.futures.ThreadPoolExecutor(int(args_threads)) as executor:
            futures = [executor.submit(download, item) for item in default_id]
            concurrent.futures.wait(futures)
            for completed_future in concurrent.futures.as_completed(futures):
                try:
                    result = completed_future.result()
                except Exception as e:
                    print(f"Exception: {e}")
        species_fa = os.listdir()
        with concurrent.futures.ThreadPoolExecutor(int(args_threads)) as executor:
            futures = [executor.submit(blastdb, item) for item in species_fa]
            concurrent.futures.wait(futures)
            for completed_future in concurrent.futures.as_completed(futures):
                try:
                    result = completed_future.result()
                except Exception as e:
                    print(f"Exception: {e}")
        os.chdir(path)
    else:
        if args_delete:
            for i in Tax_id_deletes:
                delete_fa = str(i) + '.fa'
                subprocess.run(f"rm Basic_101species_fa/{delete_fa}", shell=True, check=True)
        if args_add:
            os.chdir("Basic_101species_fa")
            with concurrent.futures.ThreadPoolExecutor(int(args_threads)) as executor:
                futures = [executor.submit(download, item) for item in Tax_id_adds]
                concurrent.futures.wait(futures)
                for completed_future in concurrent.futures.as_completed(futures):
                    try:
                        result = completed_future.result()
                    except Exception as e:
                        print(f"Exception: {e}")

        os.chdir(path)
        species_fa = os.listdir("Basic_101species_fa")
        if args_add or args_delete:
            with concurrent.futures.ThreadPoolExecutor(int(args_threads)) as executor:
                futures = [executor.submit(blastdb, item) for item in species_fa]
                concurrent.futures.wait(futures)
                for completed_future in concurrent.futures.as_completed(futures):
                    try:
                        result = completed_future.result()
                    except Exception as e:
                        print(f"Exception: {e}")

    mafft_name = query_name + "_mafft.fa"
    mafft_out = mafft_name + "s"
    if os.path.exists(mafft_name):
        os.remove(mafft_name)
    if os.path.exists(mafft_out):
        os.remove(mafft_out)
    with concurrent.futures.ThreadPoolExecutor(int(args_threads)) as executor:
        futures = [executor.submit(blastp, item, args_evalue, args_query) for item in species_fa]
        concurrent.futures.wait(futures)
        for completed_future in concurrent.futures.as_completed(futures):
            try:
                result = completed_future.result()
            except Exception as e:
                print(f"Exception: {e}")
    for i in species_fa:
        with open(path + "/Blastpout/" + i.replace("fa", "out")) as f:
            file = f.readlines()
        try:
            target_gene = file[0].split('\t')[1]
        except IndexError:
            continue
        specie_name = i.split('.')[0]
        subprocess.run(f"awk -v RS=\">\" -v seq=\"{target_gene}\" \'index($0, seq) {{print \">\" $0}}\' Basic_101species_fa/{i} | sed '/^$/d' | sed 's/>/>{specie_name}_/' >> {mafft_name}", shell=True, check=True)
    subprocess.run(f"cat {args_query} >> {mafft_name}", shell=True, check=True)
    try:
        subprocess.run(f"mafft  --globalpair --maxiterate 1000 --thread {args_threads} {mafft_name} > {mafft_out}", shell=True, check=True)
    except Exception as e:
        print(e)
    tmp_sites_csv = preprocess("Basic_101species_name.list", mafft_out)
    site_add_PS(tmp_sites_csv, args_output)
    subprocess.run(f"rm {tmp_sites_csv}", shell=True, check=True)


