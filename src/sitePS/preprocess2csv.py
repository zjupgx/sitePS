import pandas as pd
import subprocess


def preprocess(name_list, mafft_fas):
    with open(name_list) as f:
        file = f.readlines()
    file.pop(0)
    PS = {}
    specie_name = {}
    for line in file:
        line = line.rstrip("\n")
        taxid, ps_num, specie = line.split(",")
        PS[taxid] = ps_num
        specie_name[taxid] = specie

    with open(mafft_fas) as f:
        file = f.readlines()
    dictionary = {}
    bak = {}
    for line in file:
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line
            dictionary[name] = ""
        else:
            dictionary[name] += line
    for key in dictionary.keys():
        bak[key] = ",".join(list(dictionary[key]))

    write_tmp = []
    for key in bak.keys():
        write_tmp.append(PS[key.split("_")[0][1:]] + "," + specie_name[key.split("_")[0][1:]] + "," + bak[key] + "\n")
    write_tmp.sort(key=lambda x: int(x.split(",")[0]), reverse=True)
    with open(mafft_fas[:-10] + "tmp.csv", "w") as w:
        w.writelines(write_tmp)

    df = pd.read_csv(mafft_fas[:-10] + "tmp.csv", index_col=[0, 1], header=None)
    df.index.names = ["", ""]
    df.columns = [item - 1 for item in df.columns.to_list()]
    df = df.loc[:, ~((df.loc[(27, "human")] == "-") & (df.loc[(27, "human")].notna()))]

    df.columns = range(1, len(df.columns) + 1)
    subprocess.run(f"rm {mafft_fas[:-10] + 'tmp.csv'}", shell=True)
    df.to_csv(mafft_fas[:-10] + "sites.csv")
    return mafft_fas[:-10] + "sites.csv"
