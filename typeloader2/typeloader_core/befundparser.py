#!/usr/bin/env python3

from functools import reduce
from collections import defaultdict

# some of the KIR loci have 3 or maybe 4 alleles
# to be prepared, they are in the list, but probably not in the befund.csv
useGenesList = ["HLA-A_1", "HLA-A_2", "HLA-B_1", "HLA-B_2", "HLA-C_1", "HLA-C_2",
                "HLA-DRB1_1", "HLA-DRB1_2",
                "HLA-DRB3_1", "HLA-DRB3_2", "HLA-DRB4_1", "HLA-DRB4_2", "HLA-DRB5_1", "HLA-DRB5_2",
                "HLA-DQB1_1", "HLA-DQB1_2", "HLA-DPB1_1", "HLA-DPB1_2",
                "HLA-E_1", "HLA-E_2",
                "HLA-DPA1_1", "HLA-DPA1_2", "HLA-DQA1_1", "HLA-DQA1_2",
                "HLA-DMA_1", "HLA-DMA_2", "HLA-DMB_1", "HLA-DMB_2",
                "HLA-DOA_1", "HLA-DOA_2", "HLA-DOB_1", "HLA-DOB_2",
                "HLA-F_1", "HLA-F_2", "HLA-G_1", "HLA-G_2", "HLA-H_1", "HLA-H_2",
                "HLA-K_1", "HLA-K_2", "HLA-J_1", "HLA-J_2",
                "MICA", "MICB",
                "KIR2DL1-1", "KIR2DL1-2", "KIR2DL1-3", "KIR2DL1-4",
                "KIR2DL2-1", "KIR2DL2-2", "KIR2DL2-3", "KIR2DL2-4",
                "KIR2DL3-1", "KIR2DL3-2", "KIR2DL3-3", "KIR2DL3-4",
                "KIR2DL4-1", "KIR2DL4-2", "KIR2DL4-3", "KIR2DL4-4",
                "KIR2DL5A-1", "KIR2DL5A-2", "KIR2DL5A-3", "KIR2DL5A-4",
                "KIR2DL5B-1", "KIR2DL5B-2", "KIR2DL5B-3", "KIR2DL5B-4",
                "KIR2DP1-1", "KIR2DP1-2", "KIR2DP1-3", "KIR2DP1-4",
                "KIR2DS1-1", "KIR2DS1-2", "KIR2DS1-3", "KIR2DS1-4",
                "KIR2DS2-1", "KIR2DS2-2", "KIR2DS2-3", "KIR2DS2-4",
                "KIR2DS3-1", "KIR2DS3-2", "KIR2DS3-3", "KIR2DS3-4",
                "KIR2DS4-1", "KIR2DS4-2", "KIR2DS4-3", "KIR2DS4-4",
                "KIR2DS5-1", "KIR2DS5-2", "KIR2DS5-3", "KIR2DS5-4",
                "KIR3DL1-1", "KIR3DL1-2", "KIR3DL1-3", "KIR3DL1-4",
                "KIR3DL2-1", "KIR3DL2-2", "KIR3DL2-3", "KIR3DL2-4",
                "KIR3DL3-1", "KIR3DL3-2", "KIR3DL3-3", "KIR3DL3-4",
                "KIR3DP1-1", "KIR3DP1-2", "KIR3DP1-3", "KIR3DP1-4",
                "KIR3DS1-1", "KIR3DS1-2", "KIR3DS1-3", "KIR3DS1-4"
                ]
old_columns = ["A1", "A2", "B1", "B2", "C1", "C2",
               "DR1", "DR2", "DQ1", "DQ2", "DP1", "DP2"]
changeNamesFor = ["DQ", "DP", "DR"]
patientIdPos = 0
customerPos = 2


def getOtherAlleles(befundFile):
    with open(befundFile) as befundHandle:
        header = befundHandle.readline().strip()
        customer_dic = {}

        if "," in header:
            delimiter = ","
        elif ";" in header:
            delimiter = ";"
        else:
            msg = "Could not figure out delimiter of pretyping file. Please use ',' or ';'!"
            return None, msg

        parts = header.split(delimiter)
        usePos = []
        genes = []
        for pos in range(len(parts)):
            col = parts[pos].strip()
            if col in useGenesList or col in old_columns:
                genes.append(col)
                usePos.append(pos)
        befund = {}
        for line in befundHandle:
            if not (len(line.strip())):
                continue
            line = line.replace("'", "").replace('"', "")  # remove quotes
            parts = line.strip().split(delimiter)
            patient = parts[patientIdPos]
            customer = parts[customerPos]
            customer_dic[patient] = customer
            befund[patient] = defaultdict(list)

            for pos in usePos:
                befundGeneName = genes[usePos.index(pos)]
                geneName = befundGeneName.split("_")[0]
                pretyping = parts[pos].strip()
                changeName = reduce(lambda x, y: x or y,
                                    [befundGeneName.startswith(nameToChange) for nameToChange in changeNamesFor])
                if changeName:  # old class 2 columns
                    geneName = "HLA-" + befundGeneName[:2] + "B1"
                else:
                    if geneName in old_columns:
                        geneName = "HLA-" + befundGeneName[:-1]  # rename gene for old columns
                    elif geneName.startswith("KIR"):
                        geneName = geneName.split("-")[0]
                    elif geneName.startswith("MIC"):
                        pretyping = pretyping.replace("A", "").replace("B",
                                                                       "")  # MIC pretypings are given as GL strings in one col per locus, unlike HLA
                if not pretyping:
                    continue
                if pretyping == "+":
                    befund[patient][geneName].append(pretyping)
                else:
                    for value in pretyping.split("+"):
                        if value:
                            befund[patient][geneName].append(value.strip())

    return befund, customer_dic


if __name__ == '__main__':
    befundFile = r"C:\Daten\local_data\TL_issue_data\Befunde_DR_10.csv"
    befundFile = r"\\nasdd12\daten\data\Typeloader\staging\data_unittest\confirmation_file\Befunde_3DP1_1.csv"
    befund, customer_dic = getOtherAlleles(befundFile)

    for ID in befund:
        for locus in befund[ID]:
            print(locus, befund[ID][locus])
