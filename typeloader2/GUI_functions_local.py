#!/usr/bin/env python3
# -*- coding: cp1252 -*-
'''
Created on 11.03.2019

GUI_functions_local.py

contains TypeLoader functionality designed especially for use at DKMS Life Science Lab

@author: Bianca Schoene
'''

# import modules:
import os, copy
from configparser import ConfigParser
from PyQt5.QtWidgets import QMessageBox

from typeloader2 import db_internal, general, typeloader_functions, GUI_login

# ===========================================================
# parameters:

local_config_file = "config_local.ini"


# ===========================================================
# functions:

def read_local_settings(settings, log):
    """reads settings from local config file,
    returns ConfigParser object
    """
    global local_config_file
    if settings["modus"] == "staging":
        local_config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), local_config_file)
    log.info("Reading local settings from {}...".format(local_config_file))

    if os.path.exists(local_config_file):
        cf = ConfigParser()
        cf.read(local_config_file)
    else:
        return False
    return cf


def check_local(settings, log):
    """returns True if this is a DKMS user,
    else returns False
    """
    permission = False
    local_cf = read_local_settings(settings, log)
    # return False, None  # uncomment to treat all users as non-local during debugging and development

    if not local_cf:
        log.debug("No local config found")
        return False, None

    if settings["lab_of_origin"] == local_cf.get("Local", "company_name"):
        permission = True
        log.debug("This is an LSL user")
    else:
        log.debug("This is an LSL-external user")
    return permission, local_cf


def check_nonproductive(settings, log):
    """returns True if this user is not a productive user,
    else returns False
    """
    permission = False
    if settings["modus"] != "productive":
        permission = True
        log.debug("This is a TEST user")
    else:
        log.debug("This is a PROD user")
    return permission


def find_alleles_per_project(project, log, parent=None):
    """finds all alleles of a project
    """
    log.info("Finding alleles of project {}...".format(project))
    query = """select distinct sample_id_int, cell_line_old, local_name, gene, target_allele, partner_allele
            from alleles where project_name = '{}'""".format(project)
    success, data = db_internal.execute_query(query, 6, log, "getting samples from database", "DB error", parent)
    if not success:
        return False, None
    else:
        return True, data


def make_fake_ENA_file(project, log, settings, basis="local_name", parent=None):
    """creates a pseudo-ENA reply file with random ENA accession IDs and a pseudo pretypings file,
    which can be used to create fake-IPD files before ENA has assigned official accession numbers
    """
    success, data = find_alleles_per_project(project, log, parent)
    if not success:
        return False, None, None

    if not data:
        return False, "No alleles found", f"Could not find any alleles for project {project}!"

    # write fake ENA file:
    log.info("Writing fake ENA file...")
    data2 = []
    kir_contained = False

    # create fake ENA content:
    log.info("Assembling fake ENA reply content...")
    ENA_id_map = {}
    ENA_gene_map = {}

    for [sample_id_int, cell_line_old, local_name, mygene, target_allele, partner_allele] in data:
        data2.append((sample_id_int, local_name, mygene, target_allele, partner_allele))  # needed for pretyping file
        fake_ID = typeloader_functions.id_generator(size=8)
        ENA_gene_map[local_name] = mygene
        ENA_id_map[local_name] = fake_ID
        if mygene.startswith("KIR"):
            kir_contained = True

    # write fake pretyping file:
    log.info("Writing fake pretyping-file...")
    fake_file_befunde = os.path.join(settings["login_dir"], "temp", "fake_befunde.csv")
    columns = []

    gene_dic = {}
    for g in ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1",
              "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
              "HLA-DQB1", "HLA-DPB1", "HLA-E",
              "HLA-F", "HLA-G", "HLA-H", "HLA-J", "HLA-K",
              "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "HLA-DQA1", "HLA-DPA1"
              ]:
        mycolumns = [f"{g}_1", f"{g}_2"]
        gene_dic[g] = mycolumns
        columns += mycolumns
    for g in ["MICA", "MICB"]:
        gene_dic[g] = [g]
        columns.append(g)
    kir_columns = []
    for g in ['KIR2DL1', 'KIR2DL2', 'KIR2DL3', 'KIR2DL4', 'KIR2DL5A', 'KIR2DL5B',
              'KIR2DP1',
              'KIR2DS1', 'KIR2DS2', 'KIR2DS3', 'KIR2DS4', 'KIR2DS5',
              'KIR3DL1', 'KIR3DL2', 'KIR3DL3',
              'KIR3DP1', 'KIR3DS1']:
        mycolumns = []
        for i in range(1, 5):
            mycolumns.append(f'{g}-{i}')
        gene_dic[g] = mycolumns
        kir_columns += mycolumns

    if not mygene in gene_dic:  # not one of our standard-genes
        gene_dic[mygene] = ["{}_1".format(mygene), "{}_2".format(mygene)]
        columns += gene_dic[mygene]

    default_dic = {}
    if kir_contained:
        columns += kir_columns
    for col in columns:
        if col.startswith("HLA"):
            default_dic[col] = "'01:01'"
        elif col.startswith("KIR"):
            if col[-1] in ["1", "2"]:
                default_dic[col] = "'001'"
            else:
                default_dic[col] = ""
        elif col == "MICA":
            default_dic[col] = "A001+A001"
        elif col == "MICB":
            default_dic[col] = "B001+B001"

    with open(fake_file_befunde, "w") as g:
        header = "sample_ID_int,sample_id_ext,client,{}\n".format(",".join(columns))
        g.write(header)

        for (sample_id_int, cell_line, mygene, target_allele, partner_allele) in data2:
            befunde = copy.copy(default_dic)
            # overwrite pretyping of target allele:
            if not mygene in gene_dic:
                gene_dic[mygene] = ["{}_1".format(mygene), "{}_2".format(mygene)]
                columns += gene_dic[mygene]
                for col in gene_dic[mygene]:
                    befunde[col] = "01:01"
            if len(gene_dic[mygene]) == 1:
                if mygene.startswith("MIC"):
                    if not partner_allele:
                        partner_allele = "{}001".format(mygene[-1])
                    befunde[mygene] = "{}+{}".format(target_allele, partner_allele).replace("MICA*", "A")
                else:
                    log.error("Cannot generate sensible fake pretyping: {} should have 2 columns!".format(mygene))
            else:
                for i, col in enumerate(gene_dic[mygene]):
                    if i == 0:
                        befunde[col] = target_allele.split("*")[1]
                    elif i == 1:
                        if partner_allele:
                            befunde[col] = partner_allele.split("*")[1]
                    else:
                        befunde[col] = ""
            # generate row for pretyping file:
            myline = "{},{},DKMS,".format(sample_id_int, cell_line)
            pretypings = [befunde[col] for col in columns]
            myline += ",".join(pretypings) + "\n"
            g.write(myline)

    log.debug("Fake Befunde file: {}".format(fake_file_befunde))
    return True, fake_file_befunde, ENA_id_map, ENA_gene_map


def get_pretypings_from_oracledb(project, local_cf, settings, log, parent=None):
    """writes pretyping file with the pretyping results from the oracle database
    """
    try:
        import db_external
    except Exception as E:
        log.exception(E)
        log.error("Could not import db_external, probably because cx_Oracle is missing on this computer!")
        return False, None, None

    success, alleles = find_alleles_per_project(project, log, parent)
    if not success:
        return False, None, None

    sample_ids = [[row[0], row[3]] for row in alleles]  # format [[sample_id_int, locus]]
    if not sample_ids:
        msg = "No samples found for project {}!".format(project)
        log.error(msg)
        if parent:
            QMessageBox.warning(parent, "Pretypings error", msg)
        return False, None, None

    pretypings, samples, problems = db_external.get_pretypings_from_limsrep(sample_ids, local_cf, log)
    output_file = os.path.join(settings["temp_dir"], "pretypings.csv")
    db_external.write_pretypings_file(pretypings, samples, output_file, log)
    return True, output_file, problems


def compare_2_files(file1, file2):
    """compares content of 2 text files, returns all differing lines as one string
    """
    with open(file1, "r") as f:
        text1 = f.read()
    with open(file2, "r") as f:
        text2 = f.read()

    identical = True
    comp_text = ""

    query = text2.split("\n")
    ref = text1.split("\n")
    for i in range(len(ref)):
        ref_line = ref[i]
        try:
            query_line = query[i]
        except IndexError:
            query_line = ""
        if query_line != ref_line:
            if identical:
                comp_text += "Changed line(s) found:\n"
                identical = False
            comp_text += f"[ln {i + 1}] Left:\t{ref_line}\n"
            comp_text += f"[ln {i + 1}] Right:\t{query_line}\n\n"

    if identical:
        comp_text += "Files are identical! :-)"

    return comp_text


def main(log):
    user = "admin"
    project = "20190515_ADMIN_KIR_1"
    sample_id_int = "ID17040887"

    settings = GUI_login.get_settings(user, log)

    # from typeloader_GUI import create_connection, close_connection

    # db_file = settings["db_file"]
    # mydb = create_connection(log, db_file)

    file1 = r"C:\Daten\local_data\TypeLoader\albrecht\projects\20201007_AL_HLA-DRB1_DR1\ID16313610\DKMS-LSL_ID16313610_DRB1_1.ena.txt"
    file2 = r"C:\Daten\local_data\TypeLoader\albrecht\projects\20201007_AL_HLA-DRB1_DR1\ID16313610\DKMS-LSL_ID16313610_DRB1_1.ena.txt"
    comp_text = compare_2_files(file1, file2)
    print(comp_text)
    # close_connection(log, mydb)


if __name__ == '__main__':
    log = general.start_log(level="DEBUG")
    log.info("<Start>")
    main(log)
    log.info("<End>")
