#!/usr/bin/env python

import re
import os
from subprocess import run, PIPE
from collections import defaultdict
from .EMBLfunctions import fasta_generator
from .xmlfuncs import *

"""
The BLAST db has to be formatted like so:
~/blast/ncbi-blast-2.2.31+/bin/makeblastdb -dbtype nucl -in parsedhla.fa -parse_seqids -out parsedhla

The FASTA file corresponding to the GenDX file and the BLAST output file will all be written to the same location
as the GenDX XML file
"""


##########################################################

def getAlleleSequences(xmlFile, log):
    alleles = {}
    xmlData = parseXML(xmlFile)
    alleleNames = getAlleleNames(xmlData)
    data_dic = get_additional_XML_info(xmlData, log)

    for alleleName in alleleNames:
        haplotype_list = getHaplotypeIds(xmlData, alleleName)
        alleles[alleleName] = sequenceFromHaplotype(xmlData, haplotype_list)

    return alleles, data_dic


def blastSequences(inputFastaFile, parsedFasta, settings, log,
                   blastOutputFormat="5"):  # 5 corresponds to XML BLAST output
    blast = settings["blast_path"]
    database = parsedFasta
    blastXmlOutputFile = inputFastaFile.replace(".fasta", ".blast.xml").replace(".fa", ".blast.xml")
    blast_command = [blast,
                     "-query", inputFastaFile,
                     "-parse_deflines",
                     "-db", database,
                     "-dust", "no",
                     "-soft_masking", "false",
                     "-outfmt", blastOutputFormat,
                     "-out", blastXmlOutputFile]
    log.debug("Blast command:")
    log.debug(" ".join(blast_command))

    result = run(blast_command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    if result.stdout:
        log.info(result.stdout)
    if result.returncode != 0:
        log.error(result.stderr)
        log.error(f"Blast did not succeed; return code: {result.returncode}")

    if not os.path.isfile(blastXmlOutputFile):
        log.error("BlastXMLFile not generated!")
        return False
    return blastXmlOutputFile


def parse_fasta_header(fasta_header):
    """parses header of a fastq file
    """
    s = fasta_header.split(";")
    header_data = defaultdict(lambda: "")
    if len(s) > 1:  # fasta header generated by DR2S
        s2 = s[0].split()
        seq_name = s2[0]
        items = s2[1:] + s[1:]
        for item in items:
            data = item.replace('"', "").split("=")
            try:
                if data[1] in ["list()", "list(NULL)"]:
                    data[1] = ""
                header_data[data[0]] = data[1]
            except IndexError:
                print("Error: cannot distinguish key:value pair from '{}' in fasta_header".format(item))
    else:  # manual or GenXD fasta header
        seq_name = fasta_header

    return seq_name, header_data


def sanity_check_seq(seq, log):
    """checks for non-ATGC-characters in seq
    """
    log.debug("Checking sequence for non-ATGC bases...")
    ok_char = defaultdict(lambda: False)
    for char in ["A", "T", "G", "C"]:
        ok_char[char] = True

    problems = []
    for (i, char) in enumerate(seq.upper()):
        if not ok_char[char]:
            log.warning("Non-ATGC-character found: {} in position {}".format(char, i))
            problems.append((char, i))

    if problems:
        ok = False
        msg = "The uploaded allele contains the following {} non-ATGC base(s):\n".format(len(problems))
        for (char, i) in problems:
            msg += "- position {}: {}\n".format(i + 1, char)
        msg += "\nPlease fix this in the raw file and try again!"
        log.info("File rejected for non-ATGC sequence.")
    else:
        msg = "all bases are ATGC"
        ok = True

    return ok, msg


def blast_raw_seqs(input_filename, filetype, settings, log, use_given_reference=False):
    """parses raw allele file (fasta or XML)
    """
    if filetype == "XML":
        log.debug("\tConverting xml to fasta...")
        alleles, xml_data_dic = getAlleleSequences(input_filename, log)
        fastaFilename = input_filename.replace(".xml", ".fa")
        with open(fastaFilename, "w") as fastaFile:
            for alleleName in list(alleles.keys()):
                fastaFile.write(">%s\n" % alleleName)
                fastaFile.write("%s\n" % alleles[alleleName])
    else:
        fastaFilename = input_filename
        xml_data_dic = {}

    kir = settings["gene_kir"]
    hla = settings["gene_hla"]

    log.debug("\tReading fasta for sanity check...")
    header = ""
    for myfasta in fasta_generator(fastaFilename):
        (header, seq) = myfasta
        (ok, msg) = sanity_check_seq(seq, log)
        if not ok:
            return False, "Non-ATGC-Error", msg

    log.debug("\tParsing fasta header...")
    seq_name, header_data = parse_fasta_header(header)
    locus = header_data["locus"]
    if locus:  # if DRS2 fasta file
        if locus.startswith("KIR"):
            targetFamily = kir
        else:
            targetFamily = hla
    else:
        if re.search(kir, seq_name):
            targetFamily = kir
        else:
            targetFamily = hla

    if use_given_reference:
        reference_path = use_given_reference
    else:
        reference_path = os.path.join(settings["dat_path"], settings["general_dir"],
                                      settings["reference_dir"])
    parsed_kir = settings["parsed_kir"]
    parsed_hla = settings["parsed_hla"]
    hla_dat = settings["hla_dat"]
    kir_dat = settings["kir_dat"]
    hla_version = settings["hla_version"]
    kir_version = settings["kir_version"]

    if targetFamily == kir:
        parsedFasta = os.path.join(reference_path, parsed_kir)
        allelesFilename = os.path.join(reference_path, kir_dat)
        versionFilename = os.path.join(reference_path, kir_version)
    else:
        parsedFasta = os.path.join(reference_path, parsed_hla)
        allelesFilename = os.path.join(reference_path, hla_dat)
        versionFilename = os.path.join(reference_path, hla_version)

    log.debug("\tBlasting sequence...")
    try:
        BlastXMLFile = blastSequences(fastaFilename, parsedFasta, settings, log)
    except Exception as E:
        log.exception(E)
        return False, "Error while trying to BLAST raw sequence", repr("E")
    if not BlastXMLFile:
        msg = "BlastXMLFile not generated!\n"
        msg += "Please make sure the BLAST path in your user settings is correct!\n"
        msg += "(Current BLAST path: {})".format(settings["blast_path"])
        return False, "Error while trying to BLAST raw sequence", msg

    try:
        with open(versionFilename) as f:
            curr_version = f.read().strip()
            header_data["ref_version"] = curr_version
    except Exception as E:
        log.error(E)
        msg = f"Could not read the reference file:\n\n{str(E)}\n\n"
        msg += "Did you maybe specify a wrong gene family?"
        return False, "Error while opening reference file", msg
    return BlastXMLFile, targetFamily, fastaFilename, allelesFilename, header_data, xml_data_dic


if __name__ == '__main__':
    xmlFile = r"T:\nobackup\typeloader_temp\xmlfiles\0ID14893565.xml"
    getAlleleSequences(xmlFile, None)
#     header = 'DEDKM9497312_A haplotype="A";locus="KIR3DL3";ref="KIR3DL3*0140201";LIMS_DONOR_ID="ID121212121";Spendernummer="blaaatest";second="20:10:01:01";third="";note="bbbrnurneanuiarennrnrn_urnertue_neriunedute";short_read_data="yes";short_read_type="illumina";long_read_data="yes";long_read_type="pacbio";software="DR2S";version="0.0.4";date="2018-05-30"'
#     seq_name, header_data = parse_fasta_header(header)
#     print(seq_name)
#     for key in header_data:
#         print (key, "\t", header_data[key])
#     print(blastSequences(argv[1]))
# print blastGenDXSeqs(argv[1])
