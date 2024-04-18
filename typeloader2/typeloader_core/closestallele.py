#!/usr/bin/env python

# This is the TypeLoader BLAST parser for finding the closest known allele for an input sequence

"""
This routine makes the following assumptions :
    1. The sequence being queried is a full length allele
    2. The closest known allele will align as a single HSP ( http://www.ncbi.nlm.nih.gov/BLAST/tutorial/#head2 )
    3. The GenDX generated file will have 2 sequences, only 1 of them is assumed to be not present in the database
"""

from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Align
import os
import re

try:
    from . import errors
except ImportError:
    import errors


###################################################


def get_closest_known_alleles(blast_xml_filename, target_family, settings, log):
    with open(blast_xml_filename) as xmlHandle:
        xmlParser = NCBIXML.parse(xmlHandle)

        # get the associated fasta file
        query_fasta_file = blast_xml_filename.replace("blast.xml", "fa")
        closestAllelesData = parse_blast(xmlParser, target_family, query_fasta_file, settings, log)
    return closestAllelesData


def print_hsp(hsp_query, hsp_subject, hsp_match, concat_HSPs, hsp_start, hsp_align_len, query_length):
    """little debugging function
    """
    print(f"hsp_query:\t\t {hsp_query[:50]}...{hsp_query[-50:]} (len: {len(hsp_query)})")
    print(f"hsp_subject:\t {hsp_subject[:50]}...{hsp_subject[-50:]} (len: {len(hsp_subject)})")
    print(f"hsp_match:\t\t {hsp_match[:50]}...{hsp_match[-50:]} (len: {len(hsp_match)})")
    print("concat_hsps: ", concat_HSPs)
    print("hsp_start: ", hsp_start)
    print("hsp_align_len: ", hsp_align_len)
    print("query_length: ", query_length)


def parse_blast(xml_records, target_family, query_fasta_file, settings, log):
    """
    Basic description of the XML below
    For more details on parsing the BLAST XML output:
     see http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:parsing-blast

    Heirarchy of XML output -

    Record (corresponding to each query) -> Alignment (individual hits) -> HSP
    """

    closestAlleles = {}
    hsp_start = 1
    for xmlRecord in xml_records:
        queryId = xmlRecord.query_id
        output_db = xmlRecord.database
        if output_db.endswith("parsedKIR.fa"):
            output_db = os.path.join(settings["dat_path"], settings["general_dir"], settings["reference_dir"],
                                     "parsedKIR.fa")
        elif output_db.endswith("parsedhla.fa"):
            output_db = os.path.join(settings["dat_path"], settings["general_dir"], settings["reference_dir"],
                                     "parsedhla.fa")
        else:
            log.error("Unknown reference file (in closestallele.py):", output_db)
        alignments = xmlRecord.alignments
        queryLength = xmlRecord.query_length
        if not alignments:
            log.error("No alignments found: probably not a supported locus!")
            raise ValueError(
                "No fitting result found in reference. Maybe this allele belongs to a gene not supported by TypeLoader?")

        potentialClosestAlleleAlignment = alignments[0]
        hsps = potentialClosestAlleleAlignment.hsps

        closestAlleleName = potentialClosestAlleleAlignment.hit_def
        if closestAlleleName.find(target_family) == -1:
            if closestAlleleName.startswith("MIC") and target_family == "HLA":
                pass
            else:
                closestAlleleName = potentialClosestAlleleAlignment.hit_id

        try:
            ref_sequence = SeqIO.to_dict(SeqIO.parse(output_db, "fasta"))[closestAlleleName].seq
        except KeyError:
            local_name = os.path.splitext(os.path.basename(query_fasta_file))[0]
            msg = f"Could not find {closestAlleleName} in current reference db!\n" \
                  f"This was originally assigned as closest allele to {local_name}. \n\n" \
                  f"Maybe the reference db version changed since ENA submission?\n\n" \
                  f"Please consult the user manual under 'Error: database changed between submissions' " \
                  f"for instructions how to proceed from here."
            raise KeyError(msg)
        query_sequence = SeqIO.to_dict(SeqIO.parse(query_fasta_file, "fasta"))[queryId].seq
        results = puzzle_HSPs_from_first_hit(hsps, ref_sequence, query_sequence, query_fasta_file)
        hsp_query, hsp_subject, hsp_match, concatHSPS, hsp_start, hsp_align_len = results
        query_start_overhang = 0

        #print_hsp(hsp_query, hsp_subject, hsp_match, concatHSPS, hsp_start, hsp_align_len, queryLength)
        if hsp_query == "" and hsp_subject == "" and hsp_match == "":
            closestAlleles[queryId] = None
        else:
            if hsp_align_len < queryLength:  # incomplete alignment:
                log.warning("Sequence did not align fully! Probably a mismatch within 3 bp of either sequence end!")
                results = fix_incomplete_alignment(ref_sequence, query_sequence, hsp_start, hsp_align_len, queryLength,
                                                   hsp_query, hsp_subject, hsp_match, closestAlleleName, log)
                (hsp_query, hsp_subject, hsp_match, hsp_align_len, hsp_start, query_start_overhang) = results
                log.debug("fix_incomplete_alignment finished")
            closestAlleles[queryId] = closest_allele_items(hsp_query, hsp_subject, hsp_match, closestAlleleName,
                                                           concatHSPS,
                                                           hsp_start, hsp_align_len, queryLength, query_start_overhang)

    if hsp_start != 1:
        log.warning("Incomplete sequence found: first {} bp missing!".format(hsp_start - 1))

    return closestAlleles


def make_global_alignment(ref_seq, query_seq, log):
    """
    performs global alignment of the two sequences via Biopython and returns list of alignment objects
    """
    log.info("\t=> performing global alignment...")
    aligner = Align.PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2
    aligner.query_right_open_gap_score = 0
    aligner.query_left_open_gap_score = 0
    aligner.query_right_extend_gap_score = 0
    aligner.query_left_extend_gap_score = 0
    aligner.mode = "global"

    alignments = aligner.align(ref_seq, query_seq)
    log.info("\t=> done")
    return alignments


def remove_end_gaps(ref, matched, query):
    """
    removes end gaps ('-'-characters on the right side) from an alignment incl. its sequences
    """
    match = re.match('-+', query[::-1])
    if match:  # end gaps found
        n = len(match.group())  # number of end gaps
        return ref[:-1 * n], matched[:-1 * n], query[:-1 * n]
    else:  # no end gaps found => return everything unchanged
        return ref, matched, query


def fix_incomplete_alignment(ref_seq, query_seq, hsp_start, hsp_align_len, query_length,
                             hsp_query, hsp_subject, hsp_match, closest_allele_name, log):
    alignments = make_global_alignment(ref_seq, query_seq, log)
    log.debug(f"=> found {len(alignments)} alignments")
    if not alignments:
        log.error("No alignment found, sorry! Aborting...")
        return hsp_query, hsp_subject, hsp_match, hsp_align_len, hsp_start, None
    if len(alignments) > 1:
        log.warning(f"Found {len(alignments)} possible alignments, choosing one at random (might not be the best one)!")

    a = alignments[0]

    log.info("Adding info from global alignment to local alignment...")
    # find alignment positions:
    aligned_query = a.aligned[0]
    q_start = aligned_query[0][0]
    query_start_overhang = a.aligned[1][0][0]  # relevant if query has longer 5' UTR than reference

    # fix alignments:

    missing_base_num_start_ref = hsp_start - 1 - q_start  # unaligned bases of the reference
    missing_base_num_start_query = query_seq.find(hsp_query[:50])  # unaligned bases of the query
    missing_base_num_end_ref = ref_seq[::-1].find(hsp_subject[-50:][::-1])
    missing_base_num_end_query = query_seq[-60:][::-1].find(hsp_query[-50:][::-1])  # unaligned bases of the query
    start_extended = 0
    end_extended = 0

    if missing_base_num_start_ref > 0:  # problem at sequence start
        if missing_base_num_start_query > 3 or missing_base_num_start_query < 0:  # different problem; probably #138
            msg = f"{missing_base_num_start_query} unaligned bases at alignment start "
            msg += f"when aligned to closest found allele {closest_allele_name}: "
            msg += "This sequence is probably too dissimilar to all known full-length alleles.\n"
            msg += "TypeLoader currently can't handle this allele, sorry!"
            log.warning(msg)
            raise errors.DevianceError(missing_base_num_start_query, closest_allele_name)

        log.debug(
            f"{missing_base_num_start_ref} bases missing at alignment start (= {missing_base_num_start_query} bases of the query)")
        missing_bases_ref = ref_seq[q_start: q_start + missing_base_num_start_ref]
        missing_bases_query = query_seq[:missing_base_num_start_query]

        if missing_base_num_start_query == missing_base_num_start_ref:  # no InDel
            # add bases to the alignment front in reverse order, to ensure the first base turns up first
            missing_bases_ref_inverted = missing_bases_ref[::-1]
            missing_bases_query_inverted = missing_bases_query[::-1]

            for i in range(missing_base_num_start_ref):
                ref_base = missing_bases_ref_inverted[i]
                query_base = missing_bases_query_inverted[i]
                hsp_query = query_base + hsp_query
                hsp_subject = ref_base + hsp_subject

                if ref_base == query_base:
                    hsp_match = "|" + hsp_match
                else:
                    hsp_match = " " + hsp_match
            hsp_align_len += missing_base_num_start_ref
            hsp_start = hsp_start - missing_base_num_start_ref

            start_extended = missing_base_num_start_ref

        else:  # InDel
            alignments = make_global_alignment(missing_bases_ref, missing_bases_query, log)
            a = alignments[0]
            # parse missing alignment part from string form of alignment:
            [ref_aligned, match_aligned, query_aligned] = str(a).strip().split("\n")
            hsp_query = query_aligned + hsp_query
            hsp_subject = ref_aligned + hsp_subject
            hsp_match = match_aligned.replace("-", " ").replace(".", " ") + hsp_match
            hsp_align_len += len(query_aligned)
            hsp_start = hsp_start - len(query_aligned)

            start_extended = len(query_aligned)

    if start_extended:
        log.debug(f"\t=> alignment start extended by {start_extended} bp")

    if missing_base_num_end_query:  # problem at sequence end
        log.debug(
            f"{missing_base_num_end_ref} bases missing at alignment end (= {missing_base_num_end_query} bases of the query)")
        if missing_base_num_end_ref < 4:
            missing_bases_ref = ref_seq[-1 * missing_base_num_end_ref:]
        else:  # long section missing due to incomplete allele, use next 3 bp for alignment
            missing_bases_ref = ref_seq[-1 * missing_base_num_end_ref:-1 * missing_base_num_end_ref + 3]
        missing_bases_query = query_seq[-1 * missing_base_num_end_query:]

        if len(missing_bases_ref) == len(missing_bases_query):  # no InDel
            for i in range(len(missing_bases_ref)):
                ref_base = missing_bases_ref[i]
                query_base = missing_bases_query[i]
                hsp_query += query_base
                hsp_subject += ref_base

                if ref_base == query_base:
                    hsp_match += "|"
                else:
                    hsp_match += " "
            hsp_align_len += missing_base_num_end_ref

            end_extended = missing_base_num_end_ref

        else:  # InDel
            alignments = make_global_alignment(missing_bases_ref, missing_bases_query, log)
            a = alignments[0]
            # parse missing alignment part from string form of mini-alignment:
            try:
                [ref_aligned, match_aligned, query_aligned] = str(a).strip().split("\n")
            except ValueError:
                log.error("Could not parse mini-alignment!")
                log.info("Leaving alignment un-extended")
                return hsp_query, hsp_subject, hsp_match, hsp_align_len, hsp_start, query_start_overhang

            ref_aligned, match_aligned, query_aligned = remove_end_gaps(ref_aligned, match_aligned,
                                                                        query_aligned)

            hsp_query += query_aligned
            hsp_subject += ref_aligned
            hsp_match += match_aligned.replace("-", " ").replace(".", " ")
            hsp_align_len += len(query_aligned)
            # ToDo: adjust UTR boundaries afterwards to get the correct number of missing bases

            end_extended = len(query_aligned)

    if end_extended:
        log.debug(f"\t=> alignment end extended by {end_extended} bp")

    if not start_extended and not end_extended:
        log.debug("\t=> did not extend alignment after all")

    return hsp_query, hsp_subject, hsp_match, hsp_align_len, hsp_start, query_start_overhang


def closest_allele_items(hsp_query, hsp_subject, hsp_match, closest_allele_name, concat_HSPs, hsp_start, hsp_align_len,
                         query_length, query_start_overhang):
    # "-" in the query means a deletion, "-" in the hit means an insertion, a gap in the alignment is a mismatch
    deletionPositions = [pos + 1 for pos in range(len(hsp_query)) if hsp_query[pos] == "-"]
    insertionPositions = [pos + 1 for pos in range(len(hsp_subject)) if hsp_subject[pos] == "-"]
    mismatchPositions = [pos + 1 for pos in range(len(hsp_match)) if ((hsp_match[pos] == " ")
                                                                      and ((pos + 1) not in deletionPositions) and (
                                                                              (pos + 1) not in insertionPositions))]
    # bases
    deletions = [hsp_subject[deletionPos - 1] for deletionPos in deletionPositions]
    insertions = [hsp_query[insertionPos - 1] for insertionPos in insertionPositions]

    if not (len(deletionPositions) or len(insertionPositions) or len(mismatchPositions)):
        exactMatch = True
    else:
        exactMatch = False

    mismatches = list(zip([hsp_query[mismatchPos - 1] for mismatchPos in mismatchPositions],
                          [hsp_subject[mismatchPos - 1] for mismatchPos in mismatchPositions]))

    # catch cases with undetected mismatches near end: (BLAST misses these)
    if hsp_align_len < query_length - query_start_overhang:  # if not whole of query sequence could be aligned
        mismatches.append(('?', '?'))
        mismatchPositions.append(hsp_align_len + 1)  # up to this point, the alignment worked
        exactMatch = False

    differences = {'deletionPositions': deletionPositions, 'insertionPositions': insertionPositions,
                   'mismatchPositions': mismatchPositions,
                   'mismatches': mismatches, 'deletions': deletions, 'insertions': insertions}

    closestAllele = {"name": closest_allele_name, "differences": differences,
                     "exactMatch": exactMatch, "concatHSPS": concat_HSPs,
                     "hitStart": hsp_start,
                     "alignLength": hsp_align_len,
                     "queryLength": query_length,
                     "queryStartOverhang": query_start_overhang}

    return closestAllele


def puzzle_HSPs_from_first_hit(hsps, ref_sequence, query_sequence, query_fasta_file):
    # TODO: (future) Validation of HSP puzzle, until that, take first HSP
    potentialClosestHSP = hsps[0]
    hsp_query = potentialClosestHSP.query
    hsp_subject = potentialClosestHSP.sbjct
    hsp_match = potentialClosestHSP.match
    hsp_start = potentialClosestHSP.sbjct_start
    hsp_align_len = potentialClosestHSP.align_length

    concatHSPS = False

    """

    1. Iterate through each hsp in the first hit
    2. possible cuts (X):

        a:              |========|
                             XXXX
                            |============|
            result:              |=======|

        b:                      |========|
                                 XXXX
                        |============|
            result:     |=======|

        c:              |============|          |============|
                                 XXXX            XXXX
                                |====================|
            result:                  |==========|

    3. every other possibility will be ignored

    """

    """

    hspMap = {}
    bool_sequence = [0] * len(ref_sequence)

    for hsp in hsps:

        overlap_beginning_start = 0
        overlap_beginning_end = 0
        overlap_beginning_end_marker = -1
        overlap_zeros_before = False
        overlap_break = 0
        hsp_break = False

        for run in range(hsp.sbjct_start - 1, hsp.sbjct_end):
            overlap_break += 1
            if bool_sequence[run] == 1:
                # it's an overlap

                # case a:
                if hsp.sbjct_start - 1 + overlap_beginning_start == run: overlap_beginning_start += 1

                # case b:
                if hsp.sbjct_start - 1 + overlap_break - 1 == run and overlap_zeros_before:
                    if overlap_beginning_end_marker == -1: overlap_beginning_end_marker = overlap_break - 1
                    if overlap_beginning_end_marker + overlap_beginning_end != run: break
                    overlap_beginning_end += 1
            else:
                # no overlap
                if overlap_beginning_start > 0 and overlap_beginning_end > 0:
                    # case c:
                    hsp_break = True
                    break
                bool_sequence[run] = 1
                overlap_zeros_before = True

        # break if:
        #   1: if query is a completely overlap
        #   2: query is not one of the cutting cases
        if overlap_beginning_start == hsp.sbjct_end - hsp.sbjct_start + 1 or hsp_break: break

        # next: add query to hspMap:
        save_sbjct_start = hsp.sbjct_start + overlap_beginning_start
        save_sbjct_end = hsp.sbjct_start + overlap_break - overlap_beginning_end - 1
        save_sbjct = hsp.sbjct[overlap_beginning_start:overlap_break - overlap_beginning_end]

        save_query_start = hsp.query_start + overlap_beginning_start
        save_query_end = hsp.query_start + overlap_break - overlap_beginning_end - 1
        save_query = hsp.query[overlap_beginning_start:overlap_break - overlap_beginning_end]

        save_mismatch = hsp.match[overlap_beginning_start:overlap_break - overlap_beginning_end]

        print "save_sbjct_start ", save_sbjct_start, "real_start ", hsp.sbjct_start
        print "save_sbjct_end ", save_sbjct_end, "real_end ", hsp.sbjct_end
        print "save_query_start", save_query_start, "real_start ", hsp.query_start
        print "save_query_end ", save_query_end, "real_end ", hsp.query_end
        print "overlap_beginning_start " , overlap_beginning_start
        print "overlap_beginning_end " , overlap_beginning_end
        print "overlap_break " , overlap_break
        print "hsp_break ", hsp_break
        #print "\n\n"

        hspMap[save_sbjct_start] = [save_sbjct_start, save_sbjct_end, save_sbjct, save_query_start, save_query_end, save_query, save_mismatch]

    hsp_query = ""
    hsp_subject = ""
    hsp_match = ""
    iter_step = 1

    concatHSPS = True if len(hspMap) > 1 else False
    if concatHSPS: EF.write_log(err_file, "WARNING", "Used more than 1 HSP in BLAST algorithm by using " + query_fasta_file)

    # return empty string, if subject not start at position 1
    # otherwise: wrong positions!
    if bool_sequence[0] == 0:
        return (hsp_query, hsp_subject, hsp_match, concatHSPS)
    """

    """
    hspMap_Keys:
        0 = subject start
        1 = subject end
        2 = subject nucleotides
        3 = query start
        4 = query end
        5 = query nucleotides
        6 = mismatch signs
    """

    """
    for key in sorted(hspMap.iterkeys()):

        if iter_step == 1:

            # first iteration
            hsp_subject = hspMap[key][2]
            hsp_query = hspMap[key][5]
            hsp_match = hspMap[key][6]

        else:

            #print "last: " , "subS: ", hspMap[last_key][0], "subE: " ,hspMap[last_key][1], "querS: ", hspMap[last_key][3], "querE: ", hspMap[last_key][4]
            #print "now: " , "subS: ", hspMap[key][0], "subE: ", hspMap[key][1], "querS: ", hspMap[key][3], "querE: ", hspMap[key][4]

            if hspMap[last_key][4] >= hspMap[key][3]:

                #print "overlap deletion"
                hsp_starting_point = 1
                diff = hspMap[last_key][4]-hspMap[key][3]
                deleted_nucleotides = ref_sequence[hspMap[last_key][1]-diff:hspMap[key][0]]
                inserted_nucleotides = ""

            elif hspMap[last_key][1] >= hspMap[key][0]:

                #print "overlap insertion"
                hsp_starting_point = 1
                diff = hspMap[last_key][1] - hspMap[key][0]
                inserted_nucleotides = query_sequence[hspMap[last_key][4]:hspMap[key][3]]
                deleted_nucleotides = ""

            else:

                #print "no overlapping"
                hsp_starting_point = 0
                diff = 0
                deleted_nucleotides = ref_sequence[hspMap[last_key][1]:hspMap[key][0] - 1]
                inserted_nucleotides = query_sequence[hspMap[last_key][4]:hspMap[key][3] - 1]

            #print diff
            #print hsp_query
            hsp_subject = hsp_subject[0:len(hsp_subject)-diff]
            hsp_query = hsp_query[0:len(hsp_query)-diff]
            hsp_match = hsp_match[0:len(hsp_match)-diff]
            #print hsp_query

            new_subject = hspMap[key][2][hsp_starting_point:len(hspMap[key][2])]
            new_query = hspMap[key][5][hsp_starting_point:len(hspMap[key][5])]
            new_match = hspMap[key][6][hsp_starting_point:len(hspMap[key][6])]

            number_deletions = len(deleted_nucleotides)
            number_insertions = len(inserted_nucleotides)

            #print "number_deletions: ", number_deletions
            #print "deleted: ", deleted_nucleotides
            #print "number_insertions: ", number_insertions
            #print "insertions: ", inserted_nucleotides

            # deletion --> see query
            #hsp_query = hsp_query + "\n" + inserted_nucleotides + "\n" + "".join(["-"]*number_deletions) + "\n" + new_query
            # insertions --> see subject
            #hsp_subject = hsp_subject + "\n" + "".join(["-"]*number_insertions) + "\n" + deleted_nucleotides + "\n" + new_subject
            # mismatches --> see match
            #hsp_match = hsp_match + "\n" + "".join(["X"]*number_insertions) + "\n" + "".join(["X"]*number_deletions) + "\n" + new_match

            # deletion --> see query
            hsp_query = hsp_query + inserted_nucleotides + "".join(["-"]*number_deletions) + new_query
            # insertions --> see subject
            hsp_subject = hsp_subject + "".join(["-"]*number_insertions) + deleted_nucleotides + new_subject
            # mismatches --> see match
            hsp_match = hsp_match + "".join(["X"]*number_insertions) + "".join(["X"]*number_deletions) + new_match

        iter_step += 1
        last_key = key
    """

    """
    curr_file = str.split(ntpath.basename(query_fasta_file), ".")[0]
    filename = "/home/markus/playground/report/Validation_ClosestAllele/output/" + curr_file + ".query.txt"
    filebuff = open(filename,'w')
    filebuff.write(str(hsp_query))
    filebuff.close()

    filename = "/home/markus/playground/report/Validation_ClosestAllele/output/" + curr_file + ".sbjct.txt"
    filebuff = open(filename,'w')
    filebuff.write(str(hsp_subject))
    filebuff.close()

    filename = "/home/markus/playground/report/Validation_ClosestAllele/output/" + curr_file + ".mismatch.txt"
    filebuff = open(filename,'w')
    filebuff.write(str(hsp_match))
    filebuff.close()
    """

    return hsp_query, hsp_subject, hsp_match, concatHSPS, hsp_start, hsp_align_len


if __name__ == "__main__":
    pass
