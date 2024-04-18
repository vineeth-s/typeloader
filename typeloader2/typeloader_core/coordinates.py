#!/usr/bin/env python

# This routine establishes the co-ordinates for the new allele from the closest known allele
# The closest known allele is determined using getclosestKnownAllele from closestallele.py


import re
import os
from Bio import SeqIO
from collections import defaultdict
from .closestallele import get_closest_known_alleles
from .hla_embl_parser import read_dat_file
from .imgtTransform import changeToImgtCoords
from .errors import MissingUTRError, IncompleteSequenceWarning

from copy import copy
from math import ceil, floor
from pickle import load
from sys import argv

def getSpecificCodonChords(closestallele, exon1Length):

    if any(re.findall(r'KIR2DL1|KIR2DL2|KIR2DL3|KIR2DL5|KIR3DP1|KIR2DS|KIR3DL|KIR2DP|KIR3DP|', closestallele, re.IGNORECASE)):
        numExon1Coords = 21
    elif "KIR2DL4" in closestallele:
        numExon1Coords = 23
    else:
        numExon1Coords = int(floor(float(exon1Length)/3))

    return numExon1Coords


def getMismatchData(annotations):

    cdsMap = annotations["cdsMap"]
    query_start_overhang = annotations["queryStartOverhang"]
    alleleSeq = annotations["sequence"][query_start_overhang:]
    closestallele = annotations["closestAllele"]
    missing_bp = annotations["missing_bp"]

    # changed old statement to generate closestAlleleCdsSeq,
    # because pseudoexons could be in KIR
    closestAlleleCdsSeq = annotations["closestAlleleCdsSequence"]
    
    genomicExonCoords = list(cdsMap.keys())
    genomicExonCoords.sort()

    cdsExonCoords = [cdsMap[genomicExonCoord] for genomicExonCoord in genomicExonCoords]

    cdsSeq = "".join([alleleSeq[(start-1):end] for (start,end) in genomicExonCoords])

    mmCodons = []
        
    exon1Length = cdsExonCoords[0][1]
    numExon1Coords = getSpecificCodonChords(closestallele, exon1Length)

    for mmIndex in range(len(annotations["differences"]["mismatches"])):
        genomicPos = annotations["differences"]["mismatchPositions"][mmIndex]
        cdsPos = annotations["imgtDifferences_orig"]["mismatchPositions"][mmIndex][0]
        if genomicPos == cdsPos + missing_bp: # non-CDS positions
            mmCodons.append(())
            continue

        # CDS positions:
        canonicalMMCodonNum = int(ceil(float(cdsPos)/3))
        
        """
        # older version which got class I codon assignments correct, and class II assignments wrong
        # IMGT assigns the first codon in exon 2 as codon 1, and the codons in exon 1 with a -1 backwards
        if canonicalMMCodonNum > numExon1Coords: imgtMMCodonNum = canonicalMMCodonNum - numExon1Coords
        else: imgtMMCodonNum = canonicalMMCodonNum - (numExon1Coords + 1)
        """

        """
        The first codon per IMGT definition corresponds to the start of the mature protein
        For class I, this is exon 2, and for class II, this is exon 1
        """
        if "HLA-D" in closestallele: # class II 
            imgtMMCodonNum = canonicalMMCodonNum
        else:
            if canonicalMMCodonNum > numExon1Coords: 
                imgtMMCodonNum = canonicalMMCodonNum - numExon1Coords
            else: 
                imgtMMCodonNum = canonicalMMCodonNum - (numExon1Coords + 1)

        """
        Getting rid of the codon sequence madness, sort of bruting it

        #for pos in range(cdsPos - 3,cdsPos + 3):
        #    if pos == cdsPos: isCdsPos = "*"
        #    else: isCdsPos = ""
        #    #print "%s\t%s\t%s\t%s\t%s" % (pos, closestAlleleCdsSeq[pos], cdsSeq[pos], pos + 1, isCdsPos)
        ##print "-" * 10

        #cdsPos -= 1
        if not (cdsPos % 3):# coords stored are 1-based, python indices are 0-based
            #cdsPos += 1
            actualCodonSeq, closestCodonSeq = cdsSeq[cdsPos - 1  : cdsPos + 2 ], closestAlleleCdsSeq[cdsPos - 1 : cdsPos + 2] #
        elif (cdsPos % 3) == 2:
            #cdsPos += 1
            actualCodonSeq, closestCodonSeq = cdsSeq[cdsPos : cdsPos + 3], closestAlleleCdsSeq[cdsPos : cdsPos + 3]
        else:
            #cdsPos += 1
            actualCodonSeq, closestCodonSeq = cdsSeq[cdsPos - 2 : cdsPos + 1], closestAlleleCdsSeq[cdsPos - 2  : cdsPos + 1]
        """

        cdsCodonHash, closestAlleleCodonHash = {},{}
        startPos = 0

        for index in range(1,(len(cdsSeq)//3) + 1):
            cdsCodonHash[index] = cdsSeq[startPos:startPos+3]
            closestAlleleCodonHash[index] = closestAlleleCdsSeq[startPos:startPos+3]
            startPos += 3
        len_dif = len(closestAlleleCdsSeq) - len(cdsSeq)
        
        needs_shifting = defaultdict(lambda: 0)
        if len_dif > 2: # target allele contains deletion => cdsSeq does not align with complete reference seq
            startPos = len(cdsSeq) + 3
            for index in range(len(cdsSeq)//3 + 1, len(closestAlleleCdsSeq)//3 + 1):
                cdsCodonHash[index] = cdsSeq[-3:]
                mycodon = closestAlleleCdsSeq[startPos:startPos+3]
                closestAlleleCodonHash[index] = mycodon
                startPos += 3
                needs_shifting[index] = 1
        try:
            if cdsCodonHash[canonicalMMCodonNum] != closestAlleleCodonHash[canonicalMMCodonNum]:
                mmCodons.append((imgtMMCodonNum,(cdsCodonHash[canonicalMMCodonNum], closestAlleleCodonHash[canonicalMMCodonNum - needs_shifting[canonicalMMCodonNum]])))
            elif cdsCodonHash[canonicalMMCodonNum - 1] != closestAlleleCodonHash[canonicalMMCodonNum - 1]:
                mmCodons.append((imgtMMCodonNum,  (cdsCodonHash[canonicalMMCodonNum - 1], closestAlleleCodonHash[canonicalMMCodonNum - 1])))
            else:
                mmCodons.append((imgtMMCodonNum,  (cdsCodonHash[canonicalMMCodonNum + 1], closestAlleleCodonHash[canonicalMMCodonNum + 1])))
        except KeyError:
            raise KeyError("Could not calculate annotation positions: Reference allele {} has no codon number {}!".format(closestallele, canonicalMMCodonNum))

    return mmCodons


def pprint_annotations(annotations):
    """Print the relevant annotations in a readable way.

    This is just a helper function for debugging.
    """
    from pprint import pprint

    for allele in annotations:
        myannotations = copy(annotations[allele])

        # remove uninteresting keys:
        for key in ["closestAllele", "closestAlleleCdsSequence", "concatHSPS", "extraInformation", "features", "isExactMatch", "sequence"]:
            try:
                myannotations.pop(key)
            except KeyError:
                pass
        myannotations["closestAlleleSequence"] = myannotations["closestAlleleSequence"][:20] + "..."

        pprint(myannotations)


def getCoordinates(blastXmlFilename, allelesFilename, targetFamily, settings, log, isENA=True,
                   incomplete_ok=False):
    if "restricted_db" in allelesFilename:
        allelesFilename = os.path.join(settings["root_path"], settings["general_dir"],
                                       settings["reference_dir"],
                                       os.path.basename(allelesFilename))
    allAlleles, _ = read_dat_file(allelesFilename, targetFamily, isENA)
    closestAlleles = get_closest_known_alleles(blastXmlFilename, targetFamily, settings, log)
    seqsFile = blastXmlFilename.replace(".blast.xml", ".fa")

    try: 
        seqsHandle = open(seqsFile)
    except IOError:
        seqsFile = blastXmlFilename.replace(".blast.xml",".fasta")
        seqsHandle = open(seqsFile)

    seqsHash = SeqIO.to_dict(SeqIO.parse(seqsHandle, "fasta"))
    annotations = processAlleles(closestAlleles, allAlleles, seqsHash, incomplete_ok)
    # for cell_line in annotations:
    #     for key in annotations[cell_line]:
    #         item = annotations[cell_line][key]
    #         if isinstance(item, dict):
    #             print("'{}' : {}".format(key, "{"))
    #             for key2 in item:
    #                 print("\t'{}' : {}".format(key2, item[key2]))
    #             print("\t}")
    #         else:
    #             if key == 'closestAlleleSequence':
    #                 print("'{}' : {}...".format(key, annotations[cell_line][key][:20]))
    #             else:
    #                 print("'{}' : {}".format(key, annotations[cell_line][key]))

    for gendxAllele in list(annotations.keys()):
        if not annotations[gendxAllele]:
            continue
        alleleSeq = str(seqsHash[gendxAllele].seq)
        annotations[gendxAllele]["sequence"] = alleleSeq
        annotations[gendxAllele]["imgtDifferences"]["mmCodons"] = getMismatchData(annotations[gendxAllele])
        
    seqsHandle.close()

    return annotations


def shift_coordinates_for_missing_bp(missing_bp, coordinates):
    """if part of UTR5 is missing, shift coordinates in 5' direction 
    """
    #check for missing UTRs:
    UTR5_start_old = coordinates[0][1]
    UTR5_start_new = UTR5_start_old - missing_bp
    if UTR5_start_new <= 0:
        raise MissingUTRError(5)
    
    # shift coordinates:
    coordinates_new = []
    for i, (start, end) in enumerate(coordinates):
        if i == 0: # UTR5
            new = (1, UTR5_start_new)
        elif i == len(coordinates) - 1: # UTR3
            new = (start - missing_bp, end)
        else:
            new = (start - missing_bp, end - missing_bp)
        coordinates_new.append(new)

    return coordinates_new
    

def check_incomplete_utrs(missing_bp, utr3Length_orig, utr3_length_new, incomplete_ok):
    """checks whether either UTR is incomplete and issues a warning if that is the case
    (unless incomplete_ok is already set to True)
    """ 
    missing_bp_end = utr3Length_orig - utr3_length_new
    if not incomplete_ok:
        if missing_bp > 0 or missing_bp_end > 0:
            raise IncompleteSequenceWarning(missing_bp, missing_bp_end)
    return missing_bp_end
            
            
def shift_differences_for_missing_bp(missing_bp, differences, imgtDifferences, closestAlleleSequence):
    """if part of UTR5 is missing, shift differences in 3' direction 
    so they mark positions in the alignment, no the target sequence 
    """
    len_seq = len(closestAlleleSequence)
    # shift differences:
    new_differences = {}
    for key in differences:
        data = differences[key]
        new_differences[key] = []
        if data:
            if key.endswith("Positions"):
                for pos in data:
                    new_differences[key].append(pos + missing_bp)
            else:
                new_differences[key] = differences[key]
                
    # shift imgtDifferences:
    new_imgtDifferences = {}
    for key in imgtDifferences:
        data = imgtDifferences[key]
        new_imgtDifferences[key] = []
        if data:
            if key.endswith("Positions"):
                for (pos, value) in data:
                    if value: # CDS change => position should not be shifted!
                        new_imgtDifferences[key].append((pos, value))
                    else: # non-CDS
                        new_imgtDifferences[key].append((pos + missing_bp, value))
            else:
                new_imgtDifferences[key] = imgtDifferences[key]

    # remove changes behind last base:
    for (key_pos_list, key_change_list) in [("mismatchPositions", "mismatches"),
                                            ("deletionPositions", "deletions"),
                                            ("insertionPositions", "insertions")]:
        pos_list1 = new_differences[key_pos_list]
        change_list1 = new_differences[key_change_list]
        pos_list2 = new_imgtDifferences[key_pos_list]
        change_list2 = new_imgtDifferences[key_change_list]

        list_len = min([len(pos_list1), len(pos_list2), len(change_list1), len(change_list2)])  # these should be equal
        num_removeme = 0  # remove this many elements from end of each list (out of scope of reference sequence)
        for i in reversed(range(list_len)):
            pos = pos_list1[i]
            if pos > len_seq:
                num_removeme += 1

        if num_removeme:
            mylists = [pos_list1, pos_list2, change_list1]
            if change_list2 is not change_list1:  # is identical object, at least in some cases
                mylists.append(change_list2)
            for mylist in mylists:
                for i in range(num_removeme):
                    del(mylist[-1])

    return new_differences, new_imgtDifferences


def processAlleles(closestAlleles, allAlleles, hashOfQuerySequences, incomplete_ok = False):

    annotations = {}
    for alleleQuery in list(closestAlleles.keys()):
        if not closestAlleles[alleleQuery]:
            annotations[alleleQuery] = None
            continue
        closestAlleleName = closestAlleles[alleleQuery]["name"]
        differences = closestAlleles[alleleQuery]["differences"]
        isExactMatch = closestAlleles[alleleQuery]["exactMatch"]
        concatHSPS = closestAlleles[alleleQuery]["concatHSPS"]
        hitStart = closestAlleles[alleleQuery]["hitStart"]
        missing_bp = hitStart - 1
        alignLength = closestAlleles[alleleQuery]["alignLength"]
        queryLength = closestAlleles[alleleQuery]["queryLength"]
        query_start_overhang = closestAlleles[alleleQuery]["queryStartOverhang"]

        features, coordinates, extraInformation, closestAlleleCdsSequence, \
            closestAlleleSequence = calculateCoordinates(closestAlleleName, allAlleles, differences,
                                                        len(hashOfQuerySequences[alleleQuery]), missing_bp)
        
        coordinates = shift_coordinates_for_missing_bp(missing_bp, coordinates)
        
        (UTR3start, UTR3end) = coordinates[-1]
        UTR3length = UTR3end - UTR3start + 1

        if UTR3length <= 0:
            raise MissingUTRError(3)

        UTR3length_orig = len(allAlleles[closestAlleleName].UTR3)
        missing_bp_end = check_incomplete_utrs(missing_bp, UTR3length_orig, UTR3length, incomplete_ok)

        coordinates, imgtDifferences_orig, cdsMap = changeToImgtCoords(features, coordinates, differences)
        
        differences, imgtDifferences = shift_differences_for_missing_bp(missing_bp, differences, imgtDifferences_orig,
                                                                        closestAlleleSequence)
        
        annotations[alleleQuery] = {"features": features, "coordinates": coordinates,
                                    "imgtDifferences": imgtDifferences, "differences": differences,
                                    "closestAllele": closestAlleleName, "cdsMap": cdsMap,
                                    "closestAlleleCdsSequence": closestAlleleCdsSequence,
                                    "closestAlleleSequence": closestAlleleSequence, "isExactMatch": isExactMatch,
                                    "extraInformation": extraInformation, "concatHSPS": concatHSPS,
                                    "missing_bp": missing_bp, "missing_bp_end": missing_bp_end,
                                    "imgtDifferences_orig": imgtDifferences_orig,
                                    "queryLength": queryLength, "alignLength": alignLength,
                                    "queryStartOverhang": query_start_overhang}
        
        
    return annotations

def getClosestAlleleCoordinates(alleleData, queryLength):

    features = []
    coordinates = []
    exonIds = []
    extraInformation = {"pseudoexon":{}, "exon_number":{}, "intron_number":{}}
    closestAlleleCdsSequence = ""
    closestAlleleSequence = ""

    if "utr5" in alleleData.utrpos_dic:
        features.append("utr5")
        coordinates.append((alleleData.utrpos_dic["utr5"][0] + 1, alleleData.utrpos_dic["utr5"][1]))

    for exon, pseudoexon in alleleData.pseudo_exon_dic.items():
        extraInformation["pseudoexon"][exon] = pseudoexon
    for exon, exon_number in alleleData.exon_num_dic.items():
        extraInformation["exon_number"][exon] = exon_number
    for intron, intron_number in alleleData.intron_num_dic.items():
        extraInformation["intron_number"][intron] = intron_number

    #exonIds = [(exonNumber,"e") for exonNumber in alleleData.exonpos_dic.keys()]
    for exonNumber in list(alleleData.exonpos_dic.keys()):
        if extraInformation["pseudoexon"][exonNumber]:
            exonIds.append((exonNumber, "epseudo"))
        else:
            exonIds.append((exonNumber, "e"))

    intronIds = [(intronNumber, "i") for intronNumber in list(alleleData.intronpos_dic.keys())]

    cdsIds = copy(exonIds)
    cdsIds.extend(intronIds)
    cdsIds.sort()
    features.extend(cdsIds)

    for cdsId in cdsIds:
        if cdsId[1] == "e":
            coordinates.append((alleleData.exonpos_dic[cdsId[0]][0] + 1, alleleData.exonpos_dic[cdsId[0]][1]))
            closestAlleleCdsSequence += alleleData.seq[alleleData.exonpos_dic[cdsId[0]][0]:alleleData.exonpos_dic[cdsId[0]][1]]
        elif cdsId[1] == "epseudo":
            coordinates.append((alleleData.exonpos_dic[cdsId[0]][0] + 1, alleleData.exonpos_dic[cdsId[0]][1]))
        else: coordinates.append((alleleData.intronpos_dic[cdsId[0]][0] + 1, alleleData.intronpos_dic[cdsId[0]][1]))

    if "utr3" in alleleData.utrpos_dic:
        features.append("utr3")
        coordinates.append((alleleData.utrpos_dic["utr3"][0] + 1, queryLength))

    closestAlleleSequence = alleleData.seq

    return (features, coordinates, extraInformation, closestAlleleCdsSequence, closestAlleleSequence)

def calculateCoordinates(alleleName, alleles, differences, queryLength, missing_bp):

    allele = alleles[alleleName]
    features, coordinates, extraInformation, closestAlleleCdsSequence, \
        closestAlleleSequence = getClosestAlleleCoordinates(allele, queryLength)

    #features = ['utr5', (1, 'e'), (1, 'i'), (2, 'e'), (2, 'i'), ... , 'utr3']
    #coordinates = [(1, 267), ... , (10737, 10561)]

    insertions, deletions, mismatches = differences["insertionPositions"], differences["deletionPositions"], \
                                        differences["mismatchPositions"]
    
    # shift positions to alignment if missing_bp:
    myinsertions = [pos + missing_bp for pos in insertions]
    mydeletions = [pos + missing_bp for pos in deletions]
    mymismatches = [pos + missing_bp for pos in mismatches]

    """
    The logic for assigning the co-ordinates for the new allele is quite complex. The steps are listed below :

    1. Iterate through the gene model features
    2. For each feature, check the number of insertions and deletions that lie within that feature
    3. (numInsertions - numDeletions) will be the change in the 3' end of that feature, even if an insertion or deletion is at the boundary
    4. If there are mismatches within that region, these are checked individually if they are before or after an insertion or deletion and changes made
    5. Every feature following that feature has its boundaries changed by (numInsertions - numDeletions)
    6. Every mismatch in every following feature is also changed by (numInsertions - numDeletions)

    Note : We take advantage of the fact that a list is mutable, but have to keep in mind tuples are not
    """

    #Step 1
    for coordIndex in range(len(coordinates)):
        #Step 2
        regionBegin, regionEnd = coordinates[coordIndex]
        regionInsertions = [insertion for insertion in myinsertions if ((insertion >= regionBegin) and (insertion <= regionEnd))]
        regionDeletions = [deletion for deletion in mydeletions if ((deletion >= regionBegin) and (deletion <= regionEnd))]


        coordChange = len(regionInsertions) - len(regionDeletions)
        if not coordChange: continue

        # Step 3
        if coordIndex != len(coordinates) - 1:
            coordinates[coordIndex] = (coordinates[coordIndex][0], coordinates[coordIndex][1]+coordChange)
        else:
            # no position pushing in the last iteration has to be done
            # start of 3' UTR was pushed already in step5, if there are insertions/deletions
            # end of 3' UTR = sequence length
            coordinates[coordIndex] = (coordinates[coordIndex][0], coordinates[coordIndex][1])

        #Step 4
        for mismatchIndex in range(len(mymismatches)):
            mismatch = mymismatches[mismatchIndex]
            if (mismatch >= regionBegin) and (mismatch <= regionEnd):
                for regionInsertion in regionInsertions:
                    if mismatch > regionInsertion: mymismatches[mismatchIndex] += 1
                for regionDeletion in regionDeletions:
                    if mismatch > regionDeletion: mymismatches[mismatchIndex] -= 1
            #Step 6
            if mismatch > regionEnd:
                mymismatches[mismatchIndex] += coordChange

        #Step 5
        for nextFeatureIndex in range(coordIndex + 1, len(coordinates)):
            translated_region = (coordinates[nextFeatureIndex][0] + coordChange, coordinates[nextFeatureIndex][1] + coordChange)
            utr_end = coordinates[nextFeatureIndex] = (coordinates[nextFeatureIndex][0] + coordChange, coordinates[nextFeatureIndex][1])
            # push only the start of 3' UTR, if there are insertions/deletions
            # no entry if coordIndex + 1 = len(coordinates), see Step 3
            coordinates[nextFeatureIndex] = translated_region if nextFeatureIndex != len(coordinates) - 1 else utr_end

    return (features, coordinates, extraInformation, closestAlleleCdsSequence, closestAlleleSequence)

if __name__ == '__main__':

    t = getCoordinates(argv[1], argv[2], argv[3], isENA=True)
    #print t
    """
    from getAlleleSeqsAndBlast import blastSequences
    from sys import argv

    blastXmlFile = blastSequences(fastaFile)
    typeloaderOutput = getCoordinates(argv[1])

    typeloaderFile = argv[1].replace(".blast.xml",".typeloader.txt")
    typeloaderHandle = open(typeloaderFile,"w")

    for gendxAllele in typeloaderOutput.keys():
        typeloaderHandle.write("TypeLoader Query : %s\n" % gendxAllele)

        typeloaderData = typeloaderOutput[gendxAllele]
        typeloaderHandle.write("Closest Allele : %s\n\n" % typeloaderData["closestAllele"])

        for feature, coordinates in zip(typeloaderData["features"],typeloaderData["coordinates"]):
            if feature == "utr5":
                typeloaderHandle.write("UTR5\t\t")
            elif feature == "utr3":
                typeloaderHandle.write("UTR3\t\t")
            else:
                featureNum, featureCode = feature
                if featureCode == "e":
                    featureName = "Exon"
                    typeloaderHandle.write("%s %s\t\t" % (featureName, featureNum))
                else:
                    featureName = "Intron"
                    typeloaderHandle.write("%s %s\t" % (featureName, featureNum))
            typeloaderHandle.write("%s\t%s\n" % (coordinates[0], coordinates[1]))

        typeloaderHandle.write("\nDifferences to Closest Known Allele\n")
        if typeloaderData["isExactMatch"]: typeloaderHandle.write("None\n")
        else:
            diffData = typeloaderData["differences"]

            if len(diffData["mismatchPositions"]):
                mismatches, mismatchPositions = diffData["mismatches"], diffData["mismatchPositions"]
                typeloaderHandle.write("Mismatches : ")
                for mmIndex in range(len(mismatches)):
                    ##print mismatches, mismatchPositions
                    mmLocationData, mmBases = mismatchPositions[mmIndex], mismatches[mmIndex]
                    mmPos, codon = mmLocationData
                    if not codon:
                        typeloaderHandle.write("%s:%s->%s(Intronic)\t" % (mmPos,mmBases[1],mmBases[0]))
                    else:
                        typeloaderHandle.write("%s:%s->%s(Codon %s)\t" % (mmPos,mmBases[1],mmBases[0], codon))
                typeloaderHandle.write("\n")

            if len(diffData["deletionPositions"]):
                typeloaderHandle.write("Deletions : ")
                for deletion in diffData["deletionPositions"]:
                    deletionPos, codon = deletion
                    if not codon:
                        typeloaderHandle.write("%s(Intronic)\t" % deletionPos)
                    else:
                        typeloaderHandle.write("%s(Codon %s)\t" % (deletionPos, codon))
                typeloaderHandle.write("\n")

            if len(diffData["insertionPositions"]):
                typeloaderHandle.write("Insertions : ")
                for insertion in diffData["insertionPositions"]:
                    insertionPos, codon = insertion
                    if not codon:
                        typeloaderHandle.write("%s(Intronic)\t" % insertionPos)
                    else:
                        typeloaderHandle.write("%s(Codon %s)\t" % (insertionPos, codon))
                typeloaderHandle.write("\n")

        typeloaderHandle.write("\n" * 3)
    typeloaderHandle.close()
    """
