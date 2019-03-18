#!/usr/bin/env Python3
# -*- coding: cp1252 -*-
'''
test_typeloader_GUI.py

unit tests for typeloader_GUI

@author: Bianca Schoene
'''

import unittest
from unittest.mock import patch
import os, sys, re, time, platform, datetime, csv
import difflib # compare strings
import shutil
from random import randint
from configparser import ConfigParser

module_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
mypath = os.path.join(module_path)
mypath_inner = os.path.join(mypath, "src")
sys.path.append(mypath)
sys.path.append(mypath_inner)

import general, db_internal, GUI_login
from __init__ import __version__
from xml.etree import ElementTree
from collections import namedtuple

# no .pyw import possibile in linux
# deletion in Test_Clean_Stuff
shutil.copyfile(os.path.join(mypath_inner, "typeloader_GUI.pyw"), os.path.join(mypath_inner, "typeloader_GUI.py"))

import typeloader_GUI
from typeloader_core import errors, EMBLfunctions as EF, make_imgt_files as MIF, backend_make_ena as BME, imgt_text_generator as ITG
import GUI_forms_new_project as PROJECT
import GUI_forms_new_allele as ALLELE
import GUI_forms_submission_ENA as ENA
import GUI_forms_submission_IPD as IPD
import GUI_views_OVprojects, GUI_views_OValleles, GUI_views_project, GUI_views_sample
import typeloader_functions
from GUI_login import base_config_file

from PyQt5.QtWidgets import (QApplication) 
from PyQt5.QtCore import Qt


#===========================================================
# test parameters:

samples_dic =  {# samples to test 
                "sample_1" : { "input" : "1395777_A.fa",
                               "input_dir_origin" : "KIR_3DP1",
                               "local_name" : "DKMS-LSL_ID000001_3DP1_1",
                               "cell_line" : "DKMS-LSL_ID000001",
                               "gene" : "KIR3DP1",
                               "target_allele" : "KIR3DP1*0030102:new",
                               "data_unittest_dir" : "new_allele_fasta",
                               "curr_ena_file" : "DKMS-LSL_ID000001_3DP1_1.ena.txt",
                               "curr_fasta_file" : "DKMS-LSL_ID000001_3DP1_1.fa",
                               "curr_blast_file" : "DKMS-LSL_ID000001_3DP1_1.blast.xml",
                               "curr_ipd_befund_file" : "Befunde_neu_1.csv",
                               "curr_ipd_ena_acc_file" : "ENA_Accession_3DP1",
                               "id_int" : "ID000001",
                               "id_ext" : "DEDKM000001",
                               "submission_id" : "1111"},
                "sample_2" : { "input" : "ID14278154.xml",
                               "input_dir_origin" : "C_MM",
                               "local_name" : "DKMS-LSL_ID14278154_C_1",
                               "cell_line" : "DKMS-LSL_ID14278154",
                               "gene" : "HLA-C",
                               "target_allele" : "HLA-C*16:new",
                               "partner_allele" : 'HLA-C*06:new',
                               "data_unittest_dir" : "new_allele_xml",
                               "curr_ena_file" : "DKMS-LSL_ID14278154_C_1.ena.txt",
                               "curr_fasta_file" : "DKMS-LSL_ID14278154_C_1.fa",
                               "curr_blast_file" : "DKMS-LSL_ID14278154_C_1.blast.xml",
                               "curr_gendx_file" : "DKMS-LSL_ID14278154_C_1.xml",
                               "curr_ipd_befund_file" : "Befunde_C12_1.csv",
                               "curr_ipd_ena_acc_file" : "AccessionNumbers_C12_LT989974-LT990037_1",                               
                               "id_int" : "ID14278154",
                               "id_ext" : "1348480",
                               "submission_id" : "2222"},
                "sample_3" : { "input_dir_origin" : "confirmation_file",
                               "local_name" : "DKMS-LSL_ID15390636_3DP1_1",
                               "cell_line" : "DKMS-LSL_ID15390636",
                               "curr_ipd_befund_file" : "Befunde_3DP1_1.csv",
                               "curr_ipd_ena_acc_file" : "ENA_Accession_3DP1_1",
                               "blast_file_name" : "DKMS-LSL_ID15390636_3DP1_1.blast.xml",
                               "ena_file_name" : "DKMS-LSL_ID15390636_3DP1_1.ena.txt",
                               "id_int" : "ID15390636",
                               "id_ext" : "1370324_A",
                               "submission_id" : "3333",
                               "gene" : "KIR3DP1",
                               "target_allele" : 'KIR3DL3*006',
                               "partner_allele" : 'KIR3DL3*003'}                              
                }

settings_both = {"reference_dir" : "reference_data_unittest",
                "data_unittest" : "data_unittest",
                "ipd_submissions" : "IPD-submissions",
                }

# deletes database entries and project directory
delete_all_stuff_at_the_end = True
skip_other_tests = False # can be set to True to skip all tests except the one currently worked at (out-comment it there in setUpClass)

log = general.start_log(level="DEBUG")

# assemble settings for testing:
cf = ConfigParser()
cf.read(os.path.join(mypath_inner, base_config_file))
curr_settings = GUI_login.get_settings("staging", log, cf)
for key in settings_both:
    curr_settings[key] = settings_both[key]
curr_settings["embl_submission"] = curr_settings["embl_submission_test"]

embl_test_server = curr_settings["embl_submission"]
center_name = curr_settings["xml_center_name"]
mydb = typeloader_GUI.create_connection(log, curr_settings["db_file"])

# ipd_counter_lock_file = os.path.join(curr_settings["root_path"], "_general", "ipd_nr.lock")
# counter_config_file = os.path.join(curr_settings["root_path"], "_general", "counter_config.ini")
# _, (IPD_counter, _) = MIF.get_IPD_counter(counter_config_file, ipd_counter_lock_file, curr_settings, log)

project_gene = "A"
project_pool = str(randint(1,999999))
project_user = "Staging Account"
project_title = "This is an optional title information"
project_desc = "This is an optional description information"
project_name = ""  ## this will be set in create project
project_accession = "" ## this will be set in create project

app = QApplication(sys.argv)

today = datetime.datetime.now()

TargetAllele = namedtuple("TargetAllele", "gene target_allele partner_allele")

#===========================================================
# test cases:

class Test_Clean_Stuff_initial(unittest.TestCase):
    """ 
    Remove all directories and files written by previous unit tests
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping initial cleanup because skip_other_tests is set to True")
        
    @classmethod
    def tearDownClass(self):
        pass
        
    def test_clean_everything(self):
        delete_written_samples(True, "ALLELES", log)
        delete_written_samples(True, "FILES", log)
        delete_written_samples(True, "SAMPLES", log)
        delete_written_samples(True, "PROJECTS", log)
        delete_written_samples(True, "ENA_SUBMISSIONS", log)
        delete_written_samples(True, "IPD_SUBMISSIONS", log)
                
        try:
            shutil.rmtree(os.path.join(curr_settings["projects_dir"], project_name))
        except IOError:
            pass
        

class Test_Create_Project(unittest.TestCase):
    """ create project
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Create Test because skip_other_tests is set to True")
        else:
            self.form = PROJECT.NewProjectForm(log, mydb, curr_settings)
            self.form.project_dir = "" #not initially set in form
             
    @classmethod
    def tearDownClass(self):
        pass

    #@unittest.skip("demonstrating skipping")        
    def test_create_project(self):
        """
        Defines a project name and creates a project on ENA Test Server
        """                
        self.form.gene_entry.setText(project_gene)
        self.form.pool_entry.setText(project_pool)
        self.form.user_entry.setText(project_user)
        self.form.title_entry.setText(project_title)
        self.form.desc_entry.setText(project_desc)
        
        self.form.project_btn.click()
        self.form.submit_btn.click()
        self.form.accession = self.form.acc_entry.text()
        
        ## set global var for further tests
        global project_name
        global project_accession
        project_name = self.form.project_name
        project_accession = self.form.accession
    
    #@unittest.skip("demonstrating skipping")   
    def test_parse_project_name(self):
        """
        parse project name
        """
        date = general.timestamp("%Y%m%d")
        new_project = "_".join([date, "SA", project_gene, project_pool])
        self.assertEqual(self.form.project_name, new_project)
    
    #@unittest.skip("demonstrating skipping")
    def test_dir_exists(self):
        """
        If project dir is created?
        """
        self.assertTrue(os.path.exists(self.form.project_dir))
    
    #@unittest.skip("demonstrating skipping")
    def test_projectfiles_exists(self):
        """
        If project files are created successfull
        """
        self.form.project_file = os.path.join(self.form.project_dir, self.form.project_name + ".xml")
        self.form.submission_file = os.path.join(self.form.project_dir, self.form.project_name + "_sub.xml")
        self.form.output_file = os.path.join(self.form.project_dir, self.form.project_name + "_output.xml")
        self.assertTrue(os.path.exists(self.form.project_file))
        self.assertTrue(os.path.exists(self.form.submission_file))
        self.assertTrue(os.path.exists(self.form.output_file))
    
    #@unittest.skip("demonstrating skipping")
    def test_parse_project_xml(self):
        """
        Parse the written project XML file
        """        
        xml_stuff = ElementTree.parse(self.form.project_file)
        root = xml_stuff.getroot()
        
        self.assertEqual(root[0].attrib["alias"], self.form.project_name)
        self.assertEqual(root[0].attrib["center_name"], center_name)
        self.assertEqual(root[0][0].text, project_title)
        self.assertEqual(root[0][1].text, project_desc)
        
    #@unittest.skip("demonstrating skipping")
    def test_parse_submission_xml(self):
        """
        Parse the written submission XML file
        """        
        xml_stuff = ElementTree.parse(self.form.submission_file)
        root = xml_stuff.getroot()
        
        self.assertEqual(root.attrib["alias"], self.form.project_name + "_sub")
        self.assertEqual(root.attrib["center_name"], center_name)
        self.assertEqual(root[0][0][0].attrib["schema"], "project")
        self.assertEqual(root[0][0][0].attrib["source"], self.form.project_name + ".xml")

    #@unittest.skip("demonstrating skipping")
    def test_parse_output_xml(self):
        """
        Parse the written output XML file
        """        
        xml_stuff = ElementTree.parse(self.form.output_file)
        root = xml_stuff.getroot()
        
        self.assertEqual(root.attrib["submissionFile"], self.form.project_name + "_sub.xml")
        self.assertEqual(root.attrib["success"], "true")
        self.assertEqual(root[0].attrib["accession"], self.form.accession)
        self.assertEqual(root[1].attrib["alias"], self.form.project_name + "_sub")
        self.assertEqual(root[2][0].text, "Submission has been committed.")
        

class Test_Create_New_Allele(unittest.TestCase):
    """ 
    create new allele
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Create New Allele because skip_other_tests is set to True")
        else:
            self.project_name = project_name #"20180710_SA_A_1292" 
            
            self.new_sample_dir_path_1 = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"])
            self.new_ena_file_path_1 = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["curr_ena_file"])
            self.new_fasta_file_path_1 = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["curr_fasta_file"])
            self.new_blast_file_path_1 = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["curr_blast_file"])
            self.reference_file_path_1 = os.path.join(curr_settings["login_dir"], "data_unittest", samples_dic["sample_1"]["curr_ena_file"])        
            
            self.new_sample_dir_path_2 = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"])
            self.new_ena_file_path_2 = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["curr_ena_file"])
            self.new_fasta_file_path_2 = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["curr_fasta_file"])
            self.new_blast_file_path_2 = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["curr_blast_file"])
            self.reference_file_path_2 = os.path.join(curr_settings["login_dir"], "data_unittest", samples_dic["sample_2"]["curr_ena_file"])                
            
    @classmethod
    def tearDownClass(self):
        pass
        
    #@unittest.skip("skipping test_fasta_file")        
    def test_fasta_file(self):
        """
        Ceate ENA flatfile from fasta
        """ 
        self.form = ALLELE.NewAlleleForm(log, mydb, self.project_name, curr_settings, None, samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["id_ext"])
        self.form.file_widget.field.setText(os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], samples_dic["sample_1"]["data_unittest_dir"], samples_dic["sample_1"]["input"]))

        self.form.upload_btn.setEnabled(True)
        self.form.upload_btn.click()

        self.form.save_btn.click()
       
        self.assertEqual(self.form.name_lbl.text(), samples_dic["sample_1"]["target_allele"]) 
        
        new_ena_file_path = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["curr_ena_file"])
        reference_file_path = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], samples_dic["sample_1"]["data_unittest_dir"], samples_dic["sample_1"]["curr_ena_file"])
        diff_ena_files = compare_2_files(new_ena_file_path, reference_file_path)
        
        self.assertEqual(len(diff_ena_files["added_sings"]), 0)
        self.assertEqual(len(diff_ena_files["deleted_sings"]), 0)

        self.assertTrue(os.path.exists(os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"])))
        self.assertTrue(os.path.exists(os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["curr_ena_file"])))
        self.assertTrue(os.path.exists(os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["curr_fasta_file"])))
        self.assertTrue(os.path.exists(os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["curr_blast_file"])))
     
    #@unittest.skip("skipping test_fasta_file")          
    def test_xml_file(self):
        """
        Ceate ENA flatfile from xml
        """ 
        self.form = ALLELE.NewAlleleForm(log, mydb, self.project_name, curr_settings, None, samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["id_ext"])
        self.form.file_widget.field.setText(os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], samples_dic["sample_2"]["data_unittest_dir"], samples_dic["sample_2"]["input"]))
        self.form.upload_btn.setEnabled(True)
        self.form.upload_btn.click()
        
        self.assertEqual(self.form.allele2_sec.gene_field.text(), "HLA-C")
        self.assertEqual(self.form.allele2_sec.GenDX_result, "C*16:02:01-Novel-2")
        self.assertEqual(self.form.allele2_sec.name_field.text(), samples_dic["sample_2"]["target_allele"])
        self.assertEqual(self.form.allele2_sec.product_field.text(), "MHC class I antigen")
        self.assertEqual(int(self.form.allele2_sec.exon1_field.text()), 0)
        self.assertEqual(int(self.form.allele2_sec.exon2_field.text()), 0)
        
        self.form.ok_btn.click()
        self.form.save_btn.click()        
        
        new_ena_file_path = os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["curr_ena_file"])
        reference_file_path = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], samples_dic["sample_2"]["data_unittest_dir"], samples_dic["sample_2"]["curr_ena_file"])
        diff_ena_files = compare_2_files(new_ena_file_path, reference_file_path)
        
        self.assertEqual(len(diff_ena_files["added_sings"]), 0)
        self.assertEqual(len(diff_ena_files["deleted_sings"]), 0)

        self.assertTrue(os.path.exists(os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"])))
        self.assertTrue(os.path.exists(os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["curr_ena_file"])))
        self.assertTrue(os.path.exists(os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["curr_fasta_file"])))
        self.assertTrue(os.path.exists(os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["curr_blast_file"])))        
        self.assertTrue(os.path.exists(os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["curr_gendx_file"])))        

    #@unittest.skip("skipping test_fasta_alleles_entries")        
    def test_fasta_alleles_entries(self):
        """
        Test if ALLELES has the correct entry
        """ 
        query = "PRAGMA table_info (ALLELES)"
        success, data_info = execute_db_query(query, 
                                         6, 
                                         log, 
                                         "Column count at {}", 
                                         "Successful read table_info at {}", 
                                         "Can't get information from {}", 
                                         "ALLELES")
        
        query = "SELECT * from ALLELES"
        success, data_content = execute_db_query(query, 
                                         len(data_info), 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "ALLELES")        
        
        self.assertEqual(len(data_content), 2) # should be 2 - fasta & xml file
        # sample 1: fasta
        self.assertEqual(data_content[0][0], samples_dic["sample_1"]["id_int"]) # sample_id_in
        self.assertEqual(data_content[0][1], 1) # allele_nr
        self.assertEqual(data_content[0][2], self.project_name) # project_name
        self.assertEqual(data_content[0][3], 1) # project_nr
        self.assertEqual(data_content[0][4], "") # old cell_line
        self.assertEqual(data_content[0][5], "{}_{}_{}".format(curr_settings["cell_line_token"],
                                                                  samples_dic["sample_1"]["id_int"], "3DP1_1")) # local_name
        self.assertEqual(data_content[0][6], "KIR3DP1") # gene
        self.assertEqual(data_content[0][7], "novel") # goal
        self.assertEqual(data_content[0][8], "ENA-ready") # allele_status
        self.assertEqual(data_content[0][14], "completed") # lab_status
        self.assertEqual(data_content[0][24], samples_dic["sample_1"]["target_allele"]) # target_allele
        self.assertEqual(data_content[0][31], "IPD-KIR") # reference_database
        # sample 2: xml
        self.assertEqual(data_content[1][0], samples_dic["sample_2"]["id_int"]) # sample_id_in
        self.assertEqual(data_content[1][1], 1) # allele_nr
        self.assertEqual(data_content[1][2], self.project_name) # project_name
        self.assertEqual(data_content[1][3], 2) # project_nr
        self.assertEqual(data_content[1][4], "") # old cell_line
        self.assertEqual(data_content[1][5], "{}_{}_{}".format(curr_settings["cell_line_token"],
                                                                  samples_dic["sample_2"]["id_int"], "C_1")) # local_name
        
        self.assertEqual(data_content[1][6], "HLA-C") # gene
        self.assertEqual(data_content[1][7], "novel") # goal
        self.assertEqual(data_content[1][8], "ENA-ready") # allele_status
        self.assertEqual(data_content[1][14], "completed") # lab_status
        self.assertEqual(data_content[1][20], "yes") # long_read_data
        self.assertEqual(data_content[1][21], "yes") # long_read_phasing
        self.assertEqual(data_content[1][24], samples_dic["sample_2"]["target_allele"]) # target_allele
        self.assertEqual(data_content[1][25], "HLA-C*06:new") # partner_allele
        self.assertEqual(data_content[1][28], "NGSengine") # new_genotyping_software
        self.assertEqual(data_content[1][29], "2.7.0.9307") # new_software_version
        self.assertEqual(data_content[1][30], "2018-03-08") # new_genotyping_date
        self.assertEqual(data_content[1][31], "IPD-IMGT/HLA") # reference_database
        
    #@unittest.skip("skipping test_fasta_files_entries")        
    def test_fasta_files_entries(self):
        """
        Test if FILES has the correct entry
        """ 
        query = "PRAGMA table_info (FILES)"
        success, data_info = execute_db_query(query, 
                                         6, 
                                         log, 
                                         "Column count at {}", 
                                         "Successful read table_info at {}", 
                                         "Can't get information from {}", 
                                         "FILES")
        
        query = "SELECT * from FILES"
        success, data_content = execute_db_query(query, 
                                         len(data_info), 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "FILES")        

        self.assertEqual(len(data_content), 2) # should be 2 - fasta & xml file
        # sample_1: fasta
        self.assertEqual(data_content[0][0], samples_dic["sample_1"]["id_int"]) # sample_id_in
        self.assertEqual(data_content[0][1], 1) # allele_nr
        self.assertEqual(data_content[0][3], self.project_name) # project_name
        self.assertEqual(data_content[0][4], "FASTA") # raw_file_type
        self.assertEqual(data_content[0][5], samples_dic["sample_1"]["curr_fasta_file"]) # raw_file
        self.assertEqual(data_content[0][6], samples_dic["sample_1"]["curr_fasta_file"]) # fasta
        self.assertEqual(data_content[0][7], samples_dic["sample_1"]["curr_blast_file"]) # blast_xml
        self.assertEqual(data_content[0][8], samples_dic["sample_1"]["curr_ena_file"]) # ena_file
        # sample_2: xml
        self.assertEqual(data_content[1][0], samples_dic["sample_2"]["id_int"]) # sample_id_in
        self.assertEqual(data_content[1][1], 1) # allele_nr
        self.assertEqual(data_content[1][3], self.project_name) # project_name
        self.assertEqual(data_content[1][4], "XML") # raw_file_type
        self.assertEqual(data_content[1][5], samples_dic["sample_2"]["curr_gendx_file"]) # raw_file
        self.assertEqual(data_content[1][6], samples_dic["sample_2"]["curr_fasta_file"]) # fasta
        self.assertEqual(data_content[1][7], samples_dic["sample_2"]["curr_blast_file"]) # blast_xml
        self.assertEqual(data_content[1][8], samples_dic["sample_2"]["curr_ena_file"]) # ena_file        
    
    #@unittest.skip("skipping test_fasta_samples_entries")        
    def test_fasta_samples_entries(self):
        """
        Test if SAMPLES has the correct entry
        """ 
        query = "PRAGMA table_info (SAMPLES)"
        success, data_info = execute_db_query(query, 
                                         6, 
                                         log, 
                                         "Column count at {}", 
                                         "Successful read table_info at {}", 
                                         "Can't get information from {}", 
                                         "SAMPLES")
        
        query = "SELECT * from SAMPLES"
        success, data_content = execute_db_query(query, 
                                         len(data_info), 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "SAMPLES")        
        
        self.assertEqual(len(data_content), 2) # should be 2 - fasta & xml file
        # sample_1: fasta
        self.assertEqual(data_content[0][0], samples_dic["sample_1"]["id_int"]) # sample_id_in
        self.assertEqual(data_content[0][1], samples_dic["sample_1"]["id_ext"]) # sample_id_ext
        # sample_2: xml
        self.assertEqual(data_content[1][0], samples_dic["sample_2"]["id_int"]) # sample_id_in
        self.assertEqual(data_content[1][1], samples_dic["sample_2"]["id_ext"]) # sample_id_ext        
        

class Test_Send_To_ENA(unittest.TestCase):
    """ 
    Send both of the fasta and xml smaples to ENA
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Submission to ENA because skip_other_tests is set to True")
        else:
            self.project_name = project_name # "20180710_SA_A_1292" 
            self.form = ENA.ENASubmissionForm(log, mydb, self.project_name, curr_settings, parent = None)
        
    @classmethod
    def tearDownClass(self):
        pass

    #@unittest.skip("demonstrating skipping")        
    def test_submit_to_ENA(self):
        """
        Takes the 2 samples out of Test_Create_New_Allele
        """                
        # ok button clickable, if project was choosen [setUpClass]
        self.form.ok_btn.click()
        
        # they are activated by initialisation, but to be sure...
        self.form.project_files.check_dic[0].setChecked(True)
        self.form.project_files.check_dic[1].setChecked(True)
        
        self.assertEqual(self.form.project_files.item(0,1).text(), "1")
        self.assertEqual(self.form.project_files.item(0,2).text(), samples_dic["sample_1"]["id_int"])
        self.assertEqual(self.form.project_files.item(0,3).text(), samples_dic["sample_1"]["curr_ena_file"])
        self.assertEqual(self.form.project_files.item(0,4).text(), "ENA-ready")
        
        self.assertEqual(self.form.project_files.item(1,1).text(), "2")
        self.assertEqual(self.form.project_files.item(1,2).text(), samples_dic["sample_2"]["id_int"])
        self.assertEqual(self.form.project_files.item(1,3).text(), samples_dic["sample_2"]["curr_ena_file"])
        self.assertEqual(self.form.project_files.item(1,4).text(), "ENA-ready")
        
        # submit = send to ENA        
        self.form.submit_btn.click()
        # do not write in database, if close btn isn't clicked
        self.form.close_btn.click()
        
    #@unittest.skip("demonstrating skipping")
    def test_parse_ena_analysis_xml(self):
        """
        Parse the written ena analysis XML file
        """
        path = os.path.join(curr_settings["projects_dir"], self.project_name)
        # neccessary, because timestep is not known
        analysis_file = list(filter(lambda x: re.search(r'^PRJEB.*analysis.xml', x), os.listdir(path)))[0]
        file_split = analysis_file.split("_")
        analysis_file_path = os.path.join(path, analysis_file)
        curr_alias = "_".join([file_split[0], file_split[1]])
        curr_flatfile = list(filter(lambda x: re.search(r'^PRJEB.*txt.gz', x), os.listdir(path)))[0]        
        
        xml_stuff = ElementTree.parse(analysis_file_path)
        root = xml_stuff.getroot()
        
        self.assertEqual(root[0].attrib["alias"], curr_alias)
        self.assertEqual(root[0].attrib["center_name"], center_name)
        self.assertEqual(root[0][0].text, project_title)
        self.assertEqual(root[0][1].text, project_desc)
        self.assertEqual(root[0][4][0].attrib["checksum"], EF.make_md5(os.path.join(path, curr_flatfile), log))
        self.assertEqual(root[0][4][0].attrib["checksum_method"], "MD5")
        self.assertEqual(root[0][4][0].attrib["filename"], curr_flatfile)
        self.assertEqual(root[0][4][0].attrib["filetype"], "flatfile")
    
    #@unittest.skip("demonstrating skipping")    
    def test_parse_ena_submission_xml(self):
        """
        Parse the written ena submission XML file
        """        
        path = os.path.join(curr_settings["projects_dir"], self.project_name)
        # neccessary, because timestep is not known
        submission_file = list(filter(lambda x: re.search(r'^PRJEB.*submission.xml', x), os.listdir(path)))[0]
        file_split = submission_file.split("_")
        submission_file_path = os.path.join(path, submission_file)
        curr_alias = "_".join([file_split[0], file_split[1], "filesub"])
        curr_source = "_".join([file_split[0], file_split[1], "analysis.xml"])
        
        xml_stuff = ElementTree.parse(submission_file_path)
        root = xml_stuff.getroot()
        
        self.assertEqual(root.attrib["alias"], curr_alias)
        self.assertEqual(root.attrib["center_name"], center_name)
        self.assertEqual(root[0][0][0].attrib["schema"], "analysis")
        self.assertEqual(root[0][0][0].attrib["source"], curr_source)

    #@unittest.skip("demonstrating skipping")
    def test_parse_ena_output_and_db_entry(self):
        """
        Parse the written ena output XML file + the db entry
        """    
        path = os.path.join(curr_settings["projects_dir"], self.project_name)
        # neccessary, because timestep is not known
        output_file = list(filter(lambda x: re.search(r'^PRJEB.*output.xml', x), os.listdir(path)))[0]
        file_split = output_file.split("_")
        output_file_path = os.path.join(path, output_file)  
        curr_sub_id = "_".join([file_split[0], file_split[1]])
        curr_alias = "_".join([file_split[0], file_split[1], "filesub"])
        curr_submissionFile = "_".join([file_split[0], file_split[1], "submission.xml"])  
        
        query = "PRAGMA table_info (ENA_SUBMISSIONS)"
        success, data_info = execute_db_query(query, 
                                         6, 
                                         log, 
                                         "Column count at {}", 
                                         "Successful read table_info at {}", 
                                         "Can't get information from {}", 
                                         "ENA_SUBMISSIONS")
        
        query = "SELECT * from ENA_SUBMISSIONS"
        success, data_content = execute_db_query(query, 
                                         len(data_info), 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "ENA_SUBMISSIONS")        
        
        self.assertEqual(len(data_content), 1)
        self.assertEqual(data_content[0][0], self.project_name) # project_name
        self.assertEqual(data_content[0][1], curr_sub_id ) # submission_id
        self.assertEqual(data_content[0][2], 2) # nr_alleles
        self.assertEqual(data_content[0][7], "yes") # success
        
        acc_analysis = data_content[0][5]
        acc_submission = data_content[0][6]
        
        xml_stuff = ElementTree.parse(output_file_path)
        root = xml_stuff.getroot()
        
        self.assertEqual(root.attrib["submissionFile"], curr_submissionFile)
        self.assertEqual(root.attrib["success"], "true")
        self.assertEqual(root[0].attrib["accession"], acc_analysis) # analysis accession
        self.assertEqual(root[0].attrib["alias"], curr_sub_id) # alias
        self.assertEqual(root[1].attrib["accession"], acc_submission) # submission accession
        self.assertEqual(root[1].attrib["alias"], curr_alias) # alias
        self.assertEqual(root[2][0].text, "Submission has been committed.")    

 
class Test_Send_To_IMGT(unittest.TestCase):
    """ 
    Create IMGT / IPD files
    No transmission
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Submission to IPD because skip_other_tests is set to True")
        else:
            query = "SELECT project_name from projects"
            success, data_content = execute_db_query(query, 
                                         1, 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "PROJECTS")  
            self.assertTrue(success, "Could not retrieve project name from table PROJECTS")
            self.project_name = data_content[0][0]
            self.form = IPD.IPDSubmissionForm(log, mydb, self.project_name, curr_settings, parent = None)
            query = """update ALLELES set IPD_submission_nr = 'DKMS1000{}' 
            where Sample_ID_int = '{}'and allele_nr = 1""".format(samples_dic["sample_1"]["submission_id"],
                                                                  samples_dic["sample_1"]["id_int"])
            execute_db_query(query, 0, log, 
                            "updating {}.IPD_submission_nr", 
                            "Successfully updated {}.IPD_submission_nr", 
                            "Can't update {}.IPD_submission_nr", 
                            "ALLELES")

        
    @classmethod
    def tearDownClass(self):
        pass

    #@unittest.skip("demonstrating skipping")        
    def test_generate_IMGT_files(self):
        """
        Takes the 1st sample out of Test_Create_New_Allele
        """                
        # click to proceed to section 2
        self.form.ok_btn1.click()
        ENA_file = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], 
                                samples_dic["sample_1"]["input_dir_origin"], 
                                samples_dic["sample_1"]["curr_ipd_ena_acc_file"])
        log.debug("ENA-file: {}".format(ENA_file))
        self.form.ENA_file_widget.field.setText(ENA_file)
        befund_file = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], 
                                   samples_dic["sample_1"]["input_dir_origin"], 
                                   samples_dic["sample_1"]["curr_ipd_befund_file"])
        log.debug("befund-file: {}".format(befund_file))
        self.form.befund_widget.field.setText(befund_file)
        
        log.debug("Clicking ok_btn2...")
        self.form.ok_btn2.check_ready()
        self.form.ok_btn2.click()

        # they are activated by initialisation, but to be sure...
        self.form.project_files.check_dic[0].setChecked(True) # if this fails, the alleles in the ENA reply file are probably not recognized correctly
        
        self.assertEqual(self.form.project_files.item(0,1).text(), "1")
        self.assertEqual(self.form.project_files.item(0,2).text(), samples_dic["sample_1"]["id_int"])
        self.assertEqual(self.form.project_files.item(0,3).text(), samples_dic["sample_1"]["local_name"])
        self.assertEqual(self.form.project_files.item(0,4).text(), "ENA submitted")
        query = "SELECT * from ENA_SUBMISSIONS"
        success, data_content = execute_db_query(query, 
                                         2, 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "ENA_SUBMISSIONS")    
        self.assertTrue(success, "Could not retrieve data from ENA_submissions")
        self.assertEqual(self.form.project_files.item(0,5).text(), data_content[0][1]) # submissionID
        self.form.submit_btn.click()
        self.form.ok_btn.click()
    
    def test_updated_alleles_entries(self):
        """
        Test if ALLELES have been updated correctly
        """ 
        query = "SELECT * from ENA_SUBMISSIONS"
        success, data_content_ena = execute_db_query(query, 
                                         2, 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "ENA_SUBMISSIONS")   
        
        query = "SELECT * from IPD_SUBMISSIONS"
        success, data_content_ipd = execute_db_query(query, 
                                         5, 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "IPD_SUBMISSIONS")         
        
        query = "PRAGMA table_info (ALLELES)"
        success, data_info = execute_db_query(query, 
                                         6, 
                                         log, 
                                         "Column count at {}", 
                                         "Successful read table_info at {}", 
                                         "Can't get information from {}", 
                                         "ALLELES")
        
        query = "SELECT * from ALLELES"
        success, data_content = execute_db_query(query, 
                                         len(data_info), 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "ALLELES")        
        
        
        self.assertEqual(len(data_content), 2) # should be 2 - fasta & xml file
        # sample 1: fasta
        self.assertEqual(data_content[0][0], samples_dic["sample_1"]["id_int"]) # sample_id_in
        self.assertEqual(data_content[0][1], 1) # allele_nr
        self.assertEqual(data_content[0][2], self.project_name) # project_name
        self.assertEqual(data_content[0][3], 1) # project_nr
        self.assertEqual(data_content[0][4], "") # old cell_line
        self.assertEqual(data_content[0][5], "{}_{}_{}".format(curr_settings["cell_line_token"],
                                                                  samples_dic["sample_1"]["id_int"], "3DP1_1")) # local_name
        self.assertEqual(data_content[0][6], "KIR3DP1") # gene
        self.assertEqual(data_content[0][7], "novel") # goal
        self.assertEqual(data_content[0][8], "IPD submitted") # allele_status
        self.assertEqual(data_content[0][14], "completed") # lab_status
        self.assertEqual(data_content[0][24], samples_dic["sample_1"]["target_allele"]) # target_allele
        self.assertEqual(data_content[0][31], "IPD-KIR") # reference_database
        self.assertEqual(data_content[0][32], "2.7.1") # database_version
        self.assertEqual(data_content[0][36], data_content_ena[0][1]) # ena_submission_id
        self.assertEqual(data_content[0][38], "LT986596") # ena accession number: LTxxxxxx
        self.assertEqual(data_content[0][39], data_content_ipd[0][0]) # ipd_submission_id
        self.assertEqual(data_content[0][40], "DKMS1000" + samples_dic["sample_1"]["submission_id"]) # ipd_submission_nr
        
        ipd_submission_path = os.path.join(curr_settings["projects_dir"], self.project_name, curr_settings["ipd_submissions"])
        
        # test ipd submission table in db
        self.assertEqual(len(data_content_ipd), 1) # should be 1 - one submission with fasta file
        self.assertEqual(data_content_ipd[0][0], os.listdir(ipd_submission_path)[0]) # submission_id
        self.assertEqual(data_content_ipd[0][1], 1) # number of alleles
        self.assertEqual(data_content_ipd[0][4], "yes") # success
        
        ipd_submission_path = os.path.join(curr_settings["projects_dir"], self.project_name, "IPD-submissions", os.listdir(ipd_submission_path)[0])
        
        # test, if files are written to the IPD-submissions directory
        befund_file = os.path.join(ipd_submission_path, samples_dic["sample_1"]["curr_ipd_befund_file"])
        ipd_output_file = os.path.join(ipd_submission_path, "DKMS1000" + samples_dic["sample_1"]["submission_id"] + ".txt")
        ena_acc_file = os.path.join(ipd_submission_path, samples_dic["sample_1"]["curr_ipd_ena_acc_file"])
        output_zip_file = os.path.join(ipd_submission_path, data_content_ipd[0][0] + ".zip")
        
        self.assertTrue(os.path.exists(befund_file))
        self.assertTrue(os.path.exists(ipd_output_file))
        self.assertTrue(os.path.exists(ena_acc_file))
        self.assertTrue(os.path.exists(output_zip_file))
        
    #@unittest.skip("demonstrating skipping")        
    def test_written_ipd_file(self):
        """
        Test the written IPD file (KIR3DP1 in this case) in the submissions directory
        """         
        
        query = "SELECT * from IPD_SUBMISSIONS"
        success, data_content_ipd = execute_db_query(query, 
                                         1, 
                                         log, 
                                         "Get data from {}", 
                                         "Successful select * from {}", 
                                         "Can't get rows from {}", 
                                         "IPD_SUBMISSIONS")    
        
        self.assertTrue(success, "Could not select data from IPD_SUBMISSIONS")
        
        new_ipd_file_path = os.path.join(curr_settings["projects_dir"], self.project_name, curr_settings["ipd_submissions"], data_content_ipd[0][0], "DKMS1000" + samples_dic["sample_1"]["submission_id"] + ".txt")
        reference_file_path = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], samples_dic["sample_1"]["input_dir_origin"], "DKMS1000" + samples_dic["sample_1"]["submission_id"] + ".txt")
        log.debug("Reference IPD file: {}".format(reference_file_path))
        diff_ipd_files = compare_2_files(new_ipd_file_path, reference_file_path, "IPD")        
        
        self.assertEqual(len(diff_ipd_files["added_sings"]), 0)
        self.assertEqual(len(diff_ipd_files["deleted_sings"]), 0)
        

class Test_Views(unittest.TestCase):
    """ 
    tests whether all data is displayed correctly by all views
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_Views because skip_other_tests is set to True")
        else:
            self.date = general.timestamp("%Y%m%d")
            self.proj_name = "_".join([self.date, "SA", project_gene, project_pool])
            self.settings = curr_settings
            self.views = {"OVproj" : GUI_views_OVprojects.ProjectsOverview(log, mydb, parent = self),
                          "OValleles" : GUI_views_OValleles.AllelesOverview(log, mydb),
                          "projectView" : GUI_views_project.ProjectView(log, mydb, self.proj_name, parent = self),
                          "sampleView" : GUI_views_sample.SampleView(log, mydb, samples_dic["sample_1"]["id_int"], self.proj_name, parent = self)}
        
    @classmethod
    def tearDownClass(self):
        pass

    def test_OV_projects(self):
        """tests whether content of projects overview is correct
        """
        view = self.views["OVproj"]
        model = view.proxy
        num_rows = model.rowCount()
        self.assertEqual(num_rows, 1, "ProjectsOverview contains more than 1 row!")
        
        self.assertEqual(model.headerData(0, Qt.Horizontal, Qt.DisplayRole), "Project Name")
        self.assertEqual(model.data(model.index(0, 0)), self.proj_name, "Project name in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(1, Qt.Horizontal, Qt.DisplayRole), "Project Status")
        self.assertEqual(model.data(model.index(0, 1)), "Open", "Project status in ProjectsOverview unexpected (should be 'Open')")
        self.assertEqual(model.headerData(2, Qt.Horizontal, Qt.DisplayRole), "Creation Date")
        self.assertEqual(model.data(model.index(0, 2)), general.timestamp("%d.%m.%Y"), "Creation Date in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(3, Qt.Horizontal, Qt.DisplayRole), "User Name")
        self.assertEqual(model.data(model.index(0, 3)), project_user, "User name in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(4, Qt.Horizontal, Qt.DisplayRole), "Gene")
        self.assertEqual(model.data(model.index(0, 4)), project_gene, "Gene in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(5, Qt.Horizontal, Qt.DisplayRole), "Pool")
        self.assertEqual(model.data(model.index(0, 5)), project_pool, "Pool in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(6, Qt.Horizontal, Qt.DisplayRole), "Title")
        self.assertEqual(model.data(model.index(0, 6)), project_title, "Project title in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(7, Qt.Horizontal, Qt.DisplayRole), "Description")
        self.assertEqual(model.data(model.index(0, 7)), project_desc, "Project description in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(8, Qt.Horizontal, Qt.DisplayRole), "Number of Alleles")
        self.assertEqual(model.data(model.index(0, 8)), 2, "Number of alleles in ProjectsOverview unexpected")
    
    def test_OV_alleles(self):
        """tests whether content of alleles overview is correct
        """
        view = self.views["OValleles"]
        view.add_headers()
        model = view.proxy
        num_rows = model.rowCount()
        self.assertEqual(num_rows, 2, "AllelesOverviews does not contain 2 rows!")

        self.assertEqual(model.headerData(0, Qt.Horizontal, Qt.DisplayRole), "Internal Sample ID")
        self.assertEqual(model.data(model.index(0, 0)), samples_dic["sample_1"]["id_int"])
        self.assertEqual(model.data(model.index(1, 0)), samples_dic["sample_2"]["id_int"])
        self.assertEqual(model.headerData(1, Qt.Horizontal, Qt.DisplayRole), "Allele Nr. in Sample")
        self.assertEqual(model.data(model.index(0, 1)), 1)
        self.assertEqual(model.data(model.index(1, 1)), 1)
        self.assertEqual(model.headerData(2, Qt.Horizontal, Qt.DisplayRole), "Project Name")
        self.assertEqual(model.data(model.index(0, 2)), self.proj_name)
        self.assertEqual(model.data(model.index(1, 2)), self.proj_name)
        self.assertEqual(model.headerData(3, Qt.Horizontal, Qt.DisplayRole), "Nr. in Project")
        self.assertEqual(model.data(model.index(0, 3)), 1)
        self.assertEqual(model.data(model.index(1, 3)), 2)
        self.assertEqual(model.headerData(4, Qt.Horizontal, Qt.DisplayRole), "Cell Line (Old)")
        self.assertEqual(model.data(model.index(0, 4)), "")
        self.assertEqual(model.data(model.index(1, 4)), "")
        self.assertEqual(model.headerData(5, Qt.Horizontal, Qt.DisplayRole), "Allele Name")
        self.assertEqual(model.data(model.index(0, 5)), samples_dic["sample_1"]["local_name"])
        self.assertEqual(model.data(model.index(1, 5)), samples_dic["sample_2"]["local_name"])
        self.assertEqual(model.headerData(6, Qt.Horizontal, Qt.DisplayRole), "Gene")
        self.assertEqual(model.data(model.index(0, 6)), samples_dic["sample_1"]["gene"])
        self.assertEqual(model.data(model.index(1, 6)), samples_dic["sample_2"]["gene"])
        self.assertEqual(model.headerData(7, Qt.Horizontal, Qt.DisplayRole), "Goal")
        self.assertEqual(model.data(model.index(0, 7)), "novel")
        self.assertEqual(model.data(model.index(1, 7)), "novel")
        self.assertEqual(model.headerData(8, Qt.Horizontal, Qt.DisplayRole), "Allele Status")
        self.assertEqual(model.data(model.index(0, 8)), "IPD submitted")
        self.assertEqual(model.data(model.index(1, 8)), "ENA submitted")
        
        self.assertEqual(model.headerData(9, Qt.Horizontal, Qt.DisplayRole), "Original Allele #1")
        self.assertEqual(model.headerData(10, Qt.Horizontal, Qt.DisplayRole), "Original Allele #2")
        self.assertEqual(model.headerData(11, Qt.Horizontal, Qt.DisplayRole), "Software")
        self.assertEqual(model.headerData(12, Qt.Horizontal, Qt.DisplayRole), "Software Version")
        self.assertEqual(model.headerData(13, Qt.Horizontal, Qt.DisplayRole), "Genotyping Date")
        
        self.assertEqual(model.headerData(14, Qt.Horizontal, Qt.DisplayRole), "Lab Status")
        self.assertEqual(model.data(model.index(0, 14)), 'completed')
        self.assertEqual(model.data(model.index(1, 14)), 'completed')
        self.assertEqual(model.headerData(15, Qt.Horizontal, Qt.DisplayRole), "Panel")
        self.assertEqual(model.headerData(16, Qt.Horizontal, Qt.DisplayRole), "Position")
        self.assertEqual(model.headerData(17, Qt.Horizontal, Qt.DisplayRole), "Short Read Data?")
        self.assertEqual(model.headerData(18, Qt.Horizontal, Qt.DisplayRole), "SR Phased?")
        self.assertEqual(model.headerData(19, Qt.Horizontal, Qt.DisplayRole), "SR Technology")
        self.assertEqual(model.headerData(20, Qt.Horizontal, Qt.DisplayRole), "Long Read Data?")
        self.assertEqual(model.data(model.index(0, 20)), '')
        self.assertEqual(model.data(model.index(1, 20)), 'yes')
        self.assertEqual(model.headerData(21, Qt.Horizontal, Qt.DisplayRole), "LR Phased?")
        self.assertEqual(model.data(model.index(0, 21)), '')
        self.assertEqual(model.data(model.index(1, 21)), 'yes')
        self.assertEqual(model.headerData(22, Qt.Horizontal, Qt.DisplayRole), "LR Technology")
        self.assertEqual(model.headerData(23, Qt.Horizontal, Qt.DisplayRole), "Comment")
        
        self.assertEqual(model.headerData(24, Qt.Horizontal, Qt.DisplayRole), "Target Allele")
        self.assertEqual(model.data(model.index(0, 24)), samples_dic["sample_1"]["target_allele"])
        self.assertEqual(model.data(model.index(1, 24)), samples_dic["sample_2"]["target_allele"])
        self.assertEqual(model.headerData(25, Qt.Horizontal, Qt.DisplayRole), "Partner Allele")
        self.assertEqual(model.data(model.index(0, 25)), "")
        self.assertEqual(model.data(model.index(1, 25)), samples_dic["sample_2"]["partner_allele"])
        self.assertEqual(model.headerData(26, Qt.Horizontal, Qt.DisplayRole), "Mismatch Position")
        self.assertEqual(model.headerData(27, Qt.Horizontal, Qt.DisplayRole), "Null Allele?")
        self.assertEqual(model.data(model.index(0, 27)), "no")
        self.assertEqual(model.headerData(28, Qt.Horizontal, Qt.DisplayRole), "Software (new)")
        self.assertEqual(model.data(model.index(0, 28)), "")
        self.assertEqual(model.data(model.index(1, 28)), "NGSengine")
        self.assertEqual(model.headerData(29, Qt.Horizontal, Qt.DisplayRole), "Software Version")
        self.assertEqual(model.data(model.index(0, 29)), "")
        self.assertEqual(model.data(model.index(1, 29)), "2.7.0.9307")
        self.assertEqual(model.headerData(30, Qt.Horizontal, Qt.DisplayRole), "Genotyping Date")
        self.assertEqual(model.data(model.index(0, 30)), "")
        self.assertEqual(model.data(model.index(1, 30)), "2018-03-08")
        self.assertEqual(model.headerData(31, Qt.Horizontal, Qt.DisplayRole), "Reference Database")
        self.assertEqual(model.data(model.index(0, 31)), "IPD-KIR")
        self.assertEqual(model.data(model.index(1, 31)), "IPD-IMGT/HLA")
        self.assertEqual(model.headerData(32, Qt.Horizontal, Qt.DisplayRole), "Database Version") # will continually change, therefore not testing content
        self.assertEqual(model.headerData(33, Qt.Horizontal, Qt.DisplayRole), "Internal Allele Name")
        self.assertEqual(model.headerData(34, Qt.Horizontal, Qt.DisplayRole), "Official Allele Name")
        self.assertEqual(model.headerData(35, Qt.Horizontal, Qt.DisplayRole), "New or confirmed?")
        
        self.assertEqual(model.headerData(36, Qt.Horizontal, Qt.DisplayRole), "ENA Submission ID") # will continually change, therefore not testing content
        self.assertEqual(model.headerData(37, Qt.Horizontal, Qt.DisplayRole), "ENA Acception Date")
        self.assertEqual(model.data(model.index(0, 37)), "2018-07-10")
        self.assertEqual(model.data(model.index(1, 37)), "")
        self.assertEqual(model.headerData(38, Qt.Horizontal, Qt.DisplayRole), "ENA Accession Nr")
        self.assertEqual(model.data(model.index(0, 38)), "LT986596")
        self.assertEqual(model.data(model.index(1, 38)), "")
        
        self.assertEqual(model.headerData(39, Qt.Horizontal, Qt.DisplayRole), "IPD Submission ID") # will continually change, therefore not testing content
        self.assertEqual(model.headerData(40, Qt.Horizontal, Qt.DisplayRole), "IPD Submission Nr")
        self.assertEqual(model.data(model.index(0, 40)), "DKMS1000{}".format(samples_dic["sample_1"]["submission_id"]))
        self.assertEqual(model.data(model.index(1, 40)), "")
        self.assertEqual(model.headerData(41, Qt.Horizontal, Qt.DisplayRole), "HWS Submission Nr")
        self.assertEqual(model.headerData(42, Qt.Horizontal, Qt.DisplayRole), "IPD Acception Date")
        self.assertEqual(model.headerData(43, Qt.Horizontal, Qt.DisplayRole), "IPD Release")
        self.assertEqual(model.headerData(44, Qt.Horizontal, Qt.DisplayRole), "Upload Date")
        self.assertEqual(model.data(model.index(0,44)), general.timestamp("%Y-%m-%d"))
        self.assertEqual(model.headerData(45, Qt.Horizontal, Qt.DisplayRole), "Detection Date")
        
        self.assertEqual(model.headerData(46, Qt.Horizontal, Qt.DisplayRole), "SAMPLE_ID_INT")
        self.assertEqual(model.headerData(47, Qt.Horizontal, Qt.DisplayRole), "External Sample ID")
        self.assertEqual(model.headerData(48, Qt.Horizontal, Qt.DisplayRole), "Cell Line")
        self.assertEqual(model.headerData(49, Qt.Horizontal, Qt.DisplayRole), "Customer")
        self.assertEqual(model.headerData(50, Qt.Horizontal, Qt.DisplayRole), "Project Name")
        self.assertEqual(model.headerData(51, Qt.Horizontal, Qt.DisplayRole), "ENA Submission ID")
        self.assertEqual(model.headerData(52, Qt.Horizontal, Qt.DisplayRole), "Alleles in ENA Submission")
        self.assertEqual(model.headerData(53, Qt.Horizontal, Qt.DisplayRole), "Timestamp Sent (ENA Submission)")
        self.assertEqual(model.headerData(54, Qt.Horizontal, Qt.DisplayRole), "Timestamp Confirmed (ENA Submission)")
        self.assertEqual(model.headerData(55, Qt.Horizontal, Qt.DisplayRole), "Analysis Accession Nr")
        self.assertEqual(model.headerData(56, Qt.Horizontal, Qt.DisplayRole), "Submission Accession Nr")
        self.assertEqual(model.headerData(57, Qt.Horizontal, Qt.DisplayRole), "ENA Submission successful?")
        self.assertEqual(model.headerData(58, Qt.Horizontal, Qt.DisplayRole), "IPD Submission ID")
        self.assertEqual(model.headerData(59, Qt.Horizontal, Qt.DisplayRole), "Alleles in IPD Submission")
        self.assertEqual(model.headerData(60, Qt.Horizontal, Qt.DisplayRole), "Timestamp Ready (IPD Submission)")
        self.assertEqual(model.headerData(61, Qt.Horizontal, Qt.DisplayRole), "Timestamp Confirmed (IPD Submission)")
        self.assertEqual(model.headerData(62, Qt.Horizontal, Qt.DisplayRole), "IPD Submission successful?")

        # check expected empty columns:
        empty_columns = [9, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 26, 
                         33, 34, 35, 41, 42, 43, 61]
        for col in empty_columns:
            for row in [0, 1]:
                self.assertEqual(model.data(model.index(row, col)), "", 
                                 "Unexpected value in column {} ({}) row {}: expected empty cell, found '{}'".format(
                                col, model.headerData(col, Qt.Horizontal, Qt.DisplayRole), row, model.data(model.index(0, col))))

    def test_view_project1_statistics(self):
        """tests whether content of ProjectView table 'Statistics' is correct
        """
        view = self.views["projectView"].project_stats
        model = view.flipped_model
        
        self.assertEqual(model.headerData(1, Qt.Vertical, Qt.DisplayRole), "Number of alleles")
        self.assertEqual(model.data(model.index(1, 0), Qt.DisplayRole), 2)
        self.assertEqual(model.headerData(2, Qt.Vertical, Qt.DisplayRole), "Closed alleles")
        self.assertEqual(model.data(model.index(2, 0), Qt.DisplayRole), 0)
        self.assertEqual(model.headerData(3, Qt.Vertical, Qt.DisplayRole), "Submitted to ENA")
        self.assertEqual(model.data(model.index(3, 0), Qt.DisplayRole), 2)
        self.assertEqual(model.headerData(4, Qt.Vertical, Qt.DisplayRole), "Submitted to IPD")
        self.assertEqual(model.data(model.index(4, 0), Qt.DisplayRole), 1)
        self.assertEqual(model.headerData(5, Qt.Vertical, Qt.DisplayRole), "Accepted by IPD")
        self.assertEqual(model.data(model.index(5, 0), Qt.DisplayRole), 0)
        self.assertEqual(model.headerData(6, Qt.Vertical, Qt.DisplayRole), "Abandoned")
        self.assertEqual(model.data(model.index(5, 0), Qt.DisplayRole), 0)
        
    def test_view_project2_general(self):
        """tests whether content of ProjectView table 'General Information' is correct
        """
        view = self.views["projectView"].project_info
        model = view.flipped_model
#         for row in range(model.rowCount()):
#             foo = [row, model.headerData(row, Qt.Vertical, Qt.DisplayRole)]
#             for col in range(model.columnCount()):
#                 foo.append((col, model.data(model.index(row, col), Qt.DisplayRole)))
#             print (foo)
#             print('self.assertEqual(model.headerData({}, Qt.Vertical, Qt.DisplayRole), "{}")'.format(row, model.headerData(row, Qt.Vertical, Qt.DisplayRole)))
        
        self.assertEqual(model.headerData(1, Qt.Vertical, Qt.DisplayRole), "Project Status")
        self.assertEqual(model.data(model.index(1, 0), Qt.DisplayRole), "Open")
        self.assertEqual(model.headerData(2, Qt.Vertical, Qt.DisplayRole), "Created on")
        self.assertEqual(model.data(model.index(2, 0), Qt.DisplayRole), general.timestamp("%d.%m.%Y"))
        self.assertEqual(model.headerData(3, Qt.Vertical, Qt.DisplayRole), "Created by")
        self.assertEqual(model.data(model.index(3, 0), Qt.DisplayRole), project_user)
        self.assertEqual(model.headerData(4, Qt.Vertical, Qt.DisplayRole), "Gene")
        self.assertEqual(model.data(model.index(4, 0), Qt.DisplayRole), project_gene)
        self.assertEqual(model.headerData(5, Qt.Vertical, Qt.DisplayRole), "Pool")
        self.assertEqual(model.data(model.index(5, 0), Qt.DisplayRole), project_pool)
        self.assertEqual(model.headerData(6, Qt.Vertical, Qt.DisplayRole), "Title")
        self.assertEqual(model.data(model.index(6, 0), Qt.DisplayRole), project_title)
        self.assertEqual(model.headerData(7, Qt.Vertical, Qt.DisplayRole), "Description")
        self.assertEqual(model.data(model.index(7, 0), Qt.DisplayRole), project_desc)
        self.assertEqual(model.headerData(8, Qt.Vertical, Qt.DisplayRole), "ENA Project ID")
        if not skip_other_tests:
            self.assertEqual(model.data(model.index(8, 0), Qt.DisplayRole), project_accession)
        self.assertEqual(model.headerData(9, Qt.Vertical, Qt.DisplayRole), "ENA Project Submission ID") # not testing content, but if project_accession is ok, this will be ok, too

    def test_view_project3_alleles(self):
        """tests whether content of ProjectView table 'Alleles' is correct
        """
        view = self.views["projectView"].alleles
        model = view.proxy
#         for col in range(model.columnCount()):
#             foo = [col, model.headerData(col, Qt.Horizontal, Qt.DisplayRole)]
#             for row in range(model.rowCount()):
#                 foo.append((row, col, model.data(model.index(row, col), Qt.DisplayRole)))
#             print (foo)
#             print('self.assertEqual(model.headerData({}, Qt.Horizontal, Qt.DisplayRole), "{}")'.format(col, model.headerData(col, Qt.Horizontal, Qt.DisplayRole)))
        
        self.assertEqual(model.headerData(1, Qt.Horizontal, Qt.DisplayRole), "Nr")
        self.assertEqual(model.data(model.index(0, 1), Qt.DisplayRole), 1)
        self.assertEqual(model.data(model.index(1, 1), Qt.DisplayRole), 2)
        self.assertEqual(model.headerData(2, Qt.Horizontal, Qt.DisplayRole), "Target Allele")
        self.assertEqual(model.data(model.index(0, 2), Qt.DisplayRole), "{} #{} ({})".format(samples_dic["sample_1"]["id_int"], 1,
                                                                                             samples_dic["sample_1"]["gene"]))
        self.assertEqual(model.data(model.index(1, 2), Qt.DisplayRole), "{} #{} ({})".format(samples_dic["sample_2"]["id_int"], 1,
                                                                                             samples_dic["sample_2"]["gene"]))
        self.assertEqual(model.headerData(3, Qt.Horizontal, Qt.DisplayRole), "Allele Name")
        self.assertEqual(model.data(model.index(0, 3), Qt.DisplayRole), samples_dic["sample_1"]["local_name"])
        self.assertEqual(model.data(model.index(1, 3), Qt.DisplayRole), samples_dic["sample_2"]["local_name"])
        self.assertEqual(model.headerData(4, Qt.Horizontal, Qt.DisplayRole), "Allele Status")
        self.assertEqual(model.data(model.index(0, 4), Qt.DisplayRole), "IPD submitted")
        self.assertEqual(model.data(model.index(1, 4), Qt.DisplayRole), "ENA submitted")
        self.assertEqual(model.headerData(5, Qt.Horizontal, Qt.DisplayRole), "Lab Status")
        self.assertEqual(model.data(model.index(0, 5), Qt.DisplayRole), "completed")
        self.assertEqual(model.data(model.index(1, 5), Qt.DisplayRole), "completed")
        
    def test_view_sample1_general(self):
        """tests whether content of SampleView table 'General Information' is correct
        """
        view = self.views["sampleView"].sample_table
        model = view.flipped_model
        
        self.assertEqual(model.headerData(0, Qt.Vertical, Qt.DisplayRole), "Internal Donor-ID")
        self.assertEqual(model.data(model.index(0, 0), Qt.DisplayRole), samples_dic["sample_1"]["id_int"])
        self.assertEqual(model.headerData(1, Qt.Vertical, Qt.DisplayRole), "External Donor-ID")
        self.assertEqual(model.data(model.index(1, 0), Qt.DisplayRole), samples_dic["sample_1"]["id_ext"])
        self.assertEqual(model.headerData(2, Qt.Vertical, Qt.DisplayRole), "Cell Line")
        self.assertEqual(model.data(model.index(2, 0), Qt.DisplayRole), samples_dic["sample_1"]["cell_line"])
        self.assertEqual(model.headerData(3, Qt.Vertical, Qt.DisplayRole), "Customer")
        self.assertEqual(model.data(model.index(3, 0), Qt.DisplayRole), "DKMSUS")

    def test_view_sample2_alleles(self):
        """tests whether content of SampleView table 'Alleles' is correct
        """
        view = self.views["sampleView"].sample_alleles
        model = view.proxy
        
        self.assertEqual(model.headerData(2, Qt.Horizontal, Qt.DisplayRole), "Target Allele")
        self.assertEqual(model.data(model.index(0, 2), Qt.DisplayRole), "#{} ({})".format(1, samples_dic["sample_1"]["gene"]))
        self.assertEqual(model.headerData(3, Qt.Horizontal, Qt.DisplayRole), "Allele Name")
        self.assertEqual(model.data(model.index(0, 3), Qt.DisplayRole), samples_dic["sample_1"]["local_name"])
        self.assertEqual(model.headerData(4, Qt.Horizontal, Qt.DisplayRole), "Allele Status")
        self.assertEqual(model.data(model.index(0, 4), Qt.DisplayRole), "IPD submitted")
        self.assertEqual(model.headerData(5, Qt.Horizontal, Qt.DisplayRole), "Lab Status")
        self.assertEqual(model.data(model.index(0, 5), Qt.DisplayRole), "completed")
        self.assertEqual(model.headerData(6, Qt.Horizontal, Qt.DisplayRole), "Project")
        self.assertEqual(model.data(model.index(0, 6), Qt.DisplayRole), self.proj_name)
    
    def test_view_sample3_alleles(self):
        """tests whether content of SampleView widget 'Details about Allele' is correct
        """
        widget = self.views["sampleView"].allele_view
        
        def test_tab1_general(self):
            """tests whether content of SampleView widget 'Details about Allele' tab 'General' is correct
            """
            view = widget.tabs[0]
            model = view.flipped_model
            
            # test whether only right columns are shown:
            shown_columns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 14, 33, 34]
            for col in shown_columns:
                self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), False)
            for col in range(44):
                if col not in shown_columns:
                    self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), True)

            # test content and headers of columns:
            self.assertEqual(model.headerData(0, Qt.Vertical, Qt.DisplayRole), "Internal Sample ID")
            self.assertEqual(model.data(model.index(0, 0), Qt.DisplayRole), samples_dic["sample_1"]["id_int"])
            self.assertEqual(model.headerData(1, Qt.Vertical, Qt.DisplayRole), "Allele Nr. in Sample")
            self.assertEqual(model.data(model.index(1, 0), Qt.DisplayRole), 1)
            self.assertEqual(model.headerData(2, Qt.Vertical, Qt.DisplayRole), "Project Name")
            self.assertEqual(model.data(model.index(2, 0), Qt.DisplayRole), self.proj_name)
            self.assertEqual(model.headerData(3, Qt.Vertical, Qt.DisplayRole), "Nr. in Project")
            self.assertEqual(model.data(model.index(3, 0), Qt.DisplayRole), 1)
            self.assertEqual(model.headerData(4, Qt.Vertical, Qt.DisplayRole), "Cell Line (Old)")
            self.assertEqual(model.data(model.index(4, 0), Qt.DisplayRole), "")
            self.assertEqual(model.headerData(5, Qt.Vertical, Qt.DisplayRole), "Allele Name")
            self.assertEqual(model.data(model.index(5, 0), Qt.DisplayRole), samples_dic["sample_1"]["local_name"])
            self.assertEqual(model.headerData(6, Qt.Vertical, Qt.DisplayRole), "Gene")
            self.assertEqual(model.data(model.index(6, 0), Qt.DisplayRole), samples_dic["sample_1"]["gene"])
            self.assertEqual(model.headerData(7, Qt.Vertical, Qt.DisplayRole), "Goal")
            self.assertEqual(model.data(model.index(7, 0), Qt.DisplayRole), "novel")
            self.assertEqual(model.headerData(8, Qt.Vertical, Qt.DisplayRole), "Allele Status")
            self.assertEqual(model.data(model.index(8, 0), Qt.DisplayRole), "IPD submitted")
            self.assertEqual(model.headerData(14, Qt.Vertical, Qt.DisplayRole), "Lab Status")
            self.assertEqual(model.data(model.index(14, 0), Qt.DisplayRole), "completed")
            self.assertEqual(model.headerData(33, Qt.Vertical, Qt.DisplayRole), "Internal Allele Name")
            self.assertEqual(model.data(model.index(33, 0), Qt.DisplayRole), "")
            self.assertEqual(model.headerData(34, Qt.Vertical, Qt.DisplayRole), "Official Allele Name")
            self.assertEqual(model.data(model.index(34, 0), Qt.DisplayRole), "")
            
        def test_tab2_typing_old(self):
            """tests whether content of SampleView widget 'Details about Allele' tab 'Original Genotyping' is correct
            """
            view = widget.tabs[1]
            model = view.flipped_model
            
            # test whether only right columns are shown:
            shown_columns = [9, 10, 11, 12, 13]
            for col in shown_columns:
                self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), False)
            for col in range(44):
                if col not in shown_columns:
                    self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), True)
 
            # test content and headers of columns:
            self.assertEqual(model.headerData(9, Qt.Vertical, Qt.DisplayRole), "Original Allele #1")
            self.assertEqual(model.headerData(10, Qt.Vertical, Qt.DisplayRole), "Original Allele #2")
            self.assertEqual(model.headerData(11, Qt.Vertical, Qt.DisplayRole), "Software")
            self.assertEqual(model.headerData(12, Qt.Vertical, Qt.DisplayRole), "Software Version")
            self.assertEqual(model.headerData(13, Qt.Vertical, Qt.DisplayRole), "Genotyping Date")
            
            for col in shown_columns:
                value = model.data(model.index(col, 0), Qt.DisplayRole)
                msg = "Col {} ({}) should be empty, but contains {}!".format(col, 
                                                                             model.headerData(col, Qt.Vertical, Qt.DisplayRole),
                                                                             value)
                self.assertEqual(value, "", msg)
            
        def test_tab3_lab(self):
            """tests whether content of SampleView widget 'Details about Allele' tab 'Lab Processing' is correct
            """
            view = widget.tabs[2]
            model = view.flipped_model

            # test whether only right columns are shown:
            shown_columns = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
            for col in shown_columns:
                self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), False)
            for col in range(44):
                if col not in shown_columns:
                    self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), True)
            
            # test content and headers of columns:
            self.assertEqual(model.headerData(14, Qt.Vertical, Qt.DisplayRole), "Lab Status")
            self.assertEqual(model.data(model.index(14, 0), Qt.DisplayRole), "completed")
            self.assertEqual(model.headerData(15, Qt.Vertical, Qt.DisplayRole), "Panel")
            self.assertEqual(model.headerData(16, Qt.Vertical, Qt.DisplayRole), "Position")
            self.assertEqual(model.headerData(17, Qt.Vertical, Qt.DisplayRole), "Short Read Data?")
            self.assertEqual(model.headerData(18, Qt.Vertical, Qt.DisplayRole), "SR Phased?")
            self.assertEqual(model.headerData(19, Qt.Vertical, Qt.DisplayRole), "SR Technology")
            self.assertEqual(model.headerData(20, Qt.Vertical, Qt.DisplayRole), "Long Read Data?")
            self.assertEqual(model.headerData(21, Qt.Vertical, Qt.DisplayRole), "LR Phased?")
            self.assertEqual(model.headerData(22, Qt.Vertical, Qt.DisplayRole), "LR Technology")
            self.assertEqual(model.headerData(23, Qt.Vertical, Qt.DisplayRole), "Comment")
            
            for col in [15, 16, 17, 18, 19, 20, 21, 22, 23]:
                value = model.data(model.index(col, 0), Qt.DisplayRole)
                msg = "Col {} ({}) should be empty, but contains {}!".format(col, 
                                                                             model.headerData(col, Qt.Vertical, Qt.DisplayRole),
                                                                             value)
                self.assertEqual(value, "", msg)
        
        def test_tab4_typing_new(self):
            """tests whether content of SampleView widget 'Details about Allele' tab 'New Genotyping' is correct
            """
            view = widget.tabs[3]
            model = view.flipped_model
            
            # test content and headers of columns:
            self.assertEqual(model.headerData(7, Qt.Vertical, Qt.DisplayRole), "Goal")
            self.assertEqual(model.data(model.index(7, 0), Qt.DisplayRole), "novel")
            self.assertEqual(model.headerData(24, Qt.Vertical, Qt.DisplayRole), "Target Allele")
            self.assertEqual(model.data(model.index(24, 0), Qt.DisplayRole), "KIR3DP1*0030102:new")
            self.assertEqual(model.headerData(25, Qt.Vertical, Qt.DisplayRole), "Partner Allele")
            self.assertEqual(model.headerData(26, Qt.Vertical, Qt.DisplayRole), "Mismatch Position")
            self.assertEqual(model.headerData(27, Qt.Vertical, Qt.DisplayRole), "Null Allele?")
            self.assertEqual(model.data(model.index(27, 0), Qt.DisplayRole), "no")
            self.assertEqual(model.headerData(28, Qt.Vertical, Qt.DisplayRole), "Software (new)")
            self.assertEqual(model.headerData(29, Qt.Vertical, Qt.DisplayRole), "Software Version")
            self.assertEqual(model.headerData(30, Qt.Vertical, Qt.DisplayRole), "Genotyping Date")
            self.assertEqual(model.headerData(31, Qt.Vertical, Qt.DisplayRole), "Reference Database")
            self.assertEqual(model.data(model.index(31, 0), Qt.DisplayRole), "IPD-KIR")
            self.assertEqual(model.headerData(32, Qt.Vertical, Qt.DisplayRole), "Database Version") # will change, therefore not testing content
            self.assertEqual(model.headerData(33, Qt.Vertical, Qt.DisplayRole), "Internal Allele Name")
            self.assertEqual(model.headerData(34, Qt.Vertical, Qt.DisplayRole), "Official Allele Name")
            self.assertEqual(model.headerData(35, Qt.Vertical, Qt.DisplayRole), "New or confirmed?")
            
            for col in [25, 26, 28, 29, 30, 33, 34, 35]:
                value = model.data(model.index(col, 0), Qt.DisplayRole)
                msg = "Col {} ({}) should be empty, but contains {}!".format(col, 
                                                                             model.headerData(col, Qt.Vertical, Qt.DisplayRole),
                                                                             value)
                self.assertEqual(value, "", msg)
   
        def test_tab5_ena(self):
            """tests whether content of SampleView widget 'Details about Allele' tab 'ENA submission' is correct
            """
            view = widget.tabs[4]
            model = view.flipped_model
#             # get reference data:
#             shown = []
#             empty = []
#             
#             assert_text = ""
#             for row in range(model.rowCount()):
#                 if not view.table.isIndexHidden(model.index(row, 1)):
#                     assert_text += 'self.assertEqual(model.headerData({}, Qt.Vertical, Qt.DisplayRole), "{}")\n'.format(row, model.headerData(row, Qt.Vertical, Qt.DisplayRole))
#                     data = model.data(model.index(row, 0), Qt.DisplayRole)
#                     if data:
#                         assert_text += 'self.assertEqual(model.data(model.index({}, 0), Qt.DisplayRole), "{}")\n'.format(row, data)
#                     else:
#                         empty.append(row)
#                     shown.append(row)
#             
#             print("\n# test whether only right columns are shown:")
#             print("shown_columns =", shown)
#             print("for col in shown_columns:")
#             print("\tself.assertEqual(view.table.isIndexHidden(model.index(col, 1)), False)")
#             print("for col in range({}):".format(model.rowCount()))
#             print("\tif col not in shown_columns:")
#             print("\t\tself.assertEqual(view.table.isIndexHidden(model.index(col, 1)), True)")
#             
#             print("\n# test content and headers of columns:")
#             print(assert_text)
#             
#             if empty:
#                 print("for col in {}:".format(empty))
#                 print('\tself.assertEqual(model.data(model.index(col, 0), Qt.DisplayRole), "")')
                
            # test whether only right columns are shown:
            shown_columns = [3, 38, 39, 40, 41, 42, 43, 44, 45]
            for col in shown_columns:
                self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), False)
            for col in range(51):
                if col not in shown_columns:
                    self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), True)
            
            # test content and headers of columns:
            self.assertEqual(model.headerData(3, Qt.Vertical, Qt.DisplayRole), "ENA Project ID")
            if not skip_other_tests:
                self.assertEqual(model.data(model.index(3, 0), Qt.DisplayRole), project_accession)
            self.assertEqual(model.headerData(38, Qt.Vertical, Qt.DisplayRole), "ENA Submission ID")
            #TODO: carry timestamp_sent, analysis acc nr and submission acc nr from ENA submission and check for these
#             if not skip_other_tests:
#                 self.assertEqual(model.data(model.index(38, 0), Qt.DisplayRole), "{}_{}".format(project_accession, ena_timestamp))
            self.assertEqual(model.headerData(39, Qt.Vertical, Qt.DisplayRole), "Timestamp sent")
#             if not skip_other_tests:
#                 self.assertEqual(model.data(model.index(39, 0), Qt.DisplayRole), ena_timestamp)
            self.assertEqual(model.headerData(40, Qt.Vertical, Qt.DisplayRole), "Timestamp confirmed") # not testing content
            self.assertEqual(model.headerData(41, Qt.Vertical, Qt.DisplayRole), "Analysis accession nr")
#             self.assertEqual(model.data(model.index(41, 0), Qt.DisplayRole), "ERZ678736")
            self.assertEqual(model.headerData(42, Qt.Vertical, Qt.DisplayRole), "Submission accession nr")
#             self.assertEqual(model.data(model.index(42, 0), Qt.DisplayRole), "ERA1561070")
            self.assertEqual(model.headerData(43, Qt.Vertical, Qt.DisplayRole), "Submission successful?")
            self.assertEqual(model.data(model.index(43, 0), Qt.DisplayRole), "yes")
            self.assertEqual(model.headerData(44, Qt.Vertical, Qt.DisplayRole), "ENA Acception Date")
            self.assertEqual(model.data(model.index(44, 0), Qt.DisplayRole), "2018-07-10")
            self.assertEqual(model.headerData(45, Qt.Vertical, Qt.DisplayRole), "ENA Accession Nr")
            self.assertEqual(model.data(model.index(45, 0), Qt.DisplayRole), "LT986596")

        def test_tab6_IPD(self):
            """tests whether content of SampleView widget 'Details about Allele' tab 'ENA submission' is correct
            """
            view = widget.tabs[5]
            model = view.flipped_model
#             # get reference data:
#             shown = []
#             empty = []
#             
#             assert_text = ""
#             for row in range(model.rowCount()):
#                 if not view.table.isIndexHidden(model.index(row, 1)):
#                     assert_text += 'self.assertEqual(model.headerData({}, Qt.Vertical, Qt.DisplayRole), "{}")\n'.format(row, model.headerData(row, Qt.Vertical, Qt.DisplayRole))
#                     data = model.data(model.index(row, 0), Qt.DisplayRole)
#                     if data:
#                         assert_text += 'self.assertEqual(model.data(model.index({}, 0), Qt.DisplayRole), "{}")\n'.format(row, data)
#                     else:
#                         empty.append(row)
#                     shown.append(row)
#             
#             print("\n# test whether only right columns are shown:")
#             print("shown_columns =", shown)
#             print("for col in shown_columns:")
#             print("\tself.assertEqual(view.table.isIndexHidden(model.index(col, 1)), False)")
#             print("for col in range({}):".format(model.rowCount()))
#             print("\tif col not in shown_columns:")
#             print("\t\tself.assertEqual(view.table.isIndexHidden(model.index(col, 1)), True)")
#             
#             print("\n# test content and headers of columns:")
#             print(assert_text)
#             
#             if empty:
#                 print("for col in {}:".format(empty))
#                 print('\tself.assertEqual(model.data(model.index(col, 0), Qt.DisplayRole), "")')

            # test whether only right columns are shown:
            shown_columns = [39, 40, 41, 42, 43, 44, 45, 46]
            for col in shown_columns:
                self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), False)
            for col in range(47):
                if col not in shown_columns:
                    self.assertEqual(view.table.isIndexHidden(model.index(col, 1)), True)
            
            # test content and headers of columns:
            self.assertEqual(model.headerData(39, Qt.Vertical, Qt.DisplayRole), "IPD Submission ID")
            #TODO: carry timestamp_sent, analysis acc nr and submission acc nr from ENA submission and check for these
#             if not skip_other_tests:
#                 self.assertEqual(model.data(model.index(39, 0), Qt.DisplayRole), "IPD_{}".format(IPD_timestamp))
            self.assertEqual(model.headerData(40, Qt.Vertical, Qt.DisplayRole), "Timestamp Data Ready")
#             if not skip_other_tests:
#                 self.assertEqual(model.data(model.index(40, 0), Qt.DisplayRole), "{}-{}-{}".format(IPD_timestamp[:4], IPD_timestamp[4:6], IPD_timestamp[6:8]))
            self.assertEqual(model.headerData(41, Qt.Vertical, Qt.DisplayRole), "Timestamp Confirmed") # not testing content
            self.assertEqual(model.headerData(42, Qt.Vertical, Qt.DisplayRole), "Data generated successfully?")
            self.assertEqual(model.data(model.index(42, 0), Qt.DisplayRole), "yes")
            self.assertEqual(model.headerData(43, Qt.Vertical, Qt.DisplayRole), "IPD Submission Nr")
            self.assertEqual(model.data(model.index(43, 0), Qt.DisplayRole), "DKMS1000{}".format(samples_dic["sample_1"]["submission_id"]))
            self.assertEqual(model.headerData(44, Qt.Vertical, Qt.DisplayRole), "HWS Submission Nr")
            self.assertEqual(model.headerData(45, Qt.Vertical, Qt.DisplayRole), "IPD Acception Date")
            self.assertEqual(model.headerData(46, Qt.Vertical, Qt.DisplayRole), "IPD_RELEASE")
            
            for col in [41, 44, 45, 46]:
                value = model.data(model.index(col, 0), Qt.DisplayRole)
                msg = "Col {} ({}) should be empty, but contains {}!".format(col, 
                                                                             model.headerData(col, Qt.Vertical, Qt.DisplayRole),
                                                                             value)
                self.assertEqual(value, "", msg)
        
        def test_tab7_history(self):
            """tests whether content of SampleView widget 'Details about Allele' tab 'allele history' is correct
            """
            view = widget.tabs[6]
            model = view.flipped_model
            curr_date = today.strftime("%Y-%m-%d")
            
            # test content and headers of columns:
            self.assertEqual(model.headerData(0, Qt.Vertical, Qt.DisplayRole), "Original genotyping")
            self.assertEqual(model.headerData(1, Qt.Vertical, Qt.DisplayRole), "Novel allele detection")
            self.assertEqual(model.headerData(2, Qt.Vertical, Qt.DisplayRole), "New genotyping")
            self.assertEqual(model.headerData(3, Qt.Vertical, Qt.DisplayRole), "Upload of sequence")
            self.assertEqual(model.data(model.index(3, 0), Qt.DisplayRole), curr_date)
            self.assertEqual(model.headerData(4, Qt.Vertical, Qt.DisplayRole), "Submitted to ENA")
            self.assertEqual(model.data(model.index(4, 0), Qt.DisplayRole), curr_date)
            self.assertEqual(model.headerData(5, Qt.Vertical, Qt.DisplayRole), "Accepted by ENA")
#             self.assertEqual(model.data(model.index(5, 0), Qt.DisplayRole), "2018-07-10") # compare to get_file_creation_date(ENA_reply_file), if necessary
            self.assertEqual(model.headerData(6, Qt.Vertical, Qt.DisplayRole), "Submitted to IPD")
            self.assertEqual(model.data(model.index(6, 0), Qt.DisplayRole), curr_date)
            self.assertEqual(model.headerData(7, Qt.Vertical, Qt.DisplayRole), "Accepted by IPD")
            self.assertEqual(model.data(model.index(7, 0), Qt.DisplayRole), "")
            
        test_tab1_general(self)
        test_tab2_typing_old(self)
        test_tab3_lab(self)
        test_tab4_typing_new(self)
        test_tab5_ena(self)
        test_tab6_IPD(self) # #TODO: re-enable as soon as final state of IPD tab is decided upon
        test_tab7_history(self)
        
        
class Test_Make_IMGT_Files_py(unittest.TestCase):
    """ 
    Test Make_IMGT_Files in typeloader_core
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_Make_IMGT_Files because skip_other_tests is set to True")
        else:
            self.data_dir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], samples_dic["sample_3"]["input_dir_origin"])
            self.start_num = samples_dic["sample_3"]["submission_id"]
            self.IPD_filename = "DKMS1000" + self.start_num
            self.samples = [(samples_dic["sample_3"]["id_int"], samples_dic["sample_3"]["local_name"], self.IPD_filename)] 
            self.file_dic = {samples_dic["sample_3"]["local_name"] : {"blast_xml" : samples_dic["sample_3"]["blast_file_name"],
                                                                     "ena_file" : samples_dic["sample_3"]["ena_file_name"]}}
            self.ENA_id_map, self.ENA_gene_map = MIF.parse_email(os.path.join(self.data_dir, samples_dic["sample_3"]["curr_ipd_ena_acc_file"]))
            self.pretypings = os.path.join(self.data_dir, samples_dic["sample_3"]["curr_ipd_befund_file"])
            self.curr_time = time.strftime("%Y%m%d%H%M%S")
            self.subm_id = "IPD_{}".format(self.curr_time)
            self.allele_dic = {samples_dic["sample_3"]["local_name"]: 
                               TargetAllele(gene=samples_dic["sample_3"]["gene"], 
                                            target_allele=samples_dic["sample_3"]["target_allele"], 
                                            partner_allele=samples_dic["sample_3"]["partner_allele"])
                               }
            
    @classmethod
    def tearDownClass(self):
        pass
            
    def test_confirmation_file(self):    
        """
        write_imgt_files --> make_imgt_data function is included
        test it with sample 3 --> confirmation file
        delete the ipd and zip file
        """
        MIF.write_imgt_files(self.data_dir, self.samples, self.file_dic, self.allele_dic, self.ENA_id_map, 
                             self.ENA_gene_map, self.pretypings, self.subm_id, 
                             self.data_dir, curr_settings, log)
        
        self.ipd_submission_file = os.path.join(self.data_dir, self.IPD_filename + "_confirmation.txt")
        self.ipd_submission_zipfile = os.path.join(self.data_dir, self.subm_id + ".zip")
        
        log.debug(self.ipd_submission_file)
        log.debug(self.ipd_submission_zipfile)
        self.assertTrue(os.path.exists(self.ipd_submission_file))
        self.assertTrue(os.path.exists(self.ipd_submission_zipfile))        
        
        is_confirmation = False
        if 'CC   Confirmation' in open(self.ipd_submission_file).read():
            is_confirmation = True
            
        self.assertTrue(is_confirmation)
        
        os.remove(self.ipd_submission_file)
        os.remove(self.ipd_submission_zipfile)
        

class Test_EMBL_functions(unittest.TestCase):
    """ 
    Test EMBL functions
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_EMBL_functions because skip_other_tests is set to True")
        else:
            self.title = "some random title"
            self.desc = "some random description"
            self.alias = "some random alias"
            self.xml_filename = "random.xml"
            self.accession = "PR_random"
            self.concat_FF_zip = "random.txt.zip" 
            self.checksum = "r1a2n3d4o5m6"
            self.project_schema = "project" # do not change
            self.analysis_schema = "analysis" # do not change        
            self.filetype = "flatfile" # do not change
            self.checksum_method = "MD5" # do not change
        
    @classmethod
    def tearDownClass(self):
        pass
    
    def test_func_generate_project_xml(self):
        """
        Test if project_xml_object has correct entries
        """
        project_dict = {}
        xml_stuff = EF.generate_project_xml(self.title, self.desc, self.alias, center_name)

        for item in xml_stuff.iter('PROJECT'):
            project_dict['alias'] = item.attrib['alias']
            project_dict['center_name'] = item.attrib['center_name']

        for item in xml_stuff.iter('TITLE'):
            project_dict['title'] = item.text
        
        for item in xml_stuff.iter('DESCRIPTION'):    
            project_dict['desc'] = item.text
        
        self.assertEqual(project_dict['alias'], self.alias)
        self.assertEqual(project_dict['center_name'], center_name)
        self.assertEqual(project_dict['title'], self.title)
        self.assertEqual(project_dict['desc'], self.desc)

    def test_func_generate_submission_project_xml(self):
        """
        Test if project_submission_xml_object has correct entries
        """
        project_submission_dict = {}
        xml_stuff = EF.generate_submission_project_xml(
            self.alias, 
            center_name, 
            self.xml_filename)
        
        for item in xml_stuff.iter('SUBMISSION'):
            project_submission_dict['alias'] = item.attrib['alias']
            project_submission_dict['center_name'] = item.attrib['center_name']

        for item in xml_stuff.iter('ADD'):
            project_submission_dict['schema'] = item.attrib['schema']
            project_submission_dict['source'] = item.attrib['source'] 
        
        self.assertEqual(project_submission_dict['alias'], self.alias)
        self.assertEqual(project_submission_dict['center_name'], center_name)
        self.assertEqual(project_submission_dict['schema'], self.project_schema)
        self.assertEqual(project_submission_dict['source'], self.xml_filename)        
        
    def test_func_generate_analysis_xml(self):
        """
        Test if analysis_xml_object has correct entries
        """
        analysis_dict = {}        
        xml_stuff = EF.generate_analysis_xml(
            self.title, 
            self.desc, 
            self.alias, 
            self.accession, 
            center_name, 
            self.concat_FF_zip, 
            self.checksum
        )
        
        for item in xml_stuff.iter('ANALYSIS'):
            analysis_dict['alias'] = item.attrib['alias']
            analysis_dict['center_name'] = item.attrib['center_name']        
        
        for item in xml_stuff.iter('STUDY_REF'):
            analysis_dict['accession'] = item.attrib['accession']
            
        for item in xml_stuff.iter('FILE'):
            analysis_dict['checksum'] = item.attrib['checksum']        
            analysis_dict['checksum_method'] = item.attrib['checksum_method']        
            analysis_dict['filename'] = item.attrib['filename']        
            analysis_dict['filetype'] = item.attrib['filetype']
        
        self.assertEqual(analysis_dict['alias'], self.alias)
        self.assertEqual(analysis_dict['center_name'], center_name)
        self.assertEqual(analysis_dict['accession'], self.accession)
        self.assertEqual(analysis_dict['checksum'], self.checksum)
        self.assertEqual(analysis_dict['checksum_method'], self.checksum_method)
        self.assertEqual(analysis_dict['filename'], self.concat_FF_zip)
        self.assertEqual(analysis_dict['filetype'], self.filetype)

    def test_func_generate_submission_ff_xml(self):
        """
        Test if analysis_submission_xml_object has correct entries
        """        
        analysis_submission_dict = {}        
        xml_stuff = EF.generate_submission_ff_xml(
            self.alias, 
            center_name, 
            self.xml_filename
        )
        
        for item in xml_stuff.iter('SUBMISSION'):
            analysis_submission_dict['alias'] = item.attrib['alias']
            analysis_submission_dict['center_name'] = item.attrib['center_name']        
            
        for item in xml_stuff.iter('ADD'):
            analysis_submission_dict['schema'] = item.attrib['schema']
            analysis_submission_dict['source'] = item.attrib['source']                    
        
        self.assertEqual(analysis_submission_dict['alias'], self.alias)
        self.assertEqual(analysis_submission_dict['center_name'], center_name)
        self.assertEqual(analysis_submission_dict['schema'], self.analysis_schema)
        self.assertEqual(analysis_submission_dict['source'], self.xml_filename)


class Test_rejection_short_UTR3(unittest.TestCase):
    """ 
    test if TypeLoader correctly rejects sequences with incomplete UTR3
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_rejection_short_UTR3 because skip_other_tests is set to True")
        else:
            self.mydir = os.path.join(curr_settings["login_dir"], "data_unittest", "rejection")
            self.testfile_fa = os.path.join(self.mydir, "UTR3_short.fa")
            self.missing_bp_fa = "53"
            self.testfile_xml = os.path.join(self.mydir, "UTR3_short.xml")
            self.missing_bp_xml = "1"
            self.project_name = project_name
        
    @classmethod
    def tearDownClass(self):
        pass
    
    def test_reject_fasta(self):
        """test if FASTA file with incomplete UTR3 is correctly rejected
        """
        myfile = self.testfile_fa
        missing_bp = self.missing_bp_fa
        
        # upload and parse file:
        results = typeloader_functions.upload_parse_sequence_file(myfile, curr_settings, log)
        (success_upload, _, filetype, _, 
         blastXmlFile, targetFamily, fasta_filename, allelesFilename, 
         header_data) = results
        
        self.assertTrue(success_upload) # uploading and parsing should work
        
        # try to create ENA file: (should fail)
        results2 = typeloader_functions.process_sequence_file(self.project_name, 
                                                filetype, blastXmlFile, targetFamily, 
                                                fasta_filename, allelesFilename, header_data, 
                                                curr_settings, log)
        
        (success, err_type, msg) = results2
        err = "File {} should have been rejected!".format(myfile)
        self.assertFalse(success, err)
        self.assertEqual(err_type, 'Incomplete sequence', "Should have thrown an 'Incomplete sequence' error")
        ref_error = errors.IncompleteSequenceError(missing_bp)
        self.assertEqual(msg, ref_error.msg)
        
        
    def test_reject_XML(self):
        """test if XML file with incomplete UTR3 is correctly rejected
        """
        myfile = self.testfile_xml
        missing_bp = self.missing_bp_xml
        
        # upload and parse file:
        results = typeloader_functions.upload_parse_sequence_file(myfile, curr_settings, log)
        (success_upload, _, filetype, _, 
         blastXmlFile, targetFamily, fasta_filename, allelesFilename, 
         header_data) = results
        
        self.assertTrue(success_upload) # uploading and parsing should work
        
        # try to create ENA file: (should fail)
        results2 = typeloader_functions.process_sequence_file(self.project_name, 
                                                filetype, blastXmlFile, targetFamily, 
                                                fasta_filename, allelesFilename, header_data, 
                                                curr_settings, log)
        
        (success, err_type, msg) = results2
        err = "File {} should have been rejected!".format(myfile)
        self.assertFalse(success, err)
        self.assertEqual(err_type, 'Incomplete sequence', "Should have thrown an 'Incomplete sequence' error")
        ref_error = errors.IncompleteSequenceError(missing_bp)
        self.assertEqual(msg, ref_error.msg)
        

class Test_null_alleles(unittest.TestCase):
    """ 
    test if TypeLoader correctly annotates null alleles
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_null_alleles because skip_other_tests is set to True")
        else:
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "null_allele")
            self.file_dic = {"C_no_na" : ("DKMS-LSL-C-2910.fa", "DKMS-LSL-C-2910.ena.txt"),
                                          "C_na" : ("DKMS-LSL-C-2952.fa", "DKMS-LSL-C-2952.ena.txt"),                               
                                          "DBP1_no_na" : ("DKMS-LSL-DPB1-2887.pacbio.minimap.fa", "DKMS-LSL-DPB1-2887.ena.txt"),
                                          "KIR2DL5_no_na" : ("DKMS-LSL-KIR2DL5-39.fa", "DKMS-LSL-KIR2DL5-39.ena.txt"),
                                          "KIR2DL5_na" : ("DKMS-LSL-KIR2DL5-39_corrupt.fa", "DKMS-LSL-KIR2DL5-39_corrupt.ena.txt"),                               
                                          "KIR3DP1_na" : ("DKMS-LSL-KIR3DP1-15.fa", "DKMS-LSL-KIR3DP1-15.ena.txt")
                                          }            
        
    @classmethod
    def tearDownClass(self):
        pass
    
    def test_null_allele_algorithm(self):
        """test if null allele algorithm works
        """        
        custom_seq_not_null = "ACCGTTTGGTGGTGA"
        custom_seq_stop_codon_1 = "ACCGTTTGGTAATGA" # TAA
        custom_seq_stop_codon_2 = "ACCGTTTGGTGATGA" # TGA
        custom_seq_stop_codon_3 = "ACCGTTTGGTAGTGA" # TAG
        custom_seq_not_divisable = "ACCGTTTGGTGTGA" # 14 bases
        
        self.assertFalse(BME.is_null_allele(custom_seq_not_null, {"cds": {1: (1, 15)}})[0])
        self.assertTrue(BME.is_null_allele(custom_seq_stop_codon_1, {"cds": {1: (1, 15)}})[0])
        self.assertTrue(BME.is_null_allele(custom_seq_stop_codon_2, {"cds": {1: (1, 15)}})[0])
        self.assertTrue(BME.is_null_allele(custom_seq_stop_codon_3, {"cds": {1: (1, 15)}})[0])
        self.assertTrue(BME.is_null_allele(custom_seq_not_divisable, {"cds": {1: (1, 14)}})[0])
    
    def test_fasta_null_allele(self):
        """test different fasta files, 
        """
        
        for key, value in self.file_dic.items():
            
            raw_path = os.path.join(self.mydir, value[0])
            reference_path = os.path.join(self.mydir, value[1])
            
            log.info("Processing:" + key)
            log.info("raw_path:" + raw_path)
            log.info("reference_path:" + reference_path)
                             
            results = typeloader_functions.upload_parse_sequence_file(raw_path, curr_settings, log)
            success_upload, sample_name, filetype, temp_raw_file, blastXmlFile, targetFamily, fasta_filename, allelesFilename, header_data = results
            self.assertTrue(success_upload, "Sequence file was not uploaded successfully")
            
            success, myalleles, ENA_text = typeloader_functions.process_sequence_file("PROJECT_NAME", filetype, blastXmlFile, targetFamily, fasta_filename, allelesFilename, header_data, curr_settings, log)        
            self.assertTrue(success, "Sequence file was not processed successfully")
            result = compare_2_files(reference_path = reference_path, query_var = ENA_text)
            self.assertEqual(len(result["added_sings"]), 0)
            self.assertEqual(len(result["deleted_sings"]), 0)              
 

class Test_multiple_novel_alleles_part1(unittest.TestCase):
    """ 
    test if TypeLoader correctly handles multiple novel alleles in one locus
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_null_alleles because skip_other_tests is set to True")
        else:
            self.project_name = project_name
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "multiple")
            self.sample_id_int = "ID000010"
            self.sample_id_ext = "test"
            self.local_name = 'DKMS-LSL_ID000010_3DP1_1'
            self.input_file = os.path.join(self.mydir, "ID000010.fa")
            self.pretypings = os.path.join(self.mydir, "Befunde.csv")
            self.ena_response_file = os.path.join(self.mydir, "ENA_reply")
              
            log.info("Check if allele already uploaded...")
            query = "select sample_id_int from alleles where local_name = '{}'".format(self.local_name)
            _, data = execute_db_query(query, 1, log, "Allele already present?", 
                                         "Checking ALLELES for test allele {}".format(self.sample_id_int), 
                                         "Can't assess ALLELES whether {} is contained".format(self.sample_id_int), 
                                         "ALLELES")
            
            if data: # set partner allele back to None
                query2 = "update alleles set partner_allele = Null where local_name = '{}'".format(self.local_name)
                _, data = execute_db_query(query2, 0, log, "Return partner allele to empty", 
                                         "Resetting ALLELES for test allele {}".format(self.sample_id_int), 
                                         "Can't reset ALLELES.partner_allele for {}".format(self.sample_id_int), 
                                         "ALLELES")
            else:
                log.info("Uploading allele...")
                self.form1 = ALLELE.NewAlleleForm(log, mydb, self.project_name, curr_settings, None, self.sample_id_int, samples_dic["sample_1"]["id_ext"])
                self.form1.file_widget.field.setText(os.path.join(self.mydir, "ID000010.fa"))
          
                self.form1.upload_btn.setEnabled(True)
                self.form1.upload_btn.click()
          
                self.form1.save_btn.click()
                 
            # prepare IPDSubmissionForm:
            self.form = IPD.IPDSubmissionForm(log, mydb, self.project_name, curr_settings, parent = None)
              
    @classmethod
    def tearDownClass(self):
        pass
    
    @patch("GUI_forms_submission_IPD.BothAllelesNovelDialog")
    def test_multi_allelesDialogCreated(self, mock_dialog):
        """test if target allele with multiple novel alleles according to pretyping 
        but unfitting target_allele/partner_allele correctly trigger the BothAllelesNovelDialog
        """
        # choose project:
        log.info("Choosing project...")
        self.form.proj_widget.field.setText(self.project_name)
        self.form.ok_btn1.check_ready()
        self.form.ok_btn1.click()
        # add files:
        log.info("Choosing files...")
        self.form.ENA_file_widget.field.setText(self.ena_response_file)
        self.form.befund_widget.field.setText(self.pretypings)
        self.form.ok_btn2.check_ready()
        self.form.ok_btn2.click()
        # select alleles:
        log.info("Selecting alleles...")
        if not self.form.project_files.check_dic[0].isChecked():
            self.form.project_files.check_dic[0].click() # if this fails, the alleles in the ENA reply file are probably not recognized correctly
        self.form.submit_btn.check_ready()
        self.form.submit_btn.click()
          
        # check if multi_dialog is raised:
        log.info("Checking that BothAllelesNovelDialog is created...")
        mock_dialog.assert_called_once()
        log.info("=> fine")
        self.form.ok_btn.click()
        
    def test_multi_alleles_dialog(self):
        """test if BothAllelesNovelDialog input produces the intended behavior
        """
        # open MultiAlleleDialog with input data:
        problem_dic = {'DKMS-LSL_ID000010_3DP1_1': ['ID000010', 'DKMS-LSL_ID000010_3DP1_1', 
                                                    TargetAllele(gene='KIR3DP1', target_allele='KIR3DP1*0030102:new', 
                                                                 partner_allele=''), 
                                                    ['003new', '004new', '001']]}
        mydialog = IPD.BothAllelesNovelDialog(problem_dic, curr_settings, log)
        mybox = mydialog.choice_boxes[self.local_name]
        mybox.options[0].click()
        mydialog.submit_btn.click()
        
        # check correct results are now in the database:
        query = "select partner_allele from alleles where local_name = '{}'".format(self.local_name)
        _, data = execute_db_query(query, 1, log, "Partner_allele updated?", 
                                     "Checking ALLELES for partner_allele  of {}".format(self.sample_id_int), 
                                     "Can't assess ALLELES whether partner_allele for {} has been updated".format(self.sample_id_int), 
                                     "ALLELES")
        partner_allele = data[0][0]
        log.info("Partner allele updated correctly?")
        self.assertEqual(partner_allele, "KIR3DP1*004new and 001", "Partner_allele was not updated correctly!")
        log.info("\t=> yes")


class Test_multiple_novel_alleles_part2(unittest.TestCase):
    """test if BothAllelesNovelDialog is not called if target_allele and partner_allele are fine
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_null_alleles because skip_other_tests is set to True")
        else:
            self.project_name = project_name
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "multiple")
            self.sample_id_int = "ID000010"
            self.sample_id_ext = "test"
            self.local_name = 'DKMS-LSL_ID000010_3DP1_1'
            self.input_file = os.path.join(self.mydir, "ID000010.fa")
            self.pretypings = os.path.join(self.mydir, "Befunde.csv")
            self.ena_response_file = os.path.join(self.mydir, "ENA_reply")
            # prepare IPDSubmissionForm:
            self.form = IPD.IPDSubmissionForm(log, mydb, self.project_name, curr_settings, parent = None)
              
    @classmethod
    def tearDownClass(self):
        pass
 
    @patch("GUI_forms_submission_IPD.BothAllelesNovelDialog")
    def test_multi_allelesDialogNotCreated(self, mock_dialog):
        """test if target allele with multiple novel alleles according to pretyping 
        but unfitting target_allele/partner_allele correctly trigger the BothAllelesNovelDialog
        """
        self.form = IPD.IPDSubmissionForm(log, mydb, self.project_name, curr_settings, parent = None)
        # choose project:
        log.info("Choosing project...")
        self.form.proj_widget.field.setText(self.project_name)
        self.form.ok_btn1.check_ready()
        self.form.ok_btn1.click()
        # add files:
        log.info("Choosing files...")
        self.form.ENA_file_widget.field.setText(self.ena_response_file)
        self.form.befund_widget.field.setText(self.pretypings)
        self.form.ok_btn2.check_ready()
        self.form.ok_btn2.click()
        # select alleles:
        log.info("Selecting alleles...")
        if not self.form.project_files.check_dic[0].isChecked():
            self.form.project_files.check_dic[0].click() # if this fails, the alleles in the ENA reply file are probably not recognized correctly
        self.form.submit_btn.check_ready()
        self.form.submit_btn.click()
          
        # check if multi_dialog is raised:
        log.info("Checking that BothAllelesNovelDialog is NOT created...")
        mock_dialog.assert_not_called()
        log.info("=> fine")
 
 
class Test_pretyping_valid(unittest.TestCase):
    """ 
    test if TypeLoader correctly handles valid and invalid pretypings
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_pretyping_valid because skip_other_tests is set to True")
        else:
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "pretyping_check")
             
            # read samples from csv:
            log.info("Reading files...")
            log.debug("Reading samples.csv...")
            SampleObject = namedtuple("SampleObject", """name description closest_allele gene 
                                                        target_allele partner_allele 
                                                        target_family diff_text final_result error_exp""")
            self.samples = {}
            with open(os.path.join(self.mydir, "samples.csv")) as f:
                data = csv.reader(f, delimiter=",")
                for i, row in enumerate(data):
                    if i != 0:
                        if row:
                            s = SampleObject(name = row[0],
                                             description = row[1],
                                             closest_allele = row[2],
                                             gene = row[3], 
                                             target_allele = row[4],
                                             partner_allele = row[5],
                                             target_family = row[6],
                                             diff_text = row[7],
                                             final_result = row[8],
                                             error_exp = row[9])
                            self.samples[s.name] = s
             
            # read pretypings from csv:
            log.debug("Reading pretypings.csv...")
            (self.pretypings, self.topics) = MIF.getPatientBefund(os.path.join(self.mydir, "pretypings.csv"))
 
            # pepare error_dic:
            from typeloader_core import errors
            self.error_dic = {"Cannot tell which novel allele from pretyping this is" : errors.BothAllelesNovelError,
                              "no allele marked as new in pretyping" : errors.InvalidPretypingError,
                              "assigned allele name not found in pretyping" : errors.InvalidPretypingError,
                              "POS is not acceptable pretyping for a target locus": errors.InvalidPretypingError,
                              "pretyping for HLA-B missing" : errors.InvalidPretypingError}
                         
    @classmethod
    def tearDownClass(self):
        pass
      
    def test_pretypings(self):
        """test if pretypings and multiple alleles are handled correctly
        """
        for name in self.samples:
            s = self.samples[name]
            log.debug("testing {}: ({})".format(s.name, s.description))
            target_allele = TargetAllele(gene = s.gene, target_allele = s.target_allele, partner_allele = s.partner_allele)
             
            if s.error_exp: # sad path testing (are correct errors raised?)
                myerror = self.error_dic[s.error_exp]
                # make sure right error type is raised:
                self.assertRaises(myerror, ITG.make_befund_text, self.pretypings[name], s.closest_allele, target_allele, 
                                           {'gene': ['HLA', 'KIR'], 'targetFamily': s.target_family},
                                           s.diff_text, log)
                # make sure error text is correct if an error is raised:
                try:
                    ITG.make_befund_text(self.pretypings[name], s.closest_allele, target_allele, 
                                           {'gene': ['HLA', 'KIR'], 'targetFamily': s.target_family},
                                           s.diff_text, log)
                except Exception as E:
                    self.assertEqual(E.problem, s.error_exp)
                 
            else: # happy path testing (these should pass without errors)
                txt = ITG.make_befund_text(self.pretypings[name], s.closest_allele, target_allele, 
                                           {'gene': ['HLA', 'KIR'], 'targetFamily': s.target_family},
                                           s.diff_text, log)
                for line in txt.split("\n"):
                    if s.gene in line:
                        self.assertEqual(line, "FT                  /{}".format(s.final_result), 
                                         "Error in {}: {}".format(s.name, s.description)) # check result correct
            
       
class Test_Clean_Stuff(unittest.TestCase):
    """ 
    Remove all directories and files written by  all unit tests
    """
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping final cleanup because skip_other_tests is set to True")  
        
    @classmethod
    def tearDownClass(self):
        typeloader_GUI.close_connection(log, mydb)
        
    def test_clean_everything(self):
        os.remove(os.path.join(mypath_inner, "typeloader_GUI.py"))
        
        if delete_all_stuff_at_the_end:
            delete_written_samples(True, "ALLELES", log)
            delete_written_samples(True, "FILES", log)
            delete_written_samples(True, "SAMPLES", log)
            delete_written_samples(True, "PROJECTS", log)
            delete_written_samples(True, "ENA_SUBMISSIONS", log)
            delete_written_samples(True, "IPD_SUBMISSIONS", log)
                
            shutil.rmtree(os.path.join(curr_settings["projects_dir"], project_name))
        

#===========================================================
# functions:

def execute_db_query(query, num_columns, log, main_log, success_log, fail_log, format_string = ""):
    """
    Executes statements on sqlite database
    returns success, data
    """
    log.info(main_log.format(format_string))
    return db_internal.execute_query(query, num_columns, log, 
                              success_log.format(format_string), 
                              fail_log.format(format_string))

def delete_written_samples(clear_every_row, table, log, column = "", value = ""):
    """
    Deletes sample entries made by Test_XXX
    Cleaning up step
    """
    log.info("Removing unittest samples in staging database...")
    if clear_every_row:
        query = "DELETE from {}".format(table)
    else:
        query = "DELETE from {} where {} = '{}';".format(table, column, value) 
    
    db_internal.execute_query(query, 0, log, 
                              "deleting samples with {} {} in table {}".format(column, value, table), 
                              "cant delete sample with {} {} in table {}".format(column, value, table))

def compare_2_files(query_path = "", reference_path = "", filetype = "", query_var = "", reference_var = ""):
    
    result = {}
    
    if query_path == "":
        query_text = query_var
    else:
        with open(query_path , 'r') as myfile: 
            query_text = myfile.read()
        
    if reference_path == "":
        reference_text = reference_var
    else:
        with open(reference_path , 'r') as myfile: reference_text = myfile.read()            
    
    if filetype == "IPD":
        # change the date in order to compare both ipd files
        now = today.strftime("%d/%m/%Y") 
        reference_text = re.sub('DT.*Submitted\)\n.*Release\)', 'DT   {} (Submitted)\nDT   {} (Release)'.format(now, now), reference_text)
        reference_text = reference_text.replace("{TL-VERSION}", __version__) # replace TL-version part of reference file with current version

    diffInstance = difflib.Differ()
    diffList = list(diffInstance.compare(query_text.strip(), reference_text.strip()))
    
    result["added_sings"] = ' '.join(x[2:] for x in diffList if x.startswith('+ '))
    result["deleted_sings"]= ' '.join(x[2:] for x in diffList if x.startswith('- '))    
    if result["added_sings"] != "" or result["deleted_sings"] != "":
        log.error("Differences found!")
        log.debug("New file: {}".format(query_path))
        log.debug("Reference file: {}".format(reference_path))
        log.warning("Differences (added): " + result["added_sings"])
        log.warning("Differences (deleted): " + result["deleted_sings"])
    return result

# functions to order Tests:
# from: https://gist.github.com/catb0t/304ececa6c55f6e3788d
# case_factory: gets all Test_Classes and tests
# suiteFacotry: orders the testcases

def suiteFactory(
        *testcases,
        testSorter   = None,
        suiteMaker   = unittest.makeSuite,
        newTestSuite = unittest.TestSuite
    ):
    
    """
    make a test suite from test cases, or generate test suites from test cases.
    *testcases     = TestCase subclasses to work on
    testSorter     = sort tests using this function over sorting by line number
    suiteMaker     = should quack like unittest.makeSuite.
    newTestSuite   = should quack like unittest.TestSuite.
    """

    if testSorter is None:
        ln         = lambda f:    getattr(tc, f).__code__.co_firstlineno
        testSorter = lambda a, b: ln(a) - ln(b)

    test_suite = newTestSuite()
    for tc in testcases:
        test_suite.addTest(suiteMaker(tc, sortUsing=testSorter))

    return test_suite


def caseFactory(
        scope        = globals().copy(),
        caseSorter   = lambda f: __import__("inspect").findsource(f)[1],
        caseSuperCls = unittest.TestCase,
        caseMatches  = __import__("re").compile("^Test")
    ):
    
    """
    get TestCase-y subclasses from frame "scope", filtering name and attribs
    scope        = iterable to use for a frame; preferably a hashable (dictionary).
    caseMatches  = regex to match function names against; blank matches every TestCase subclass
    caseSuperCls = superclass of test cases; unittest.TestCase by default
    caseSorter   = sort test cases using this function over sorting by line number
    """
    return sorted(
        [
            scope[obj] for obj in scope
                if re.match(caseMatches, obj)
                and issubclass(scope[obj], caseSuperCls)
        ],
        key=caseSorter
    )

def log_uncaught_exceptions(cls, exception, tb):
    """reimplementation of sys.excepthook;
    catches uncaught exceptions, logs them and exits the app
    """
    import traceback
    from PyQt5.QtCore import QCoreApplication
    log.critical('{0}: {1}'.format(cls, exception))
    log.exception(msg = "Uncaught Exception", exc_info = (cls, exception, tb))
    #TODO: (future) maybe find a way to display the traceback only once, both in console and logfile?
    sys.__excepthook__(cls, exception, traceback)
    QCoreApplication.exit(1)
pass
#===========================================================
# main:

if __name__ == "__main__":
    sys.excepthook = log_uncaught_exceptions 
    cases = suiteFactory(*caseFactory())
    runner = unittest.TextTestRunner(verbosity=2, failfast=True)
    runner.run(cases)    
    

