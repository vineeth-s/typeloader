#!/usr/bin/env Python3
# -*- coding: cp1252 -*-
'''
test_typeloader_GUI.py

unit tests for typeloader_GUI

@author: Bianca Schoene
'''
import gzip
import unittest
from unittest.mock import patch
import os, sys, re, time, platform, datetime, csv
import difflib  # compare strings
import shutil
import copy
from pathlib import Path
from random import randint
from configparser import ConfigParser

mypath_inner = Path(__file__).parent.parent
mypath = mypath_inner.parent

from typeloader2 import general, db_internal, GUI_login
from xml.etree import ElementTree
from collections import namedtuple

# no .pyw import possible in linux
# deletion in Test_Clean_Stuff
shutil.copyfile(mypath_inner / "typeloader_GUI.pyw", mypath_inner / "typeloader_GUI.py")

from typeloader2 import typeloader_GUI
from typeloader2.typeloader_core import errors, EMBLfunctions as EF, make_imgt_files as MIF, backend_make_ena as BME, \
    imgt_text_generator as ITG, closestallele as CA, getAlleleSeqsAndBlast as GASB
from typeloader2 import GUI_forms_new_project as PROJECT
from typeloader2 import GUI_forms_new_allele as ALLELE
from typeloader2 import GUI_forms_new_allele_bulk as BULK
from typeloader2 import GUI_forms_submission_ENA as ENA
from typeloader2 import GUI_forms_submission_IPD as IPD
from typeloader2 import GUI_views_OVprojects, GUI_views_OValleles, GUI_views_project, GUI_views_sample, GUI_forms, GUI_download_files
from typeloader2 import GUI_mini_dialogs
from typeloader2 import typeloader_functions
from typeloader2.GUI_login import base_config_file, check_update_needed
from typeloader2 import db_external

from PyQt5.QtWidgets import (QApplication)
from PyQt5.QtCore import Qt, QTimer, QModelIndex

# ===========================================================
# test parameters:

delete_all_stuff_at_the_end = False  # deletes database entries and project directory
skip_other_tests = False  # set to True to skip all tests except the current WiP (out-comment it there in setUpClass)
skip_make_project = False  # set to True to skip initial cleanup and new project creation
project_name = ""  # this will be set in create project

samples_dic = {  # samples to test
    "sample_1": {"input": "input.fa",
                 "input_dir_origin": "KIR_2DL4",
                 "local_name": "DKMS-LSL_ID000001_2DL4_1",
                 "cell_line": "DKMS-LSL_ID000001",
                 "gene": "KIR2DL4",
                 "target_allele": "KIR2DL4*0080401:new",
                 "data_unittest_dir": "KIR_2DL4",
                 "curr_ena_file": "DKMS-LSL_ID000001_2DL4_1.ena.txt",
                 "curr_fasta_file": "DKMS-LSL_ID000001_2DL4_1.fa",
                 "curr_blast_file": "DKMS-LSL_ID000001_2DL4_1.blast.xml",
                 "curr_ipd_befund_file": "Befunde_neu_1.csv",
                 "ena_acc": "4ST3TSE6",
                 "id_int": "ID000001",
                 "id_ext": "DEDKM000001",
                 "submission_id": "1111",
                 "customer": "DKMS",
                 "provenance": "United Kingdom",
                 "collection_date": "2023"},
    "sample_2": {"input": "5597571 A.xml",
                 "input_dir_origin": "A_MM",
                 "local_name": "DKMS-LSL_ID14278154B_A_1",
                 "cell_line": "DKMS-LSL_ID14278154B",
                 "gene": "HLA-A",
                 "target_allele": "HLA-A*01:new",
                 "partner_allele": 'HLA-A*33:new',
                 "data_unittest_dir": "new_allele_xml",
                 "curr_ena_file": "DKMS-LSL_ID14278154B_A_1.ena.txt",
                 "curr_fasta_file": "DKMS-LSL_ID14278154B_A_1.fa",
                 "curr_blast_file": "DKMS-LSL_ID14278154B_A_1.blast.xml",
                 "curr_gendx_file": "DKMS-LSL_ID14278154B_A_1.xml",
                 "curr_ipd_befund_file": "Befunde_A_WDH.csv",
                 "id_int": "ID14278154B",
                 "id_ext": "1348480",
                 "submission_id": "2222",
                 "provenance": "missing: third party data",
                 "collection_date": "missing: third party data",
                 "customer": ""
                 },
    "sample_3": {"input_dir_origin": "confirmation_file",
                 "local_name": "DKMS-LSL_ID15390636_3DP1_1",
                 "cell_line": "DKMS-LSL_ID15390636",
                 "gene": "HLA-DPB1",
                 "curr_ipd_befund_file": "Befunde_3DP1_1.csv",
                 "blast_file_name": "DKMS-LSL_ID15390636_3DP1_1.blast.xml",
                 "ena_file_name": "DKMS-LSL_ID15390636_3DP1_1.ena.txt",
                 "ena_acc": "LT986495",
                 "id_int": "ID15390636",
                 "id_ext": "1370324_A",
                 "submission_id": "3333",
                 "gene": "KIR3DP1",
                 "target_allele": 'KIR3DL3*006',
                 "partner_allele": 'KIR3DL3*003'}
}

settings_both = {"reference_dir": "reference_data_unittest",
                 "data_unittest": "data_unittest",
                 "ipd_submissions": "IPD-submissions",
                 }

log = general.start_log(level="DEBUG", info_to_file="typeloader_tests.log")
__version__ = general.read_package_variable("__version__")
log.info(f"Running TypeLoader tests for version {__version__}...")
# assemble settings for testing:
cf = ConfigParser()
cf.read(os.path.join(mypath_inner, base_config_file))
curr_settings = GUI_login.get_settings("staging", log, cf)
for key in settings_both:
    curr_settings[key] = settings_both[key]
curr_settings["embl_submission"] = curr_settings["embl_submission_test"]
curr_settings["running_modus"] = "automatic tests"
typeloader_functions.update_curr_versions(curr_settings, log)

embl_test_server = curr_settings["embl_submission"]
center_name = curr_settings["xml_center_name"]
mydb = typeloader_GUI.create_connection(log, curr_settings["db_file"])

# ipd_counter_lock_file = os.path.join(curr_settings["root_path"], "_general", "ipd_nr.lock")
# counter_config_file = os.path.join(curr_settings["root_path"], "_general", "counter_config.ini")
# _, (IPD_counter, _) = MIF.get_IPD_counter(counter_config_file, ipd_counter_lock_file, curr_settings, log)

project_gene = "X"
project_pool = str(randint(1, 999999))
project_user = "Staging Account"
project_title = "This is an optional title information"
project_desc = "This is an optional description information"
project_accession = ""  ## this will be set in create project

KIR_VERSION = "290"  # the KIR version with which the reference KIR files were created
HLA_VERSION = "3412"  # this line is purely informational; the HLA version is not changed during any tests

app = QApplication(sys.argv)

today = datetime.datetime.now()

TargetAllele = namedtuple("TargetAllele", "gene target_allele partner_allele")


# ===========================================================
# test cases:

class Test_0_Clean_Stuff_initial(unittest.TestCase):
    """
    Remove all directories and files written by previous unit tests
    """

    @classmethod
    def setUpClass(self):
        if skip_make_project:
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


class TestResetReferenceManually(unittest.TestCase):
    """test whether manual reference resets work
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping reference reset test because skip_other_tests is set to True")

        self.form = GUI_mini_dialogs.ResetReferenceDialog(curr_settings, log, parent=None, testing=True)
        self.target = "KIR"
        self.target_version = "2.7.1"
        self.target_version_after_update = "271"
        self.reference_local_path = os.path.join(curr_settings["root_path"],
                                                 curr_settings["general_dir"],
                                                 curr_settings["reference_dir"])

    @classmethod
    def tearDownClass(self):
        pass

    def test01_version_handling(self):
        """test that version input is handled correctly
         """
        scenarios = [("", False, "Missing version", "Please insert a version number, e.g., 3.39.0!"),
                     ("39b", False, "Invalid version",
                      "Version should only contain digits and up to 2 dots, like '2.7.0' or '3.39'!"),
                     ("2.7.0", True, "270", None),
                     ("27", True, "270", None),
                     ("2.7.1", True, "271", None),
                     ("2.7", True, "270", None),
                     ("2-7", True, "270", None),
                     ("2_7", True, "270", None),
                     ]
        for (version_org, success, exp_version, exp_msg) in scenarios:
            version_ok, new_version, msg = self.form.handle_version_input(version_org)
            self.assertEqual(version_ok, success)
            self.assertEqual(new_version, exp_version)
            self.assertEqual(msg, exp_msg)

    def test02_confirmation_handling(self):
        """test that confirmation handling is correct
         """
        # TODO: check actual input
        proceed, msg = self.form.check_proceed(self.target)
        self.assertTrue(proceed)
        self.assertEqual(msg, None)

    def test03_happy_path(self):
        self.form.version_field.setText(self.target_version)
        self.form.btn_dic[self.target].click()
        self.assertEqual(self.form.updated, self.target_version_after_update)


class TestUpdateReference(unittest.TestCase):
    """test whether reference updates work
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping reference update test because skip_other_tests is set to True")

        self.form = GUI_mini_dialogs.RefreshReferenceDialog(curr_settings, log, None)
        self.target = "KIR"
        self.reference_local_path = os.path.join(curr_settings["root_path"],
                                                 curr_settings["general_dir"],
                                                 curr_settings["reference_dir"])

    @classmethod
    def tearDownClass(self):
        """reset database to the version with which the reference files were created"""
        success, err, update_msg = typeloader_functions.perform_reference_update(self.target,
                                                                                 self.reference_local_path,
                                                                                 curr_settings["blast_path"],
                                                                                 curr_settings["proxy"],
                                                                                 log,
                                                                                 version=KIR_VERSION)
        if not success:
            log.error(err)
            log.info(update_msg)

        self.assertTrue(success, "version setback not successful!")
        typeloader_functions.update_curr_versions(curr_settings, log)

    def test01_ref_update_kir_noerrors(self):
        """test that KIR is updated successfully;
        only tests that no errors occured
         """
        self.form.btn_dic[self.target].click()
        self.assertEqual(self.form.updated, [self.target])

    def test02_ref_update_kir_check_ok(self):
        """test that after the update, the .dat file checksum is equal to that in the reference repo
         """
        update_me = check_update_needed(self.reference_local_path, curr_settings["proxy"], log,
                                        skip_if_updated_today=False)
        self.assertTrue("KIR" not in update_me)


class Test_1_Create_Project(unittest.TestCase):
    """ create project
    """

    @classmethod
    def setUpClass(self):
        if skip_make_project:
            self.skipTest(self, "Skipping Create Test because skip_make_project is set to True")
        else:
            self.form = PROJECT.NewProjectForm(log, mydb, curr_settings)
            self.form.project_dir = ""  # not initially set in form

            # format: (field, default_value, whitespace_ok)
            self.myvalues = [(self.form.user_entry, project_user, True),
                             (self.form.gene_entry, project_gene, False),
                             (self.form.pool_entry, project_pool, False),
                             (self.form.title_entry, project_title, True),
                             (self.form.desc_entry, project_desc, True)]
            for (field, value, _) in self.myvalues:
                field.setText(value)

    @classmethod
    def tearDownClass(self):
        pass

    @unittest.skipIf(skip_other_tests, "skip_other_tests is set to True")
    def test_1_reject_invalid(self):
        """make sure invalid entries are rejected properly
        """
        invalid = ['"', "'", "_", "#", "%"]  # disallowed characters without spaces

        for (field, orig_value, whitespace_ok) in self.myvalues:
            if whitespace_ok:
                invalid_chars = invalid
            else:
                invalid_chars = invalid + [" "]  # spaces should be rejected, too

            # check that default value is ok:
            field.setText(orig_value)
            self.form.get_values()
            invalid_msg = self.form.check_all_fields_valid()
            self.assertFalse(invalid_msg)

            # check that invalid characters are rejected:
            for char in invalid_chars:
                mytext = orig_value + char
                field.setText(mytext)
                self.form.get_values()
                invalid_msg = self.form.check_all_fields_valid()
                self.assertTrue(invalid_msg)
            field.setText(orig_value)

    def test_2_create_project_success(self):
        """
        Defines a project name and creates a project on ENA Test Server
        """
        for (field, value, _) in self.myvalues:
            field.setText(value)

        self.form.project_btn.click()
        self.form.submit_btn.click()
        self.form.accession = self.form.acc_entry.text()

        ## set global var for further tests
        global project_name
        global project_accession
        project_name = self.form.project_name
        project_accession = self.form.accession

    @unittest.skipIf(skip_other_tests, "skip_other_tests is set to True")
    def test_3_parse_project_name(self):
        """
        parse project name
        """
        date = general.timestamp("%Y%m%d")
        new_project = "_".join([date, "SA", project_gene, project_pool])
        self.assertEqual(self.form.project_name, new_project)

    @unittest.skipIf(skip_other_tests, "skip_other_tests is set to True")
    def test_4_dir_exists(self):
        """
        If project dir is created?
        """
        self.assertTrue(os.path.exists(self.form.project_dir))

    @unittest.skipIf(skip_other_tests, "skip_other_tests is set to True")
    def test_5_projectfiles_exists(self):
        """
        If project files are created successfull
        """
        self.form.project_file = os.path.join(self.form.project_dir, self.form.project_name + ".xml")
        self.form.submission_file = os.path.join(self.form.project_dir, self.form.project_name + "_sub.xml")
        self.form.output_file = os.path.join(self.form.project_dir, self.form.project_name + "_output.xml")
        self.assertTrue(os.path.exists(self.form.project_file))
        self.assertTrue(os.path.exists(self.form.submission_file))
        self.assertTrue(os.path.exists(self.form.output_file))

    @unittest.skipIf(skip_other_tests, "skip_other_tests is set to True")
    def test_6_parse_project_xml(self):
        """
        Parse the written project XML file
        """
        xml_stuff = ElementTree.parse(self.form.project_filename)
        root = xml_stuff.getroot()

        self.assertEqual(root[0].attrib["alias"], self.form.project_name)
        self.assertEqual(root[0].attrib["center_name"], center_name)
        self.assertEqual(root[0][0].text, project_title)
        self.assertEqual(root[0][1].text, project_desc)

    @unittest.skipIf(skip_other_tests, "skip_other_tests is set to True")
    def test_7_parse_submission_xml(self):
        """
        Parse the written submission XML file
        """
        xml_stuff = ElementTree.parse(self.form.submission_file)
        root = xml_stuff.getroot()

        self.assertEqual(root.attrib["alias"], self.form.project_name + "_sub")
        self.assertEqual(root.attrib["center_name"], center_name)
        self.assertEqual(root[0][0][0].attrib["schema"], "project")
        self.assertEqual(root[0][0][0].attrib["source"], self.form.project_name + ".xml")

    @unittest.skipIf(skip_other_tests, "skip_other_tests is set to True")
    def test_8_parse_output_xml(self):
        """
        Parse the written output XML file
        """
        xml_stuff = ElementTree.parse(self.form.output_file)
        root = xml_stuff.getroot()

        self.assertEqual(root.attrib["submissionFile"], self.form.project_name + "_sub.xml")
        self.assertEqual(root.attrib["success"], "true")
        self.assertEqual(root[0].attrib["accession"], self.form.accession)
        self.assertEqual(root[1].attrib["alias"], self.form.project_name + "_sub")
        self.assertEqual(root[2][-1].text, "This submission is a TEST submission and will be discarded within 24 hours")


class Test_2_ProjectStatus(unittest.TestCase):
    """ check toggling of project status
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping ProjectStatus Test because skip_other_tests is set to True")
        else:
            self.project = project_name

    @classmethod
    def tearDownClass(self):
        pass

    def test1_get_status(self):
        """test if we can get the status of our new project correctly
        """
        project_open = GUI_forms.check_project_open(self.project, log)
        self.assertEqual(project_open, True)

    def test2_get_status_nonsense(self):
        """test if we get status 'Open' for a non-existing project
        (policy: when in doubt, treat as 'Open')
        """
        project_open = GUI_forms.check_project_open("bla", log)
        self.assertEqual(project_open, True)

    def test3_toggle_status_to_closed(self):
        """test if we can toggle the status of our project to 'Closed'
        """
        success, new_status, new_index = typeloader_functions.toggle_project_status(self.project, "Open", log)
        self.assertEqual(success, True)
        self.assertEqual(new_status, "Closed")
        self.assertEqual(new_index, 1)

    def test4_get_status_closed(self):
        """test if we can get the status correctly if project is closed
        """
        project_open = GUI_forms.check_project_open(self.project, log)
        self.assertEqual(project_open, False)

    def test5_toggle_status_to_open(self):
        """test if we can toggle the project status back to 'Open'
        """
        success, new_status, new_index = typeloader_functions.toggle_project_status(self.project, "Closed", log)
        self.assertEqual(success, True)
        self.assertEqual(new_status, "Open")
        self.assertEqual(new_index, 0)


class Test_Create_New_Allele(unittest.TestCase):
    """
    create new allele
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Create New Allele because skip_other_tests is set to True")
        else:
            self.project_name = project_name  # "20180710_SA_A_1292"

            self.new_sample_dir_path_1 = os.path.join(curr_settings["projects_dir"], self.project_name,
                                                      samples_dic["sample_1"]["id_int"])
            self.new_ena_file_path_1 = os.path.join(curr_settings["projects_dir"], self.project_name,
                                                    samples_dic["sample_1"]["id_int"],
                                                    samples_dic["sample_1"]["curr_ena_file"])
            self.new_fasta_file_path_1 = os.path.join(curr_settings["projects_dir"], self.project_name,
                                                      samples_dic["sample_1"]["id_int"],
                                                      samples_dic["sample_1"]["curr_fasta_file"])
            self.new_blast_file_path_1 = os.path.join(curr_settings["projects_dir"], self.project_name,
                                                      samples_dic["sample_1"]["id_int"],
                                                      samples_dic["sample_1"]["curr_blast_file"])
            self.reference_file_path_1 = os.path.join(curr_settings["login_dir"], "data_unittest",
                                                      samples_dic["sample_1"]["curr_ena_file"])

            self.new_sample_dir_path_2 = os.path.join(curr_settings["projects_dir"], self.project_name,
                                                      samples_dic["sample_2"]["id_int"])
            self.new_ena_file_path_2 = os.path.join(curr_settings["projects_dir"], self.project_name,
                                                    samples_dic["sample_2"]["id_int"],
                                                    samples_dic["sample_2"]["curr_ena_file"])
            self.new_fasta_file_path_2 = os.path.join(curr_settings["projects_dir"], self.project_name,
                                                      samples_dic["sample_2"]["id_int"],
                                                      samples_dic["sample_2"]["curr_fasta_file"])
            self.new_blast_file_path_2 = os.path.join(curr_settings["projects_dir"], self.project_name,
                                                      samples_dic["sample_2"]["id_int"],
                                                      samples_dic["sample_2"]["curr_blast_file"])
            self.reference_file_path_2 = os.path.join(curr_settings["login_dir"], "data_unittest",
                                                      samples_dic["sample_2"]["curr_ena_file"])

    @classmethod
    def tearDownClass(self):
        pass

    # @unittest.skip("skipping test_fasta_file")
    def test_fasta_file(self):
        """
        Ceate ENA flatfile from fasta
        """
        self.form = ALLELE.NewAlleleForm(log, mydb, self.project_name, curr_settings, None,
                                         samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["id_ext"])
        input_file = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                  samples_dic["sample_1"]["data_unittest_dir"],
                                  samples_dic["sample_1"]["input"])
        self.form.file_widget.field.setText(input_file)

        self.form.upload_btn.setEnabled(True)
        self.form.upload_btn.click()

        self.form.save_btn.click()

        self.assertEqual(self.form.name_lbl.text(), samples_dic["sample_1"]["target_allele"])

        new_ena_file_path = os.path.join(curr_settings["projects_dir"], self.project_name,
                                         samples_dic["sample_1"]["id_int"], samples_dic["sample_1"]["curr_ena_file"])
        reference_file_path = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                           samples_dic["sample_1"]["data_unittest_dir"],
                                           samples_dic["sample_1"]["curr_ena_file"])
        diff_ena_files = compare_2_files(new_ena_file_path, reference_file_path)
        self.assertEqual(len(diff_ena_files["added_sings"]), 0)
        self.assertEqual(len(diff_ena_files["deleted_sings"]), 0)

        self.assertTrue(os.path.exists(
            os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"])))
        self.assertTrue(os.path.exists(
            os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"],
                         samples_dic["sample_1"]["curr_ena_file"])))
        self.assertTrue(os.path.exists(
            os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"],
                         samples_dic["sample_1"]["curr_fasta_file"])))
        self.assertTrue(os.path.exists(
            os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_1"]["id_int"],
                         samples_dic["sample_1"]["curr_blast_file"])))

    def test_xml_file(self):
        """
        Create ENA flatfile from xml
        """
        self.form = ALLELE.NewAlleleForm(log, mydb, self.project_name, curr_settings, None,
                                         samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["id_ext"])
        xml_file = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                samples_dic["sample_2"]["data_unittest_dir"], samples_dic["sample_2"]["input"])
        log.info(f"XML raw file: {xml_file}")
        self.form.file_widget.field.setText(xml_file)
        self.form.upload_btn.setEnabled(True)
        self.form.upload_btn.click()

        self.form.allele1_sec.checkbox.setChecked(True)  # choose second allele

        self.assertEqual(self.form.allele1_sec.gene_field.text(), "HLA-A")
        self.assertEqual(self.form.allele1_sec.GenDX_result, "A*01:01:01:01")
        self.assertEqual(self.form.allele1_sec.name_field.text(), samples_dic["sample_2"]["target_allele"])
        self.assertEqual(self.form.allele1_sec.product_field.text(), "MHC class I antigen")

        self.form.ok_btn.click()
        self.form.save_btn.click()

        new_ena_file_path = os.path.join(curr_settings["projects_dir"], self.project_name,
                                         samples_dic["sample_2"]["id_int"], samples_dic["sample_2"]["curr_ena_file"])
        reference_file_path = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                           samples_dic["sample_2"]["data_unittest_dir"],
                                           samples_dic["sample_2"]["curr_ena_file"])
        print(reference_file_path)
        diff_ena_files = compare_2_files(new_ena_file_path, reference_file_path)
        self.assertEqual(len(diff_ena_files["added_sings"]), 0)
        self.assertEqual(len(diff_ena_files["deleted_sings"]), 0)

        self.assertTrue(os.path.exists(
            os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"])))
        self.assertTrue(os.path.exists(
            os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"],
                         samples_dic["sample_2"]["curr_ena_file"])))
        self.assertTrue(os.path.exists(
            os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"],
                         samples_dic["sample_2"]["curr_fasta_file"])))
        self.assertTrue(os.path.exists(
            os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"],
                         samples_dic["sample_2"]["curr_blast_file"])))
        self.assertTrue(os.path.exists(
            os.path.join(curr_settings["projects_dir"], self.project_name, samples_dic["sample_2"]["id_int"],
                         samples_dic["sample_2"]["curr_gendx_file"])))

    # @unittest.skip("skipping test_fasta_alleles_entries")
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

        query = "select * from alleles"
        success, data_content = execute_db_query(query,
                                                 len(data_info),
                                                 log,
                                                 "Get data from {}",
                                                 "Successful select * from {}",
                                                 "Can't get rows from {}",
                                                 "ALLELES")

        self.assertEqual(len(data_content), 2)  # should be 2 - fasta & xml file
        # sample 1: fasta
        self.assertEqual(data_content[0][0], samples_dic["sample_1"]["id_int"])  # sample_id_in
        self.assertEqual(data_content[0][1], 1)  # allele_nr
        self.assertEqual(data_content[0][2], self.project_name)  # project_name
        self.assertEqual(data_content[0][3], 1)  # project_nr
        self.assertEqual(data_content[0][4], "")  # old cell_line
        self.assertEqual(data_content[0][5], "{}_{}_{}".format(curr_settings["cell_line_token"],
                                                               samples_dic["sample_1"]["id_int"],
                                                               "2DL4_1"))  # local_name
        self.assertEqual(data_content[0][6], "KIR2DL4")  # gene
        self.assertEqual(data_content[0][7], "novel")  # goal
        self.assertEqual(data_content[0][8], "ENA-ready")  # allele_status
        self.assertEqual(data_content[0][14], "completed")  # lab_status
        self.assertEqual(data_content[0][24], samples_dic["sample_1"]["target_allele"])  # target_allele
        self.assertEqual(data_content[0][31], "IPD-KIR")  # reference_database
        # sample 2: xml
        self.assertEqual(data_content[1][0], samples_dic["sample_2"]["id_int"])  # sample_id_in
        self.assertEqual(data_content[1][1], 1)  # allele_nr
        self.assertEqual(data_content[1][2], self.project_name)  # project_name
        self.assertEqual(data_content[1][3], 2)  # project_nr
        self.assertEqual(data_content[1][4], "")  # old cell_line
        self.assertEqual(data_content[1][5], "{}_{}_{}".format(curr_settings["cell_line_token"],
                                                               samples_dic["sample_2"]["id_int"], "A_1"))  # local_name

        self.assertEqual(data_content[1][6], "HLA-A")  # gene
        self.assertEqual(data_content[1][7], "novel")  # goal
        self.assertEqual(data_content[1][8], "ENA-ready")  # allele_status
        self.assertEqual(data_content[1][14], "completed")  # lab_status
        self.assertEqual(data_content[1][20], "yes")  # long_read_data
        self.assertEqual(data_content[1][21], "")  # long_read_phasing
        self.assertEqual(data_content[1][24], samples_dic["sample_2"]["target_allele"])  # target_allele
        self.assertEqual(data_content[1][25], samples_dic["sample_2"]["partner_allele"])  # partner_allele
        self.assertEqual(data_content[1][28], "NGSengine")  # new_genotyping_software
        self.assertEqual(data_content[1][29], "")  # new_software_version
        self.assertEqual(data_content[1][30], "")  # new_genotyping_date
        self.assertEqual(data_content[1][31], "IPD-IMGT/HLA")  # reference_database

    # @unittest.skip("skipping test_fasta_files_entries")
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

        query = "select * from files"
        success, data_content = execute_db_query(query,
                                                 len(data_info),
                                                 log,
                                                 "Get data from {}",
                                                 "Successful select * from {}",
                                                 "Can't get rows from {}",
                                                 "FILES")

        self.assertEqual(len(data_content), 2)  # should be 2 - fasta & xml file
        # sample_1: fasta
        self.assertEqual(data_content[0][0], samples_dic["sample_1"]["id_int"])  # sample_id_in
        self.assertEqual(data_content[0][1], 1)  # allele_nr
        self.assertEqual(data_content[0][3], self.project_name)  # project_name
        self.assertEqual(data_content[0][4], "FASTA")  # raw_file_type
        self.assertEqual(data_content[0][5], samples_dic["sample_1"]["curr_fasta_file"])  # raw_file
        self.assertEqual(data_content[0][6], samples_dic["sample_1"]["curr_fasta_file"])  # fasta
        self.assertEqual(data_content[0][7], samples_dic["sample_1"]["curr_blast_file"])  # blast_xml
        self.assertEqual(data_content[0][8], samples_dic["sample_1"]["curr_ena_file"])  # ena_file
        # sample_2: xml
        self.assertEqual(data_content[1][0], samples_dic["sample_2"]["id_int"])  # sample_id_in
        self.assertEqual(data_content[1][1], 1)  # allele_nr
        self.assertEqual(data_content[1][3], self.project_name)  # project_name
        self.assertEqual(data_content[1][4], "XML")  # raw_file_type
        self.assertEqual(data_content[1][5], samples_dic["sample_2"]["curr_gendx_file"])  # raw_file
        self.assertEqual(data_content[1][6], samples_dic["sample_2"]["curr_fasta_file"])  # fasta
        self.assertEqual(data_content[1][7], samples_dic["sample_2"]["curr_blast_file"])  # blast_xml
        self.assertEqual(data_content[1][8], samples_dic["sample_2"]["curr_ena_file"])  # ena_file

    # @unittest.skip("skipping test_fasta_samples_entries")
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

        query = "select * from samples"
        success, data_content = execute_db_query(query,
                                                 len(data_info),
                                                 log,
                                                 "Get data from {}",
                                                 "Successful select * from {}",
                                                 "Can't get rows from {}",
                                                 "SAMPLES")

        self.assertEqual(len(data_content), 2)  # should be 2 rows - fasta & xml file
        # sample_1: fasta
        self.assertEqual(data_content[0][0], samples_dic["sample_1"]["id_int"])  # sample_id_in
        self.assertEqual(data_content[0][1], samples_dic["sample_1"]["id_ext"])  # sample_id_ext
        self.assertEqual(data_content[0][2], samples_dic["sample_1"]["cell_line"])
        self.assertEqual(data_content[0][3], samples_dic["sample_1"]["customer"])
        self.assertEqual(data_content[0][4], samples_dic["sample_1"]["provenance"])
        self.assertEqual(data_content[0][5], samples_dic["sample_1"]["collection_date"])
        # sample_2: xml
        self.assertEqual(data_content[1][0], samples_dic["sample_2"]["id_int"])  # sample_id_in
        self.assertEqual(data_content[1][1], samples_dic["sample_2"]["id_ext"])  # sample_id_ext
        self.assertEqual(data_content[1][2], samples_dic["sample_2"]["cell_line"])  # sample_id_ext
        for field in [data_content[1][3], data_content[1][4], data_content[1][5]]:  # these are not specified in XML
            self.assertFalse(field)


class Test_Reject_Invalid_Fastas(unittest.TestCase):
    """ tests whether FASTA files with invalid format are rejected correctly
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_Reject_Invalid_Fastas because skip_other_tests is set to True")
        else:
            self.project_name = project_name
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "new_allele_fasta")

    @classmethod
    def tearDownClass(self):
        pass

    def test_fasta_rejected_no_header(self):
        """reject FASTA file without header
        """
        myfile = os.path.join(self.mydir, "invalid1_noHeader.fasta")
        success, msg = typeloader_functions.upload_new_allele_complete(self.project_name, "X", "test", myfile,
                                                                       "DKMS", "Germany", "2024", curr_settings, mydb,
                                                                       log)
        self.assertFalse(success)
        self.assertEqual(msg, "Problem with the FASTA file: FASTA files should have a header starting with >")

    def test_fasta_rejected_no_sequence(self):
        """reject FASTA file without sequence
        """
        myfile = os.path.join(self.mydir, "invalid2_noSeq.fasta")
        success, msg = typeloader_functions.upload_new_allele_complete(self.project_name, "X", "test", myfile,
                                                                       "DKMS", "Germany", "2024", curr_settings, mydb,
                                                                       log)
        self.assertFalse(success)
        self.assertEqual(msg,
                         "Problem with the FASTA file: FASTA files must contain a valid nucleotide sequence after the header!")

    def test_fasta_rejected_empty_header(self):
        """reject FASTA file with empty header
        """
        myfile = os.path.join(self.mydir, "invalid3_emptyHeader.fasta")
        success, msg = typeloader_functions.upload_new_allele_complete(self.project_name, "X", "test", myfile,
                                                                       "DKMS", "Germany", "2024", curr_settings, mydb,
                                                                       log)
        self.assertFalse(success)
        self.assertEqual(msg,
                         "Problem with the FASTA file: This input FASTA file has an empty header! Please put something after the '>'!")


class Test_Send_To_ENA(unittest.TestCase):
    """
    Send both of the fasta and xml samples to ENA
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Submission to ENA because skip_other_tests is set to True")
        else:
            self.project_name = project_name  # "20180710_SA_A_1292"
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                      "ENA_submission")
            self.donor = samples_dic["sample_1"]["id_int"]
            self.donor2 = samples_dic["sample_2"]["id_int"]
            self.provenance = samples_dic["sample_1"]["provenance"]
            self.collection_date = samples_dic["sample_1"]["collection_date"]
            self.customer = samples_dic["sample_1"]["customer"]

            self.form = ENA.ENASubmissionForm(log, mydb, self.project_name, curr_settings, parent=None)

    @classmethod
    def tearDownClass(self):
        pass

    # @unittest.skip("demonstrating skipping")
    def test_submit_to_ENA(self):
        """
        Takes the 2 samples out of Test_Create_New_Allele
        """
        # ok button clickable, if project was chosen [setUpClass]
        self.form.ok_btn.click()

        # they are activated by initialisation, but to be sure...
        self.form.project_files.check_dic[0].setChecked(True)
        self.form.project_files.check_dic[1].setChecked(True)

        self.assertEqual(self.form.project_files.item(0, 1).text(), "1")
        self.assertEqual(self.form.project_files.item(0, 2).text(), samples_dic["sample_1"]["id_int"])
        self.assertEqual(self.form.project_files.item(0, 3).text(), samples_dic["sample_1"]["curr_ena_file"])
        self.assertEqual(self.form.project_files.item(0, 4).text(), "ENA-ready")

        self.assertEqual(self.form.project_files.item(1, 1).text(), "2")
        self.assertEqual(self.form.project_files.item(1, 2).text(), samples_dic["sample_2"]["id_int"])
        self.assertEqual(self.form.project_files.item(1, 3).text(), samples_dic["sample_2"]["curr_ena_file"])
        self.assertEqual(self.form.project_files.item(1, 4).text(), "ENA-ready")

        # submit = send to ENA
        self.form.submit_btn.click()
        # do not write in database, if close btn isn't clicked
        self.form.close_btn.click()

    def test_spatiotemporal_data_string(self):
        """Ensure that the displayed text of update_spatiotemporal_data() is correct."""
        msg_exp = f"""The following 1 samples had already defined values, which were kept:
 - {self.donor}:
    provenance: '{self.provenance}' (database: no data)
    collection date: '{self.collection_date}' (database: no data)
    customer: '{self.customer}' (database: no data)


The following 1 samples had missing values in the oracle database:
 - {self.donor2}: customer, provenance (customer in database: ''), collection date
=> Using 'missing: third party data' (Customer fields for unknown left empty)"""

        self.assertEqual(self.form.spatiotemporal_msg.strip(), msg_exp.strip())

    def test_parse_ena_manifest_and_flatfile(self):
        """Parse the written ena manifest file and flatfile
        """
        path = os.path.join(curr_settings["projects_dir"], self.project_name)
        # neccessary because timestamp is not known:
        manifest_files = list(filter(lambda x: re.search(r'^PRJEB.*manifest.txt', x), os.listdir(path)))
        self.assertTrue(len(manifest_files) > 0)
        submission_file = manifest_files[0]
        file_split = submission_file.split("_")
        submission_file_path = os.path.join(path, submission_file)
        project_id = file_split[0]
        curr_alias = "_".join([file_split[0], file_split[1], "filesub"])
        flatfile = "_".join([file_split[0], file_split[1], "flatfile.txt.gz"])

        # check manifest:
        manifest_ref_file = os.path.join(self.mydir, "manifest.txt")
        log.debug(f"Reference file: {manifest_ref_file}")
        replacements = {"{STUDY}": project_id,
                        "{ALIAS}": curr_alias,
                        "{FLATFILE}": flatfile,
                        "{TL-VERSION}": __version__}
        with open(manifest_ref_file) as f:
            manifest_ref_text = f.read()
            for mykey in replacements:
                manifest_ref_text = manifest_ref_text.replace(mykey, replacements[mykey])
        result = compare_2_files(query_path=submission_file_path, reference_var=manifest_ref_text)
        self.assertEqual(len(result["added_sings"]), 0)
        self.assertEqual(len(result["deleted_sings"]), 0)

        # check flatfile:
        ref_flatfile = os.path.join(self.mydir, "flatfile.txt")
        new_flatfile = os.path.join(path, flatfile)
        with gzip.open(new_flatfile, "rt") as f:
            query_text = f.read()

        result = compare_2_files(query_var=query_text, reference_path=ref_flatfile)
        self.assertEqual(len(result["added_sings"]), 0)
        self.assertEqual(len(result["deleted_sings"]), 0)

    def test_parse_ena_output_and_db_entry(self):
        """
        Parse the written ena output XML file + the db entry
        """
        path = os.path.join(curr_settings["projects_dir"], self.project_name)
        # neccessary, because timestep is not known
        manifest_file = list(filter(lambda x: re.search(r'^PRJEB.*manifest.txt', x), os.listdir(path)))[0]
        file_split = manifest_file.split("_")
        curr_sub_id = "_".join([file_split[0], file_split[1]])
        webin_report = os.path.join(path, f"{curr_sub_id}_webin-cli.report")

        query = "PRAGMA table_info (ENA_SUBMISSIONS)"
        success, data_info = execute_db_query(query,
                                              6,
                                              log,
                                              "Column count at {}",
                                              "Successful read table_info at {}",
                                              "Can't get information from {}",
                                              "ENA_SUBMISSIONS")
        self.assertTrue(success)
        query = "select * from ena_submissions"
        success, data_content = execute_db_query(query,
                                                 len(data_info),
                                                 log,
                                                 "Get data from {}",
                                                 "Successful select * from {}",
                                                 "Can't get rows from {}",
                                                 "ENA_SUBMISSIONS")

        self.assertTrue(success)
        self.assertEqual(len(data_content), 1)
        self.assertEqual(data_content[0][0], self.project_name)  # project_name
        self.assertEqual(data_content[0][1], curr_sub_id)  # submission_id
        self.assertEqual(data_content[0][2], 2)  # nr_alleles
        self.assertEqual(data_content[0][7], "yes")  # success

        acc_analysis = data_content[0][5]
        acc_submission = data_content[0][6]
        self.assertEqual(acc_submission, '')  # is no longer used

        # read webin_cli:
        with open(webin_report, "r") as f:
            text = f.read()
            s = [line for line in text.split("\n") if line]
            if "This was a TEST submission" in s[-1]:
                s.pop(-1)  # remove very last line
            # check penultimate line:
            self.assertTrue("Files have been uploaded to webin" in s[-2] and "ebi.ac.uk." in s[-2])
            # check last line:
            s2 = s[-1].split(
                "The submission has been completed successfully. The following analysis accession was assigned to the submission:")
            self.assertEqual(len(s2), 2)
            submission_id = s2[-1].strip()
            self.assertEqual(submission_id, acc_analysis)


class Test_Send_To_IMGT(unittest.TestCase):
    """Create IMGT / IPD files
    No transmission
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Submission to IPD because skip_other_tests is set to True")
        else:
            log.debug("setting up test...")
            query = "select project_name from projects"
            success, data_content = execute_db_query(query,
                                                     1,
                                                     log,
                                                     "Get data from {}",
                                                     "Successful select * from {}",
                                                     "Can't get rows from {}",
                                                     "PROJECTS")
            self.assertTrue(success, "Could not retrieve project name from table PROJECTS")
            self.project_name = data_content[0][0]
            self.form = IPD.IPDSubmissionForm(log, mydb, self.project_name, curr_settings, parent=None)
            query = f"""update ALLELES set IPD_submission_nr = 'DKMS1000{samples_dic["sample_1"]["submission_id"]}'
                        where Sample_ID_int = '{samples_dic["sample_1"]["id_int"]}' and allele_nr = 1"""
            success, data = execute_db_query(query, 0, log,
                                             "updating {}.IPD_submission_nr",
                                             "Successfully updated {}.IPD_submission_nr",
                                             "Can't update {}.IPD_submission_nr",
                                             "ALLELES")
            if success:
                log.debug("\t=> setup successful")
            else:
                raise Exception("setUpClass did not work! Can't update ALLELES.IPD_submission_nr!")

    @classmethod
    def tearDownClass(self):
        pass

    # @unittest.skip("demonstrating skipping")
    def test_generate_IMGT_files(self):
        """
        Takes the 1st sample out of Test_Create_New_Allele
        """
        # click to proceed to section 2
        self.form.ok_btn1.click()
        befund_file = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                   samples_dic["sample_1"]["input_dir_origin"],
                                   samples_dic["sample_1"]["curr_ipd_befund_file"])
        log.debug("befund-file: {}".format(befund_file))
        self.form.befund_widget.field.setText(befund_file)

        # mock ENA accession retrieval:
        self.form.ENA_id_map = {samples_dic["sample_1"]["local_name"]: samples_dic["sample_1"]["ena_acc"]}
        self.form.ENA_gene_map = {samples_dic["sample_1"]["local_name"]: samples_dic["sample_1"]["gene"]}

        # create IPD files:
        log.debug("Clicking ok_btn2...")
        self.form.ok_btn2.check_ready()
        self.form.ok_btn2.click()

        # they are activated by initialisation, but to be sure...
        self.form.project_files.check_dic[0].setChecked(
            True)  # if this fails, the alleles in the ENA reply file are probably not recognized correctly

        self.assertEqual(self.form.project_files.item(0, 1).text(), "1")
        self.assertEqual(self.form.project_files.item(0, 2).text(), samples_dic["sample_1"]["id_int"])
        self.assertEqual(self.form.project_files.item(0, 3).text(), samples_dic["sample_1"]["local_name"])
        self.assertEqual(self.form.project_files.item(0, 4).text(), "ENA submitted")
        query = "select * from ena_submissions"
        success, data_content = execute_db_query(query,
                                                 2,
                                                 log,
                                                 "Get data from {}",
                                                 "Successful select * from {}",
                                                 "Can't get rows from {}",
                                                 "ENA_SUBMISSIONS")
        self.assertTrue(success, "Could not retrieve data from ENA_submissions")
        self.assertEqual(self.form.project_files.item(0, 5).text(), data_content[0][1])  # submissionID
        self.form.submit_btn.click()
        self.form.ok_btn.click()

    def test_updated_alleles_entries(self):
        """
        Test if ALLELES have been updated correctly
        """
        query = "select * from ena_submissions"
        success, data_content_ena = execute_db_query(query,
                                                     2,
                                                     log,
                                                     "Get data from {}",
                                                     "Successful select * from {}",
                                                     "Can't get rows from {}",
                                                     "ENA_SUBMISSIONS")

        query = "select * from ipd_submissions"
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

        query = "select * from alleles"
        success, data_content = execute_db_query(query,
                                                 len(data_info),
                                                 log,
                                                 "Get data from {}",
                                                 "Successful select * from {}",
                                                 "Can't get rows from {}",
                                                 "ALLELES")

        self.assertEqual(len(data_content), 2)  # should be 2 - fasta & xml file
        # sample 1: fasta
        self.assertEqual(data_content[0][0], samples_dic["sample_1"]["id_int"])  # sample_id_in
        self.assertEqual(data_content[0][1], 1)  # allele_nr
        self.assertEqual(data_content[0][2], self.project_name)  # project_name
        self.assertEqual(data_content[0][3], 1)  # project_nr
        self.assertEqual(data_content[0][4], "")  # old cell_line
        self.assertEqual(data_content[0][5], "{}_{}_{}".format(curr_settings["cell_line_token"],
                                                               samples_dic["sample_1"]["id_int"],
                                                               "2DL4_1"))  # local_name
        self.assertEqual(data_content[0][6], "KIR2DL4")  # gene
        self.assertEqual(data_content[0][7], "novel")  # goal
        self.assertEqual(data_content[0][8], "IPD submitted")  # allele_status
        self.assertEqual(data_content[0][14], "completed")  # lab_status
        self.assertEqual(data_content[0][24], samples_dic["sample_1"]["target_allele"])  # target_allele
        self.assertEqual(data_content[0][31], "IPD-KIR")  # reference_database
        self.assertEqual(data_content[0][32], "2.9.0")  # database_version
        self.assertEqual(data_content[0][36], data_content_ena[0][1])  # ena_submission_id
        self.assertEqual(data_content[0][38], "4ST3TSE6")  # ena accession number: LTxxxxxx
        self.assertEqual(data_content[0][39], data_content_ipd[0][0])  # ipd_submission_id
        self.assertEqual(data_content[0][40],
                         "DKMS1000" + samples_dic["sample_1"]["submission_id"])  # ipd_submission_nr

        ipd_submission_path = os.path.join(curr_settings["projects_dir"], self.project_name,
                                           curr_settings["ipd_submissions"])

        # test ipd submission table in db
        self.assertEqual(len(data_content_ipd), 1)  # should be 1 - one submission with fasta file
        self.assertEqual(data_content_ipd[0][0], os.listdir(ipd_submission_path)[0])  # submission_id
        self.assertEqual(data_content_ipd[0][1], 1)  # number of alleles
        self.assertEqual(data_content_ipd[0][4], "yes")  # success

        ipd_submission_path = os.path.join(curr_settings["projects_dir"], self.project_name, "IPD-submissions",
                                           os.listdir(ipd_submission_path)[0])

        # test, if files are written to the IPD-submissions directory
        befund_file = os.path.join(ipd_submission_path, samples_dic["sample_1"]["curr_ipd_befund_file"])
        ipd_output_file = os.path.join(ipd_submission_path,
                                       "DKMS1000" + samples_dic["sample_1"]["submission_id"] + ".txt")
        output_zip_file = os.path.join(ipd_submission_path, data_content_ipd[0][0] + ".zip")

        self.assertTrue(os.path.exists(befund_file))
        self.assertTrue(os.path.exists(ipd_output_file))
        self.assertTrue(os.path.exists(output_zip_file))

    # @unittest.skip("demonstrating skipping")
    def test_written_ipd_file(self):
        """
        Test the written IPD file (KIR3DP1 in this case) in the submissions directory
        """

        query = "select * from ipd_submissions"
        success, data_content_ipd = execute_db_query(query,
                                                     1,
                                                     log,
                                                     "Get data from {}",
                                                     "Successful select * from {}",
                                                     "Can't get rows from {}",
                                                     "IPD_SUBMISSIONS")

        self.assertTrue(success, "Could not select data from IPD_SUBMISSIONS")

        new_ipd_file_path = os.path.join(curr_settings["projects_dir"], self.project_name,
                                         curr_settings["ipd_submissions"], data_content_ipd[0][0],
                                         "DKMS1000" + samples_dic["sample_1"]["submission_id"] + ".txt")
        reference_file_path = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                           samples_dic["sample_1"]["input_dir_origin"],
                                           "DKMS1000" + samples_dic["sample_1"]["submission_id"] + ".txt")
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
            self.views = {"OVproj": GUI_views_OVprojects.ProjectsOverview(log, mydb, parent=self),
                          "OValleles": GUI_views_OValleles.AllelesOverview(log, mydb),
                          "projectView": GUI_views_project.ProjectView(log, mydb, self.proj_name, parent=self),
                          "sampleView": GUI_views_sample.SampleView(log, mydb, samples_dic["sample_1"]["id_int"],
                                                                    self.proj_name, parent=self)}

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
        self.assertEqual(model.data(model.index(0, 1)), "Open",
                         "Project status in ProjectsOverview unexpected (should be 'Open')")
        self.assertEqual(model.headerData(2, Qt.Horizontal, Qt.DisplayRole), "Creation Date")
        self.assertEqual(model.data(model.index(0, 2)), general.timestamp("%d.%m.%Y"),
                         "Creation Date in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(3, Qt.Horizontal, Qt.DisplayRole), "User Name")
        self.assertEqual(model.data(model.index(0, 3)), project_user, "User name in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(4, Qt.Horizontal, Qt.DisplayRole), "Gene")
        self.assertEqual(model.data(model.index(0, 4)), project_gene, "Gene in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(5, Qt.Horizontal, Qt.DisplayRole), "Pool")
        self.assertEqual(model.data(model.index(0, 5)), project_pool, "Pool in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(6, Qt.Horizontal, Qt.DisplayRole), "Title")
        self.assertEqual(model.data(model.index(0, 6)), project_title, "Project title in ProjectsOverview unexpected")
        self.assertEqual(model.headerData(7, Qt.Horizontal, Qt.DisplayRole), "Description")
        self.assertEqual(model.data(model.index(0, 7)), project_desc,
                         "Project description in ProjectsOverview unexpected")
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
        self.assertEqual(model.data(model.index(1, 21)), '')
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
        self.assertEqual(model.data(model.index(0, 27)), "yes")
        self.assertEqual(model.headerData(28, Qt.Horizontal, Qt.DisplayRole), "Software (new)")
        self.assertEqual(model.data(model.index(0, 28)), "")
        self.assertEqual(model.data(model.index(1, 28)), "NGSengine")
        self.assertEqual(model.headerData(29, Qt.Horizontal, Qt.DisplayRole), "Software Version")
        self.assertEqual(model.data(model.index(0, 29)), "")
        self.assertEqual(model.data(model.index(1, 29)), "")
        self.assertEqual(model.headerData(30, Qt.Horizontal, Qt.DisplayRole), "Genotyping Date")
        self.assertEqual(model.data(model.index(0, 30)), "")
        self.assertEqual(model.data(model.index(1, 30)), "")
        self.assertEqual(model.headerData(31, Qt.Horizontal, Qt.DisplayRole), "Reference Database")
        self.assertEqual(model.data(model.index(0, 31)), "IPD-KIR")
        self.assertEqual(model.data(model.index(1, 31)), "IPD-IMGT/HLA")
        self.assertEqual(model.headerData(32, Qt.Horizontal, Qt.DisplayRole),
                         "Database Version")  # will continually change, therefore not testing content
        self.assertEqual(model.headerData(33, Qt.Horizontal, Qt.DisplayRole), "Internal Allele Name")
        self.assertEqual(model.headerData(34, Qt.Horizontal, Qt.DisplayRole), "Official Allele Name")
        self.assertEqual(model.headerData(35, Qt.Horizontal, Qt.DisplayRole), "New or confirmed?")

        self.assertEqual(model.headerData(36, Qt.Horizontal, Qt.DisplayRole),
                         "ENA Submission ID")  # will continually change, therefore not testing content
        self.assertEqual(model.headerData(37, Qt.Horizontal, Qt.DisplayRole), "ENA Acception Date")
        # self.assertEqual(model.data(model.index(0, 37)), "2020-05-26")
        self.assertEqual(model.data(model.index(1, 37)), "")
        self.assertEqual(model.headerData(38, Qt.Horizontal, Qt.DisplayRole), "ENA Accession Nr")
        self.assertEqual(model.data(model.index(0, 38)), "4ST3TSE6")
        self.assertEqual(model.data(model.index(1, 38)), "")

        self.assertEqual(model.headerData(39, Qt.Horizontal, Qt.DisplayRole),
                         "IPD Submission ID")  # will continually change, therefore not testing content
        self.assertEqual(model.headerData(40, Qt.Horizontal, Qt.DisplayRole), "IPD Submission Nr")
        self.assertEqual(model.data(model.index(0, 40)), "DKMS1000{}".format(samples_dic["sample_1"]["submission_id"]))
        self.assertEqual(model.data(model.index(1, 40)), "")
        self.assertEqual(model.headerData(41, Qt.Horizontal, Qt.DisplayRole), "HWS Submission Nr")
        self.assertEqual(model.headerData(42, Qt.Horizontal, Qt.DisplayRole), "IPD Acception Date")
        self.assertEqual(model.headerData(43, Qt.Horizontal, Qt.DisplayRole), "IPD Release")
        self.assertEqual(model.headerData(44, Qt.Horizontal, Qt.DisplayRole), "Upload Date")
        self.assertEqual(model.data(model.index(0, 44)), general.timestamp("%Y-%m-%d"))
        self.assertEqual(model.headerData(45, Qt.Horizontal, Qt.DisplayRole), "Detection Date")

        self.assertEqual(model.headerData(46, Qt.Horizontal, Qt.DisplayRole), "SAMPLE_ID_INT")
        self.assertEqual(model.headerData(47, Qt.Horizontal, Qt.DisplayRole), "External Sample ID")
        self.assertEqual(model.headerData(48, Qt.Horizontal, Qt.DisplayRole), "Cell Line")
        self.assertEqual(model.headerData(49, Qt.Horizontal, Qt.DisplayRole), "Customer")
        self.assertEqual(model.headerData(50, Qt.Horizontal, Qt.DisplayRole), "Provenance")
        self.assertEqual(model.headerData(51, Qt.Horizontal, Qt.DisplayRole), "Collection Date")
        self.assertEqual(model.headerData(52, Qt.Horizontal, Qt.DisplayRole), "Project Name")
        self.assertEqual(model.headerData(53, Qt.Horizontal, Qt.DisplayRole), "ENA Submission ID")
        self.assertEqual(model.headerData(54, Qt.Horizontal, Qt.DisplayRole), "Alleles in ENA Submission")
        self.assertEqual(model.headerData(55, Qt.Horizontal, Qt.DisplayRole), "Timestamp Sent (ENA Submission)")
        self.assertEqual(model.headerData(56, Qt.Horizontal, Qt.DisplayRole), "Timestamp Confirmed (ENA Submission)")
        self.assertEqual(model.headerData(57, Qt.Horizontal, Qt.DisplayRole), "Analysis Accession Nr")
        self.assertEqual(model.headerData(58, Qt.Horizontal, Qt.DisplayRole), "Submission Accession Nr")
        self.assertEqual(model.headerData(59, Qt.Horizontal, Qt.DisplayRole), "ENA Submission successful?")
        self.assertEqual(model.headerData(60, Qt.Horizontal, Qt.DisplayRole), "IPD Submission ID")
        self.assertEqual(model.headerData(61, Qt.Horizontal, Qt.DisplayRole), "Alleles in IPD Submission")
        self.assertEqual(model.headerData(62, Qt.Horizontal, Qt.DisplayRole), "Timestamp Ready (IPD Submission)")
        self.assertEqual(model.headerData(63, Qt.Horizontal, Qt.DisplayRole), "Timestamp Confirmed (IPD Submission)")
        self.assertEqual(model.headerData(64, Qt.Horizontal, Qt.DisplayRole), "IPD Submission successful?")

        # check expected empty columns:
        empty_columns = [9, 10, 11, 12, 13, 15, 17, 18, 19, 22, 23, 26,
                         33, 34, 35, 41, 42, 43, 63]
        for col in empty_columns:
            for row in [0, 1]:
                self.assertEqual(model.data(model.index(row, col)), "",
                                 "Unexpected value in column {} ({}) row {}: expected empty cell, found '{}'".format(
                                     col, model.headerData(col, Qt.Horizontal, Qt.DisplayRole), row,
                                     model.data(model.index(0, col))))

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
        self.assertEqual(model.headerData(9, Qt.Vertical, Qt.DisplayRole),
                         "ENA Project Submission ID")  # not testing content, but if project_accession is ok, this will be ok, too

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
        self.assertEqual(model.data(model.index(0, 2), Qt.DisplayRole),
                         "{} #{} ({})".format(samples_dic["sample_1"]["id_int"], 1,
                                              samples_dic["sample_1"]["gene"]))
        self.assertEqual(model.data(model.index(1, 2), Qt.DisplayRole),
                         "{} #{} ({})".format(samples_dic["sample_2"]["id_int"], 1,
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
        self.assertEqual(model.data(model.index(0, 2), Qt.DisplayRole),
                         "#{} ({})".format(1, samples_dic["sample_1"]["gene"]))
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
            shown_columns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 14, 23, 33, 34]
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
            self.assertEqual(model.headerData(23, Qt.Vertical, Qt.DisplayRole), "Comment")
            self.assertEqual(model.data(model.index(23, 0), Qt.DisplayRole), "")
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
                                                                             model.headerData(col, Qt.Vertical,
                                                                                              Qt.DisplayRole),
                                                                             value)
                self.assertEqual(value, "", msg)

        def test_tab3_lab(self):
            """tests whether content of SampleView widget 'Details about Allele' tab 'Lab Processing' is correct
            """
            view = widget.tabs[2]
            model = view.flipped_model

            # test whether only right columns are shown:
            shown_columns = [14, 15, 16, 17, 18, 19, 20, 21, 22]
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

            for col in [15, 16, 17, 18, 19, 20, 21, 22]:
                value = model.data(model.index(col, 0), Qt.DisplayRole)
                msg = "Col {} ({}) should be empty, but contains {}!".format(col,
                                                                             model.headerData(col, Qt.Vertical,
                                                                                              Qt.DisplayRole),
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
            self.assertEqual(model.data(model.index(24, 0), Qt.DisplayRole), "KIR2DL4*0080401:new")
            self.assertEqual(model.headerData(25, Qt.Vertical, Qt.DisplayRole), "Partner Allele")
            self.assertEqual(model.headerData(26, Qt.Vertical, Qt.DisplayRole), "Mismatch Position")
            self.assertEqual(model.headerData(27, Qt.Vertical, Qt.DisplayRole), "Null Allele?")
            self.assertEqual(model.data(model.index(27, 0), Qt.DisplayRole), "yes")
            self.assertEqual(model.headerData(28, Qt.Vertical, Qt.DisplayRole), "Software (new)")
            self.assertEqual(model.headerData(29, Qt.Vertical, Qt.DisplayRole), "Software Version")
            self.assertEqual(model.headerData(30, Qt.Vertical, Qt.DisplayRole), "Genotyping Date")
            self.assertEqual(model.headerData(31, Qt.Vertical, Qt.DisplayRole), "Reference Database")
            self.assertEqual(model.data(model.index(31, 0), Qt.DisplayRole), "IPD-KIR")
            self.assertEqual(model.headerData(32, Qt.Vertical, Qt.DisplayRole),
                             "Database Version")  # will change, therefore not testing content
            self.assertEqual(model.headerData(33, Qt.Vertical, Qt.DisplayRole), "Internal Allele Name")
            self.assertEqual(model.headerData(34, Qt.Vertical, Qt.DisplayRole), "Official Allele Name")
            self.assertEqual(model.headerData(35, Qt.Vertical, Qt.DisplayRole), "New or confirmed?")

            for col in [25, 26, 28, 29, 30, 33, 34, 35]:
                value = model.data(model.index(col, 0), Qt.DisplayRole)
                msg = "Col {} ({}) should be empty, but contains {}!".format(col,
                                                                             model.headerData(col, Qt.Vertical,
                                                                                              Qt.DisplayRole),
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
            # TODO: carry timestamp_sent, analysis acc nr and submission acc nr from ENA submission and check for these
            #             if not skip_other_tests:
            #                 self.assertEqual(model.data(model.index(38, 0), Qt.DisplayRole), "{}_{}".format(project_accession, ena_timestamp))
            self.assertEqual(model.headerData(39, Qt.Vertical, Qt.DisplayRole), "Timestamp sent")
            #             if not skip_other_tests:
            #                 self.assertEqual(model.data(model.index(39, 0), Qt.DisplayRole), ena_timestamp)
            self.assertEqual(model.headerData(40, Qt.Vertical, Qt.DisplayRole),
                             "Timestamp confirmed")  # not testing content
            self.assertEqual(model.headerData(41, Qt.Vertical, Qt.DisplayRole), "Analysis accession nr")
            #             self.assertEqual(model.data(model.index(41, 0), Qt.DisplayRole), "ERZ678736")
            self.assertEqual(model.headerData(42, Qt.Vertical, Qt.DisplayRole), "Submission accession nr")
            #             self.assertEqual(model.data(model.index(42, 0), Qt.DisplayRole), "ERA1561070")
            self.assertEqual(model.headerData(43, Qt.Vertical, Qt.DisplayRole), "Submission successful?")
            self.assertEqual(model.data(model.index(43, 0), Qt.DisplayRole), "yes")
            self.assertEqual(model.headerData(44, Qt.Vertical, Qt.DisplayRole), "ENA Acception Date")
            # self.assertEqual(model.data(model.index(44, 0), Qt.DisplayRole), "2020-05-26")
            self.assertEqual(model.headerData(45, Qt.Vertical, Qt.DisplayRole), "ENA Accession Nr")
            self.assertEqual(model.data(model.index(45, 0), Qt.DisplayRole), "4ST3TSE6")

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
            # TODO: carry timestamp_sent, analysis acc nr and submission acc nr from ENA submission and check for these
            #             if not skip_other_tests:
            #                 self.assertEqual(model.data(model.index(39, 0), Qt.DisplayRole), "IPD_{}".format(IPD_timestamp))
            self.assertEqual(model.headerData(40, Qt.Vertical, Qt.DisplayRole), "Timestamp Data Ready")
            #             if not skip_other_tests:
            #                 self.assertEqual(model.data(model.index(40, 0), Qt.DisplayRole), "{}-{}-{}".format(IPD_timestamp[:4], IPD_timestamp[4:6], IPD_timestamp[6:8]))
            self.assertEqual(model.headerData(41, Qt.Vertical, Qt.DisplayRole),
                             "Timestamp Confirmed")  # not testing content
            self.assertEqual(model.headerData(42, Qt.Vertical, Qt.DisplayRole), "Data generated successfully?")
            self.assertEqual(model.data(model.index(42, 0), Qt.DisplayRole), "yes")
            self.assertEqual(model.headerData(43, Qt.Vertical, Qt.DisplayRole), "IPD Submission Nr")
            self.assertEqual(model.data(model.index(43, 0), Qt.DisplayRole),
                             "DKMS1000{}".format(samples_dic["sample_1"]["submission_id"]))
            self.assertEqual(model.headerData(44, Qt.Vertical, Qt.DisplayRole), "HWS Submission Nr")
            self.assertEqual(model.headerData(45, Qt.Vertical, Qt.DisplayRole), "IPD Acception Date")
            self.assertEqual(model.headerData(46, Qt.Vertical, Qt.DisplayRole), "IPD_RELEASE")

            for col in [41, 44, 45, 46]:
                value = model.data(model.index(col, 0), Qt.DisplayRole)
                msg = "Col {} ({}) should be empty, but contains {}!".format(col,
                                                                             model.headerData(col, Qt.Vertical,
                                                                                              Qt.DisplayRole),
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
        test_tab6_IPD(self)  # #TODO: re-enable as soon as final state of IPD tab is decided upon
        test_tab7_history(self)


class TestLogFileDialog(unittest.TestCase):
    """tests whether logfile can be downloaded
    """
    logfile = None
    target_path = None

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping TestLogFileDialog because skip_other_tests is set to True")
        else:
            self.dialog = GUI_download_files.LogFileDialog(curr_settings, log)
            self.logfile = os.path.join(curr_settings["recovery_dir"], "testlog.log")
            self.target_path = "mytestlog.log"
            self.log_text = "This is a test logfile"
            with open(self.logfile, "w") as g:
                g.write(self.log_text)

    @classmethod
    def tearDownClass(self):
        try:
            for myfile in [self.logfile, self.target_path]:
                os.remove(myfile)
        except FileNotFoundError:
            pass

    def test1_path_rejected(self):
        """test whether chosen logfile path is rejected if not in the right directory or wrong filetype
        """
        msg = self.dialog.get_file("test.log")
        self.assertEqual(msg, "This file is not a logfile of the user you started this dialog from!")

        msg = self.dialog.get_file(os.path.join(curr_settings["recovery_dir"], "test.db"))
        self.assertEqual(msg, "This is not a log file! Please choose a file that ends with .log!")

    def test2_path_caught(self):
        """test whether valid logfile path is caught correctly
        """
        self.dialog.get_file(self.logfile)
        self.assertEqual(self.dialog.file, self.logfile)

    def test3_download_logfile(self):
        """test if logfile can be downloaded successfully and correctly
        """
        # remove target file if existant and assure clean slate:
        try:
            os.remove(self.target_path)
        except FileNotFoundError:
            pass
        self.assertFalse(os.path.isfile(self.target_path))

        # download logfile:
        self.dialog.download_file(self.target_path, suppress_messagebox=True)

        # check downloaded logfile:
        self.assertTrue(os.path.isfile(self.target_path))

        with open(self.target_path, "r") as f:
            text = f.read()
            self.assertEqual(text, self.log_text)

        # check that warning is created:
        warning = self.dialog.warning
        self.assertTrue(warning.startswith("Before sending this file to anyone, please open it in a text editor"))
        # TODO: access the QMessagebox, ensure it is actually shown, and then close it


class Test_Make_IMGT_Files_py(unittest.TestCase):
    """
    Test Make_IMGT_Files in typeloader_core
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_Make_IMGT_Files because skip_other_tests is set to True")
        else:
            self.data_dir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                         samples_dic["sample_3"]["input_dir_origin"])
            self.start_num = samples_dic["sample_3"]["submission_id"]
            self.IPD_filename = "DKMS1000" + self.start_num
            self.samples = [
                (samples_dic["sample_3"]["id_int"], samples_dic["sample_3"]["local_name"], self.IPD_filename)]
            self.file_dic = {
                samples_dic["sample_3"]["local_name"]: {"blast_xml": samples_dic["sample_3"]["blast_file_name"],
                                                        "ena_file": samples_dic["sample_3"]["ena_file_name"]}}
            self.ENA_id_map = {samples_dic["sample_3"]["local_name"]: samples_dic["sample_3"]["ena_acc"]}
            self.ENA_gene_map = {samples_dic["sample_3"]["local_name"]: samples_dic["sample_3"]["gene"]}
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


class TestMakeIMGTFilesWith5PrimeOverhang(unittest.TestCase):
    """Test correct handling of fasta file containing a sequence that starts before the reference sequence
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "TestMakeIMGTFilesWith5PrimeOverhang because skip_other_tests is set to True")
        else:
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "start_overhang")
            self.fasta_file = os.path.join(self.mydir, "sequence_starting_28bp_before_reference.fa")
            self.submission_id = "DKMS900001"

            self.project_name = project_name
            self.project_dir = os.path.join(curr_settings["projects_dir"], self.project_name)
            self.sample_id_int = "testIMGT2"
            self.local_name = f'DKMS-LSL_{self.sample_id_int}_2DS3_1'
            self.pretypings = os.path.join(self.mydir, "fake_befunde.csv")

            self.samples = [(self.sample_id_int, self.local_name, '')]
            self.file_dic = {self.local_name: {'blast_xml': f'{self.local_name}.blast.xml',
                                               'ena_file': f'{self.local_name}.ena.txt'}}
            self.allele_dic = {self.local_name: TargetAllele(gene='KIR2DS3',
                                                             target_allele='KIR2DS3*0020103:new',
                                                             partner_allele='KIR2DS3*0010301')}

            self.ENA_id_map = {self.local_name: '15368J48'}
            self.ENA_gene_map = {self.local_name: 'KIR2DS3'}

            self.diff_string = "KIR2DS3*002new differs from KIR2DS3*0020103 like so : Mismatches = pos 447 in codon 128 (TGC -> TGA);pos 6267 (G -> T);pos 6303 (T -> C);pos 6307 (T -> G);pos 6390 (G -> A);pos 6392 (A -> G);. Deletions = pos 15069 (T). Insertions = pos 15038 (C)."

            success, msg = typeloader_functions.upload_new_allele_complete(self.project_name, self.sample_id_int, "bla",
                                                            self.fasta_file, "DKMS", "Germany", "2024",
                                                            curr_settings, mydb, log)
            self.assertTrue(success, "Could not upload sequence_starting_28bp_before_reference.fa")
            log.info(f"{self.fasta_file} uploaded successfully as {self.sample_id_int}")

    @classmethod
    def tearDownClass(self):
        pass

    def test_diff_string_ok(self):
        """test whether correct IPD diff string is generated, including correct codons
        """
        imgt_data, _, _ = MIF.make_imgt_data(self.project_dir, self.samples, self.file_dic, self.allele_dic,
                                             self.ENA_id_map, self.ENA_gene_map, self.pretypings,
                                             curr_settings, log)
        self.assertTrue(imgt_data)
        cell_line = list(imgt_data.keys())[0]
        diff_string = imgt_data[cell_line].split("CC")[1].split("XX")[0].strip()
        self.assertEqual(diff_string, self.diff_string)


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
            self.project_schema = "project"  # do not change
            self.analysis_schema = "analysis"  # do not change
            self.filetype = "flatfile"  # do not change
            self.checksum_method = "MD5"  # do not change

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


class Test_BulkUpload(unittest.TestCase):
    """
    test bulk uploading of fasta sequences
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_BulkUpload because skip_other_tests is set to True")
        else:
            self.project_name = project_name
            self.bulk_file = os.path.join(curr_settings["login_dir"], "data_unittest", "bulk", "bulk_upload.csv")
            self.form = BULK.NewAlleleBulkForm(log, mydb, self.project_name, curr_settings)

    @classmethod
    def tearDownClass(self):
        pass

    def test_bulk_upload(self):
        """
        upload sequences from test file
        """
        self.form.file_widget.field.setText(self.bulk_file)
        self.form.upload_btn.check_ready()
        done = self.form.perform_bulk_upload(auto_confirm=True)
        self.assertTrue(done)

    def test_success(self):
        """make sure expected results occur
        """
        expected_result = """Successfully uploaded 2 of 4 alleles:
  - #2: DKMS-LSL_ID2_bulk_KIR_2DL5B_1
  - #4: DKMS-LSL_ID3_bulk_HLAshort_2DL1_1

Encountered problems in 2 of 4 alleles:
  - #1: Incomplete sequence: This sequence misses the last 400 bp (3' end)!
  - #3: Incomplete sequence: This sequence misses the first 53 bp (5' end) and the last 471 bp (3' end)!

The problem-alleles were NOT added. Please fix them and try again!"""
        result = self.form.report_txt.toPlainText().strip()
        self.assertEqual(result, expected_result)


class Test_provenance_and_collection_date(unittest.TestCase):
    """
    test correct recognition of entries for spatiotemporal data as ok or bad
    """

    @classmethod
    def setUpClass(self):
        oracle_results = {}
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_provenance_and_collection_date because skip_other_tests is set to True")
        else:
            self.ok_dates = ["2024", "2021-01", "2023-12-01", "2222-22-22"]
            self.bad_dates = ["None", "1.1.2002", "01.01.2002", "01. May 2002", "May 2002"]
            self.ok_provenances = ["Germany", "USA", "United Kingdom"]
            self.bad_provenances = ["Blue", "United States", "UK"]
            self.sample_dic = {"ID12230832": {"collection_date": "2015",  # happy
                                              "customer": "DKMS",
                                              "country": "Germany"},
                               "ID18958619": {"collection_date": "2020",  # happy
                                              "customer": "DKMSPL",
                                              "country": "Poland"},
                               "ID18819935": {"collection_date": "2020",  # customer already defined
                                              "customer": "DKMS-BMST",
                                              "country": "India"},
                               "ID17080773": {"collection_date": "2018",  # customer and country already defined
                                              "customer": "X_Test_02",
                                              "country": "missing: third party data"},
                               "ID1": {"customer": "",  # sample not in database
                                       "collection_date": "missing: third party data",
                                       "country": "missing: third party data"}
                               }
            self.already_defined_dic = copy.deepcopy(self.sample_dic)
            self.already_defined_dic["ID17080773"]["customer"] = "Elrond"
            self.already_defined_dic["ID17080773"]["country"] = "Rivendell"
            self.already_defined_dic["ID18819935"]["customer"] = "BMST"

    @classmethod
    def tearDownClass(self):
        pass

    def test_dates(self):
        """
        test that date checks work as expected
        """
        for mydate in self.ok_dates:
            ok, msg = typeloader_functions.check_date(mydate)
            assert ok
            assert not msg

        for mydate in self.bad_dates:
            ok, msg = typeloader_functions.check_date(mydate)
            assert not ok
            assert msg.startswith("Dates need to be formatted ISO8601 compliant")

    def test_provenances_happy(self):
        """test that all approved provenances are ok, and that all 'missing' options are listed last"""
        ok_provenances = typeloader_functions.assemble_country_list(curr_settings, log)

        missing_found = False
        missing_options_listed_last = True

        for country in ok_provenances:
            if country.startswith("missing"):
                missing_found = True
            else:
                if missing_found:
                    missing_options_listed_last = False
        assert missing_options_listed_last

        ok, msg, not_ok_indices = typeloader_functions.check_countries_ok(ok_provenances, curr_settings, log)
        assert ok
        assert not msg
        assert not not_ok_indices

    def test_provenances_sad(self):
        """test correct refusal of incorrect provenance options"""
        ok, msg, not_ok_indices = typeloader_functions.check_countries_ok(self.bad_provenances, curr_settings, log)
        assert not ok
        assert msg == "All items must be exactly spelled like in the official list!\n" \
                      "The following items do not match: \n" \
                      "\t- 'Blue'\t- 'United States'\t- 'UK'"
        assert not_ok_indices == [0, 1, 2]

    def test_oracle_db_results(self):
        samples = list(self.sample_dic.keys())
        results = db_external.get_countries_and_dates_from_oracle_db(samples, log)
        self.assertTrue(results)

        for sample in self.sample_dic:
            self.assertTrue(sample in results)
            (country, year, customer) = results[sample]
            self.assertEqual(country, self.sample_dic[sample]["country"])
            self.assertEqual(year, self.sample_dic[sample]["collection_date"])
            self.assertEqual(customer, self.sample_dic[sample]["customer"])
        self.__class__.oracle_results = results

    def test_integration(self):
        oracle_results = self.__class__.oracle_results
        results = typeloader_functions.integrate_spatiotemporal_data(oracle_results, self.already_defined_dic)
        missing, already_defined, update_queries, result_dic = results
        report = typeloader_functions.report_spatiotemporal_updates(missing, already_defined)

        for id in self.sample_dic:
            self.assertEqual(result_dic[id], self.already_defined_dic[id])

        self.assertEqual(sorted(list(missing.keys())), ['ID1'])
        self.assertEqual(missing['ID1'], ['customer',
                                          "provenance (customer in database: '')",
                                          'collection date'])

        self.assertEqual(sorted(list(already_defined.keys())), ['ID17080773', 'ID18819935'])
        self.assertEqual(already_defined['ID17080773'], [('provenance', 'Rivendell', 'missing: third party data'),
                                                         ('customer', 'Elrond', 'X_Test_02')])
        self.assertEqual(already_defined['ID18819935'], [('customer', 'BMST', 'DKMS-BMST')])

        queries_exp = []
        for sample in sorted(self.already_defined_dic.keys()):
            country = self.already_defined_dic[sample]["country"]
            year = self.already_defined_dic[sample]["collection_date"]
            customer = self.already_defined_dic[sample]["customer"]
            query = f"""update SAMPLES
                        set country = '{country}', collection_date = '{year}', customer = '{customer}'
                        where sample_id_int = '{sample}'
                    """
            queries_exp.append(query)

        self.assertEqual(len(update_queries), len(queries_exp))
        for i in range(len(queries_exp)):
            query = " ".join(update_queries[i].split())  # set all whitespace to ' '
            query_exp = " ".join(queries_exp[i].split())  # set all whitespace to ' '
            self.assertEqual(query, query_exp)

        report_exp = """The following 2 samples had already defined values, which were kept:
 - ID17080773:
    provenance: 'Rivendell' (database: no data)
    customer: 'Elrond' (database: 'X_Test_02')
 - ID18819935: customer: 'BMST' (database: 'DKMS-BMST')


The following 1 samples had missing values in the oracle database:
 - ID1: customer, provenance (customer in database: ''), collection date
=> Using 'missing: third party data' (Customer fields for unknown left empty)
"""
        print(report)
        print(report_exp)
        self.assertEqual(report.strip(), report_exp.strip())


class Test_missing_spatiotemporal_data(unittest.TestCase):
    """
    test correct handling of missing spatiotemporal data
    """

    @classmethod
    def setUpClass(self):
        oracle_results = {}
        if skip_other_tests:
            self.skipTest(self, "Test_missing_spatiotemporal_data because skip_other_tests is set to True")
        else:
            self.sample_dic = {"ID12230832": {"collection_date": "2015",  # happy
                                              "customer": "DKMS",
                                              "country": "Germany"},
                               "ID18958619": {"collection_date": "2020",  # happy
                                              "customer": "DKMSPL",
                                              "country": "Poland"},
                               "ID18819935": {"collection_date": "2020",  # customer already defined
                                              "customer": "DKMS-BMST",
                                              "country": "India"},
                               "ID17080773": {"collection_date": "2018",  # customer and country already defined
                                              "customer": "X_Test_02",
                                              "country": "missing: third party data"},
                               "ID1": {"customer": "",  # sample not in database
                                       "collection_date": "missing: third party data",
                                       "country": "missing: third party data"}
                               }
            self.already_defined_dic = copy.deepcopy(self.sample_dic)
            self.already_defined_dic["ID17080773"]["customer"] = "Elrond"
            self.already_defined_dic["ID17080773"]["country"] = "Rivendell"
            self.already_defined_dic["ID18819935"]["customer"] = "BMST"

    @classmethod
    def tearDownClass(self):
        pass

    def test_dates(self):
        pass


class TestIncompleteSequences(unittest.TestCase):
    """
    test if TypeLoader correctly handles sequences with incomplete UTR3
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_rejection_short_UTR3 because skip_other_tests is set to True")
        else:
            self.mydir = os.path.join(curr_settings["login_dir"], "data_unittest", "incomplete_UTR")
            self.project_name = project_name
            testcase_file = os.path.join(self.mydir, "bulk_upload_incompletes.csv")

            log.debug(f"Reading testcases from file {testcase_file}...")
            TestCase = namedtuple("TestCase", "nr file_name sample_id_int incomplete should_succeed exp_error")
            self.testcases = []
            with open(testcase_file, "r") as f:
                data = csv.reader(f, delimiter=",")
                for row in data:
                    if row:
                        if row[0] != "nr":
                            case = TestCase(nr=row[0], file_name=row[2], sample_id_int=row[3], incomplete=row[6],
                                            should_succeed=row[7], exp_error=row[8])
                            self.testcases.append(case)
            log.debug(f"\t=> found {len(self.testcases)} testcases")

    @classmethod
    def tearDownClass(self):
        pass

    def test_various_cases(self):
        """testing various combinations of incomplete or missing UTRs and the parameter "incomplete_ok"
        """
        for case in self.testcases:
            log.debug(f"\ttesting {case.file_name}...")
            raw_path = os.path.join(self.mydir, case.file_name)
            if case.incomplete == "ok":
                incomplete_ok = True
            else:
                incomplete_ok = False

            if case.should_succeed == "True":
                exp_success = True
            else:
                exp_success = False

            success, msg = typeloader_functions.upload_new_allele_complete(self.project_name, case.sample_id_int,
                                                                           "test", raw_path, "DKMS", "Germany", "2024",
                                                                           curr_settings,
                                                                           mydb, log, incomplete_ok=incomplete_ok)
            self.assertEqual(success, exp_success)
            if not success:
                log.debug(f"Expected error should start with '{case.exp_error}'!")
                log.debug(f"Error encountered: '{msg}'")
                self.assertTrue(msg.startswith(case.exp_error))


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
            self.file_dic = {"C_no_na": ("DKMS-LSL-C-2910.fa", "DKMS-LSL-C-2910.ena.txt"),
                             "C_na": ("DKMS-LSL-C-2952.fa", "DKMS-LSL-C-2952.ena.txt"),
                             "DBP1_no_na": ("DKMS-LSL-DPB1-2887.pacbio.minimap.fa", "DKMS-LSL-DPB1-2887.ena.txt"),
                             "KIR2DL5_no_na": ("DKMS-LSL-KIR2DL5-39.fa", "DKMS-LSL-KIR2DL5-39.ena.txt"),
                             "KIR2DL5_na": ("DKMS-LSL-KIR2DL5-39_corrupt.fa", "DKMS-LSL-KIR2DL5-39_corrupt.ena.txt"),
                             "KIR3DP1_na": ("DKMS-LSL-KIR3DP1-15.fa", "DKMS-LSL-KIR3DP1-15.ena.txt")
                             }

    @classmethod
    def tearDownClass(self):
        pass

    def test_null_allele_algorithm(self):
        """test if null allele algorithm works
        """
        custom_seq_not_null = "ACCGTTTGGTGGTGA"
        custom_seq_stop_codon_1 = "ACCGTTTGGTAATGA"  # TAA
        custom_seq_stop_codon_2 = "ACCGTTTGGTGATGA"  # TGA
        custom_seq_stop_codon_3 = "ACCGTTTGGTAGTGA"  # TAG
        custom_seq_not_divisable = "ACCGTTTGGTGTGA"  # 14 bases

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

            success, myalleles, ENA_text = typeloader_functions.process_sequence_file("PROJECT_NAME", filetype,
                                                                                      blastXmlFile, targetFamily,
                                                                                      fasta_filename, allelesFilename,
                                                                                      header_data, curr_settings, log,
                                                                                      incomplete_ok=True)
            self.assertTrue(success, "Sequence file was not processed successfully")
            result = compare_2_files(reference_path=reference_path, query_var=ENA_text)
            self.assertEqual(len(result["added_sings"]), 0)
            self.assertEqual(len(result["deleted_sings"]), 0)


class Test_MIC(unittest.TestCase):
    """
    Test correct handling of MIC alleles
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_MIC because skip_other_tests is set to True")
        else:
            self.sample_id_int = "MIC001"
            self.local_name = "DKMS-LSL_MIC001_MICA_1"
            self.fasta_name = "MICA_001-7.fa"
            self.target_allele = "MICA*001:new"
            self.partner_allele = ""
            self.gene = "MICA"
            self.ena_acc = "CT5H1R7V"

            self.data_dir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                         "MIC")  # input and reference files
            raw_path = os.path.join(self.data_dir, self.fasta_name)
            self.ENA_id_map = {self.local_name: self.ena_acc}
            self.ENA_gene_map = {self.local_name: self.gene}
            self.pretypings = os.path.join(self.data_dir, "fake_befunde.csv")
            self.ref_ipd_file = os.path.join(self.data_dir, "DKMS10005555.txt")
            log.debug("Pretypings file: {}".format(self.pretypings))

            self.output_dir = os.path.join(curr_settings["projects_dir"], project_name, self.sample_id_int)
            self.xml_file = "{}.blast.xml".format(self.local_name)
            self.ena_file = "{}.ena.txt".format(self.local_name)

            success, msg = typeloader_functions.upload_new_allele_complete(project_name, self.sample_id_int, "test",
                                                                           raw_path,
                                                                           "DKMS", "Germany", "2024",
                                                                           curr_settings, mydb, log,
                                                                           incomplete_ok=True)
            self.assertTrue(success, msg)

    @classmethod
    def tearDownClass(self):
        pass

    def test_file_uploaded_ok(self):
        """tests if MIC file was uploaded successfully and the right data is in the database
        """
        query = """select project_name, local_name, gene, target_allele, null_allele
        from alleles where sample_id_int = '{}'""".format(self.sample_id_int)
        success, data = execute_db_query(query, 5, log, "Get data from ALLELES", "Successful select from {}",
                                         "Can't get rows from {}",
                                         "ALLELES")
        self.assertTrue(success)
        [[myproject, local_name, gene, target_allele, null_allele]] = data
        self.assertEqual(myproject, project_name)
        self.assertEqual(local_name, self.local_name)
        self.assertEqual(gene, "MICA")
        self.assertEqual(target_allele, "MICA*001:new")
        self.assertEqual(null_allele, "yes")

    def test_IPD_file_ok(self):
        """tests whether IPD file can be generated and is correct
        """
        self.start_num = "5555"
        self.IPD_filename = "DKMS1000" + self.start_num
        self.samples = [(self.sample_id_int, self.local_name, self.IPD_filename)]
        self.file_dic = {self.local_name: {"blast_xml": self.xml_file,
                                           "ena_file": self.ena_file}}
        self.allele_dic = {self.local_name:
                               TargetAllele(self.gene,
                                            target_allele=self.target_allele,
                                            partner_allele=self.partner_allele)
                           }

        results = MIF.write_imgt_files(os.path.dirname(self.output_dir),
                                       self.samples, self.file_dic, self.allele_dic, self.ENA_id_map,
                                       self.ENA_gene_map, self.pretypings, self.start_num,
                                       self.output_dir, curr_settings, log)

        self.assertNotEqual(results[0], False)
        (_, cell_lines, _, _, _, success, error) = results

        self.assertEqual(cell_lines[self.local_name], self.IPD_filename)
        self.assertTrue(success)
        self.assertEqual(error, None)

        ipd_submission_file = os.path.join(self.output_dir, self.IPD_filename + ".txt")
        log.debug("  IPD submission file: {}".format(ipd_submission_file))
        self.assertTrue(os.path.exists(ipd_submission_file))

        diff = compare_2_files(ipd_submission_file, self.ref_ipd_file, filetype="IPD")
        self.assertEqual(len(diff["added_sings"]), 0)
        self.assertEqual(len(diff["deleted_sings"]), 0)


class Test_multiple_novel_alleles_part1(unittest.TestCase):
    """
    test if TypeLoader correctly handles multiple novel alleles in one locus
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping Test_multiple_novel_alleles_part1 because skip_other_tests is set to True")
        else:
            self.project_name = project_name
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "multiple")
            self.sample_id_int = "ID000010"
            self.sample_id_ext = "test"
            self.local_name = 'DKMS-LSL_ID000010_3DP1_1'
            self.input_file = os.path.join(self.mydir, "ID000010.fa")
            self.pretypings = os.path.join(self.mydir, "Befunde.csv")
            self.gene = "KIR3DP1"
            self.ena_acc = "LT986596"

            log.info("Check if allele already uploaded...")
            query = "select sample_id_int from alleles where local_name = '{}'".format(self.local_name)
            _, data = execute_db_query(query, 1, log, "Allele already present?",
                                       "Checking ALLELES for test allele {}".format(self.sample_id_int),
                                       "Can't assess ALLELES whether {} is contained".format(self.sample_id_int),
                                       "ALLELES")

            if data:  # set partner allele back to None
                self.log.debug("Setting partner allele back to Null")
                query2 = "update alleles set partner_allele = Null where local_name = '{}'".format(self.local_name)
                _, data = execute_db_query(query2, 0, log, "Return partner allele to empty",
                                           "Resetting ALLELES for test allele {}".format(self.sample_id_int),
                                           "Can't reset ALLELES.partner_allele for {}".format(self.sample_id_int),
                                           "ALLELES")
            else:
                log.info("Uploading allele...")
                self.form1 = ALLELE.NewAlleleForm(log, mydb, self.project_name, curr_settings, None, self.sample_id_int,
                                                  samples_dic["sample_1"]["id_ext"])
                self.form1.file_widget.field.setText(os.path.join(self.mydir, "ID000010.fa"))

                self.form1.upload_btn.setEnabled(True)
                self.form1.upload_btn.click()

                self.form1.save_btn.click()

            # prepare IPDSubmissionForm:
            self.form = IPD.IPDSubmissionForm(log, mydb, self.project_name, curr_settings, parent=None)
            self.form.ENA_id_map = {self.local_name: self.ena_acc}
            self.form.ENA_gene_map = {self.local_name: self.gene}

    @classmethod
    def tearDownClass(self):
        pass

    @patch("typeloader2.GUI_forms_submission_IPD.BothAllelesNovelDialog")
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
        self.form.befund_widget.field.setText(self.pretypings)
        self.form.ok_btn2.check_ready()
        self.form.ok_btn2.click()
        # select alleles:
        log.info("Selecting alleles...")
        if not self.form.project_files.check_dic[0].isChecked():
            self.form.project_files.check_dic[
                0].click()  # if this fails, the alleles in the ENA reply file are probably not recognized correctly
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
        problem_dic = {self.local_name: [self.sample_id_int, self.local_name,
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
                                   "Can't assess ALLELES whether partner_allele for {} has been updated".format(
                                       self.sample_id_int),
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
            self.ena_acc = "LT986596"
            self.gene = "KIR3DP1"
            # prepare IPDSubmissionForm:
            self.form = IPD.IPDSubmissionForm(log, mydb, self.project_name, curr_settings, parent=None)
            self.form.ENA_id_map = {self.local_name: self.ena_acc}
            self.form.ENA_gene_map = {self.local_name: self.gene}

    @classmethod
    def tearDownClass(self):
        pass

    @patch("typeloader2.GUI_forms_submission_IPD.BothAllelesNovelDialog")
    def test_multi_allelesDialogNotCreated(self, mock_dialog):
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
        self.form.befund_widget.field.setText(self.pretypings)
        self.form.ok_btn2.check_ready()
        self.form.ok_btn2.click()

        # select alleles:
        log.info("Selecting alleles...")
        if not self.form.project_files.check_dic[0].isChecked():
            self.form.project_files.check_dic[
                0].click()  # if this fails, the alleles in the ENA reply file are probably not recognized correctly
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
            myfile = os.path.join(self.mydir, "samples.csv")
            log.debug("Reading samples.csv from {}...".format(myfile))
            SampleObject = namedtuple("SampleObject", """name description closest_allele gene
                                                            target_allele partner_allele
                                                            target_family diff_text final_result error_exp""")
            self.samples = {}
            with open(myfile) as f:
                data = csv.reader(f, delimiter=",")
                for i, row in enumerate(data):
                    if i != 0:
                        if row:
                            s = SampleObject(name=row[0],
                                             description=row[1],
                                             closest_allele=row[2],
                                             gene=row[3],
                                             target_allele=row[4],
                                             partner_allele=row[5],
                                             target_family=row[6],
                                             diff_text=row[7],
                                             final_result=row[8],
                                             error_exp=row[9])
                            self.samples[s.name] = s

            # read pretypings from csv:
            log.debug("Reading pretypings.csv...")
            (self.pretypings, self.topics) = MIF.getPatientBefund(os.path.join(self.mydir, "pretypings.csv"))

            # pepare error_dic:
            from typeloader2.typeloader_core import errors
            self.error_dic = {"Cannot tell which novel allele from pretyping this is": errors.BothAllelesNovelError,
                              "no allele marked as new in pretyping": errors.InvalidPretypingError,
                              "assigned allele name not found in pretyping": errors.InvalidPretypingError,
                              "POS is not acceptable pretyping for a target locus": errors.InvalidPretypingError,
                              "pretyping for HLA-B missing": errors.InvalidPretypingError,
                              "Pretyping contains '|'! GL-Strings are only accepted for KIR!": errors.InvalidPretypingError}

    @classmethod
    def tearDownClass(self):
        pass

    def test_pretypings(self):
        """test if pretypings and multiple alleles are handled correctly
        """
        for name in self.samples:
            s = self.samples[name]
            log.debug("testing {}: ({})".format(s.name, s.description))
            target_allele = TargetAllele(gene=s.gene, target_allele=s.target_allele, partner_allele=s.partner_allele)
            self_name = s.target_allele.split("*")[1]
            befund = self.pretypings[name]
            geneMap = {'gene': ['HLA', 'KIR'], 'targetFamily': s.target_family}
            if s.error_exp:  # sad path testing (are correct errors raised?)
                myerror = self.error_dic[s.error_exp]
                # make sure right error type is raised:
                self.assertRaises(myerror, ITG.make_befund_text, befund, self_name, target_allele, s.closest_allele,
                                  geneMap, s.diff_text, log)
                # make sure error text is correct if an error is raised:
                try:
                    ITG.make_befund_text(befund, self_name, target_allele, s.closest_allele,
                                         geneMap, s.diff_text, log)
                except Exception as E:
                    self.assertEqual(E.problem, s.error_exp)

            else:  # happy path testing (these should pass without errors)
                txt = ITG.make_befund_text(befund, self_name, target_allele, s.closest_allele,
                                           geneMap, s.diff_text, log)
                for line in txt.split("\n"):
                    if s.gene in line:
                        self.assertEqual(line, "FT                  /{}".format(s.final_result),
                                         "Error in {}: {}".format(s.name, s.description))  # check result correct

    def test_check_all_required_loci_present_KIR(self):
        KIR_gene = "KIR3DL2"
        KIR_allele_name = "023new"
        KIR_target_allele = TargetAllele(gene='KIR3DL2', target_allele='KIR3DL2*023:new', partner_allele='')
        KIR_befund_ok = '''
            FT                  /KIR2DL1*="001+001"
            FT                  /KIR2DL2*="001+001"
            FT                  /KIR2DL3*="001+001"
            FT                  /KIR2DL4*="001+001"
            FT                  /KIR2DL5A*="001+001"
            FT                  /KIR2DL5B*="001+001"
            FT                  /KIR2DS1*="001+001"
            FT                  /KIR2DS2*="NEG"
            FT                  /KIR2DS3*="POS"
            FT                  /KIR2DS4*="001+001"
            FT                  /KIR2DS5*="001+001"
            FT                  /KIR3DL1*="001+001"
            FT                  /KIR3DL2*="023:new+001"
            FT                  /KIR3DL3*="001+001"
            FT                  /KIR3DS1*="001+001"
            '''
        KIR_befund_bad = '''
            FT                  /KIR2DL1*="001+001"
            FT                  /KIR2DL2*="001+001"
            FT                  /KIR2DL3*="001+001"
            FT                  /KIR2DL4*="001+001"
            FT                  /KIR2DL5*="A001+B001"
            FT                  /KIR2DS1*="001+001"
            FT                  /KIR2DS2*="NEG"
            FT                  /KIR2DS3*="POS"
            FT                  /KIR2DS4*="001+001"
            FT                  /KIR2DS5*="001+001"
            FT                  /KIR3DL1*="001+001"
            FT                  /KIR3DL2*="023:new+001"
            FT                  /KIR3DL3*="001+001"
            FT                  /KIR3DS1*="001+001"
            '''
        KIR_befund_missing = '''
            FT                  /KIR2DL1*="001+001"
            FT                  /KIR2DL4*="001+001"
            FT                  /KIR2DL5A*="001+001"
            FT                  /KIR2DL5B*="001+001"
            FT                  /KIR2DP1*="001+001"
            FT                  /KIR2DS1*="001+001"
            FT                  /KIR2DS2*="001+001"
            FT                  /KIR2DS3*="001+001"
            FT                  /KIR2DS4*="001+001"
            FT                  /KIR2DS5*="001+001"
            FT                  /KIR3DL1*="001+001"
            FT                  /KIR3DL2*="023:new+001"
            FT                  /KIR3DL3*="001+001"
            FT                  /KIR3DP1*="001+001"
            FT                  /KIR3DS1*="001+001"
            '''

        HLA_befund_ok = '''
            FT                  /HLA-A*="01:01:01:01,03:01:02"
            FT                  /HLA-B*="01:01,01:01"
            FT                  /HLA-DRB1*="01:01,01:01"
            FT                  /HLA-DQB1*="01:01,01:01"
            FT                  /HLA-DPB1*="01:01,01:01"
            FT                  /HLA-E*="01:01,01:01"
            FT                  /MICA*="001,001"
            FT                  /MICB*="001,001"
            '''
        HLA_befund_missing = '''
            FT                  /HLA-A*="01:01:01:01,03:01:02"
            FT                  /HLA-DQB1*="01:01,01:01"
            FT                  /HLA-DPB1*="01:01,01:01"
            FT                  /HLA-E*="01:01,01:01"
            FT                  /MICA*="001,001"
            FT                  /MICB*="001,001"
            '''
        HLA_gene = "HLA-A"
        HLA_allele_name = "01:new"
        HLA_target_allele = TargetAllele(gene='HLA-A', target_allele='HLA-A*01:new', partner_allele='HLA-A*03:01:02')

        MIC_befund_ok = '''
                FT                  /HLA-A*="01:01:01:01,03:01:02"
                FT                  /HLA-B*="01:01,01:01"
                FT                  /HLA-DRB1*="01:01,01:01"
                FT                  /HLA-DQB1*="01:01,01:01"
                FT                  /HLA-DPB1*="01:01,01:01"
                FT                  /HLA-E*="01:01,01:01"
                FT                  /MICA*="001:01,003:01"
                FT                  /MICB*="001,001"
                '''
        MIC_befund_missing = '''
                FT                  /HLA-A*="01:01:01:01,03:01:02"
                FT                  /HLA-B*="01:01,01:01"
                FT                  /HLA-DQB1*="01:01,01:01"
                FT                  /HLA-DPB1*="01:01,01:01"
                FT                  /HLA-E*="01:01,01:01"
                FT                  /MICB*="001,001"
                '''
        MIC_gene = "MICA"
        MIC_allele_name = "001:new"
        MIC_target_allele = TargetAllele(gene='MICA', target_allele='MICA*001:new', partner_allele='MICA*003:01')

        param_dic = {"KIR": (KIR_gene, KIR_allele_name, KIR_target_allele),
                     "HLA": (HLA_gene, HLA_allele_name, HLA_target_allele),
                     "MIC": (MIC_gene, MIC_allele_name, MIC_target_allele)}

        for (target, befund_text, exp_err) in [("KIR", KIR_befund_ok, False),
                                               ("KIR", KIR_befund_missing, "pretyping for KIR2DL2 and KIR2DL3 missing"),
                                               ("KIR", KIR_befund_bad, "pretyping for KIR2DL5A and KIR2DL5B missing"),
                                               ("HLA", HLA_befund_ok, False),
                                               ("HLA", HLA_befund_missing, "pretyping for HLA-B and HLA-DRB1 missing"),
                                               ("MIC", MIC_befund_ok, False),
                                               ("MIC", MIC_befund_missing, "pretyping for HLA-DRB1 and MICA missing")
                                               ]:
            (gene, allele_name, target_allele) = param_dic[target]

            if exp_err:
                with self.assertRaises(errors.InvalidPretypingError):
                    ITG.check_all_required_loci(befund_text, gene, target_allele, [], allele_name, log)
                try:
                    ITG.check_all_required_loci(befund_text, gene, target_allele, [], allele_name, log)
                except errors.InvalidPretypingError as E:
                    self.assertEqual(E.problem, exp_err)

            else:
                ok = ITG.check_all_required_loci(befund_text, gene, target_allele, [], allele_name, log)
                self.assertTrue(ok)


class TestEdgecases(unittest.TestCase):
    """
    Test whether differences within the first or last 3 bp of a sequence are caught correctly (#124)
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping TestEdgecases because skip_other_tests is set to True")
        else:
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "outlier_bases")
            self.project_dir = os.path.join(curr_settings["login_dir"], "projects", project_name)
            self.sample_id_int = "EC01"

            # read testcases from csv file:
            self.testcases = []
            TestCase = namedtuple("TestCase", """nr desc filename exp locus target_family closest_allele exact_match
                                                 hit_start align_len query_len del_pos ins_pos mm_pos dels inss mms""")
            with open(os.path.join(self.mydir, "testcases.csv")) as f:
                data = csv.reader(f, delimiter=",")
                for row in data:
                    if row:
                        if row[0] != "nr":
                            nr = row[0]
                            desc = row[1]
                            filename = os.path.join(self.mydir, row[2])
                            exp = row[3]
                            locus = row[4]
                            closest = row[5]
                            exact = True if row[6] == "True" else False  # convert to bool
                            start = int(row[7])
                            align_len = int(row[8])
                            query_len = int(row[9])
                            del_pos = [] if not row[10] else row[10].replace('"', "").split("|")
                            del_pos = [int(x) for x in del_pos]
                            ins_pos = [] if not row[11] else row[11].replace('"', "").split("|")
                            ins_pos = [int(x) for x in ins_pos]
                            mm_pos = [] if not row[12] else row[12].replace('"', "").split("|")
                            mm_pos = [int(x) for x in mm_pos]
                            dels = [] if not row[13] else row[13].replace('"', "").split("|")
                            inss = [] if not row[14] else row[14].replace('"', "").split("|")
                            mms_rough = [] if not row[15] else row[15].replace('"', "").split("|")
                            mms = []
                            if mms_rough:
                                for mystring in mms_rough:  # transform from string to proper tuple
                                    mms.append(tuple(mystring.replace("('", "").replace("')", "").split("', '")))
                            target_fam = row[16]
                            mycase = TestCase(nr=nr, desc=desc, filename=filename, exp=exp, locus=locus,
                                              hit_start=start,
                                              target_family=target_fam, closest_allele=closest, exact_match=exact,
                                              align_len=align_len, query_len=query_len, del_pos=del_pos,
                                              ins_pos=ins_pos,
                                              mm_pos=mm_pos, dels=dels, inss=inss, mms=mms)
                            self.testcases.append(mycase)

            log.info(f"Established {len(self.testcases)} testcases from file.")

    @classmethod
    def tearDownClass(self):
        pass

    def test_edgecases(self):
        """
        testing various cases of changes within 3 bp of either sequence boundary:
        are mm, ins and del positions in these regions identified and located correctly?
        """
        for case in self.testcases:
            log.info(f"Testing case {case.nr}:{case.desc}...")
            self.assertTrue(os.path.exists(case.filename))

            mydic = CA.get_closest_known_alleles(case.filename, case.target_family, curr_settings, log)

            mykey = list(mydic.keys())[0]
            self.assertEqual(mykey, case.closest_allele)
            d = mydic[mykey]
            self.assertEqual(d["name"], case.closest_allele)
            self.assertEqual(d["exactMatch"], case.exact_match)
            self.assertEqual(d["differences"]['deletions'], case.dels)
            self.assertEqual(d["differences"]['insertions'], case.inss)
            self.assertEqual(d["differences"]['mismatches'], case.mms)
            self.assertEqual(d["differences"]['deletionPositions'], case.del_pos)
            self.assertEqual(d["differences"]['insertionPositions'], case.ins_pos)
            self.assertEqual(d["differences"]['mismatchPositions'], case.mm_pos)


class TestDeleteOtherAllele(unittest.TestCase):
    """
    Test whether deletion of non-chosen partner allele of an XML file in subsequent files works
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping TestDeleteOtherAllele because skip_other_tests is set to True")
        else:
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "delete_partner")

            self.blast_xml_file = os.path.join(self.mydir, "DKMS-LSL_ID1_E_1.blast.xml")
            self.fasta_file = os.path.join(self.mydir, "DKMS-LSL_ID1_E_1.fa")
            self.partner_allele = "E*01:03:02:01-Novel-2"

            self.fasta_ref_file = os.path.join(self.mydir, "DKMS-LSL_ID1_E_1_test.fa")
            self.blast_ref_file = os.path.join(self.mydir, "DKMS-LSL_ID1_E_1_test.blast.xml")

    @classmethod
    def tearDownClass(self):
        pass

    def test_deletion_of_unchosen_partner_allele(self):
        """run remove_partner_allele() and check whether partner allele was removed
        """
        from itertools import zip_longest
        typeloader_functions.remove_other_allele(self.blast_xml_file, self.fasta_file, self.partner_allele, log,
                                                 replace=False)
        blast_file_out = self.blast_xml_file + "1"
        fasta_file_out = self.fasta_file + "1"

        log.debug(f"Checking {fasta_file_out}...")
        self.assertTrue(os.path.isfile(fasta_file_out))
        diff = compare_2_files(fasta_file_out, self.fasta_ref_file)
        self.assertEqual(len(diff["added_sings"]), 0)
        self.assertEqual(len(diff["deleted_sings"]), 0)

        log.debug(f"Checking {blast_file_out}...")
        self.assertTrue(os.path.isfile(blast_file_out))
        # compare_2_files() gets stuck in an infinite loop for the blast_xml_files, probably due to XML-ish nature
        with open(blast_file_out) as f1, open(self.blast_ref_file) as f2:
            for (line_new, line_ref) in zip_longest(f1, f2):
                self.assertEqual(line_new, line_ref)

        for new_file in [fasta_file_out, blast_file_out]:
            os.remove(new_file)
            self.assertFalse(os.path.isfile(new_file))


class TestRejectionDeviance(unittest.TestCase):
    """Test whether alleles too different from all known full length alleles are rejected (#138)
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping TestRejectionDeviance because skip_other_tests is set to True")
        else:
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "deviance")
            self.myfile = os.path.join(self.mydir, "2DS1_deviant.fa")

    @classmethod
    def tearDownClass(self):
        pass

    def test_rejection(self):
        """testing whether known deviant allele is rejected correctly
        """
        self.assertTrue(os.path.isfile(self.myfile))
        (_, _, _, _, blastXmlFile, target_family,
         _, _, _) = typeloader_functions.upload_parse_sequence_file(self.myfile, curr_settings, log)
        with self.assertRaises(errors.DevianceError):
            closestAlleles = CA.get_closest_known_alleles(blastXmlFile, target_family, curr_settings, log)


class TestRejectXMLFiles(unittest.TestCase):
    """Test whether unparseable XML files are rejected with proper message (#191)
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping TestRejectXMLFiles because skip_other_tests is set to True")
        else:
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "reject_xml")

    @classmethod
    def tearDownClass(self):
        pass

    def test_rejection_unsuitable(self):
        """testing whether an XML file that is not an NGSEngine file is rejected
        """
        myfile = os.path.join(self.mydir, "unsuitable.xml")
        self.assertTrue(os.path.isfile(myfile))
        with self.assertRaises(errors.UnknownXMLFormatError) as E:
            GASB.getAlleleSequences(myfile, log)
            msg = E.exception.msg
            self.assertTrue(msg.endswith("Hopefully, they can teach TypeLoader to handle it."))

    def test_rejection_multi_loci(self):
        """testing whether an NGSENgine XML file containing > 1 locus is rejected with the correct msg
        """
        myfile = os.path.join(self.mydir, "multiloci.xml")
        self.assertTrue(os.path.isfile(myfile))
        with self.assertRaises(errors.UnknownXMLFormatError) as E:
            GASB.getAlleleSequences(myfile, log)
            msg = E.exception.msg
            self.assertTrue(msg.endswith("input files restricted to one locus and try again with these."))


class TestRestrictedReference(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping TestRestrictedReference because skip_other_tests is set to True")
        else:
            self.project_name = project_name
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"], "restricted_db")
            self.myfile = os.path.join(self.mydir, "problem_allele.xml")
            self.filetype = "XML"
            self.sample_id_int = "restrict_me"
            self.sample_id_ext = "please"
            self.customer = "me"
            self.ref_alleles = ["HLA-B*07:386N", "HLA-B*35:03:01:01"]

            self.restricted_db_good = os.path.join(self.mydir,
                                                   "restricted_db_good")
            self.restricted_db_wrong_target = os.path.join(self.mydir,
                                                           "restricted_db_wrong_target_family")
            self.restricted_db_wrong_alleles = os.path.join(self.mydir,
                                                            "restricted_db_bad")

    @classmethod
    def tearDownClass(self):
        pass

    def test1_reject_at_first(self):
        """test that problem allele is rejected at first because it's too divergent
        """
        success, results1 = typeloader_functions.handle_new_allele_parsing(project_name,
                                                                           self.sample_id_int,
                                                                           self.sample_id_ext,
                                                                           self.myfile,
                                                                           self.customer,
                                                                           curr_settings, log)
        self.assertTrue(success)

        (header_data, filetype, sample_name, targetFamily,
         temp_raw_file, blastXmlFile, fasta_filename, allelesFilename) = results1

        results2 = typeloader_functions.process_sequence_file(project_name, filetype, blastXmlFile,
                                                              targetFamily, fasta_filename,
                                                              allelesFilename, header_data,
                                                              curr_settings, log)

        success = results2[0]
        msg = results2[1]

        self.assertFalse(success)
        self.assertTrue(msg == "Allele too divergent")

    def test2_fail_with_bad_restricted_ref_wrong_target(self):
        """test that allele does not work with a restricted db containing only alleles from wrong
        target family
        """
        myref = self.restricted_db_wrong_target  # contains only KIR
        success, results1 = typeloader_functions.handle_new_allele_parsing(project_name,
                                                                           self.sample_id_int,
                                                                           self.sample_id_ext,
                                                                           self.myfile,
                                                                           self.customer,
                                                                           curr_settings, log,
                                                                           myref)
        self.assertFalse(success)
        self.assertTrue("Did you maybe specify a wrong gene family?" in results1)

    def test3_fail_with_bad_restricted_ref_right_target(self):
        """test that allele does not work with a restricted db containing only non-fitting alleles
        """
        myref = self.restricted_db_wrong_alleles
        success, results1 = typeloader_functions.handle_new_allele_parsing(project_name,
                                                                           self.sample_id_int,
                                                                           self.sample_id_ext,
                                                                           self.myfile,
                                                                           self.customer,
                                                                           curr_settings, log,
                                                                           myref)
        self.assertTrue(success)

        (header_data, filetype, sample_name, targetFamily,
         temp_raw_file, blastXmlFile, fasta_filename, allelesFilename) = results1

        results2 = typeloader_functions.process_sequence_file(project_name, filetype, blastXmlFile,
                                                              targetFamily, fasta_filename,
                                                              allelesFilename, header_data,
                                                              curr_settings, log)

        success = results2[0]
        msg = results2[1]

        self.assertFalse(success)
        self.assertTrue(msg in ['Too many possible alignments', "Allele too divergent"])

    def test4_success_with_good_restricted_db(self):
        """test that allele works when given restricted db with correct alleles
        """
        myref = self.restricted_db_good
        success, results1 = typeloader_functions.handle_new_allele_parsing(project_name,
                                                                           self.sample_id_int,
                                                                           self.sample_id_ext,
                                                                           self.myfile,
                                                                           self.customer,
                                                                           curr_settings, log,
                                                                           myref)
        self.assertTrue(success)

        (header_data, filetype, sample_name, targetFamily,
         temp_raw_file, blastXmlFile, fasta_filename, allelesFilename) = results1

        results2 = typeloader_functions.process_sequence_file(project_name, filetype, blastXmlFile,
                                                              targetFamily, fasta_filename,
                                                              allelesFilename, header_data,
                                                              curr_settings, log)
        self.assertTrue(results2[0])  # success
        [allele1, allele2] = results2[1]
        self.assertEqual(allele1.gendx_result, "B*07:386N-Novel-1")
        self.assertEqual(allele2.gendx_result, "B*35:03:01:01-Existing-2")

    def test5_test_UI(self):
        """test whether restricted db handling via the GUI works correctly
        (happy path test: walk through NewAlleleForm with one file that will not work without
        restricted db but will work with the correct restricted db; check everythins works
        as expected)
        """
        self.form = ALLELE.NewAlleleForm(log, mydb, self.project_name, curr_settings, None,
                                         self.sample_id_int, self.sample_id_ext, testing=True)
        self.dialog = self.form.dialog
        self.assertFalse(self.dialog)  # dialog does not exist, yet

        # enter test file:
        self.form.file_widget.field.setText(self.myfile)
        self.form.upload_btn.setEnabled(True)
        self.form.upload_btn.click()

        # now dialog should be open:
        self.dialog = self.form.dialog
        self.assertTrue(self.dialog)  # dialog is called

        # use dialog to create appropriate restricted reference:
        self.dialog.proceed_btn1.click()
        self.dialog.hla_btn.click()
        self.dialog.proceed_btn2.click()
        for i in range(self.dialog.allele_table.table.rowCount()):
            txt = self.dialog.allele_table.table.item(i, 0).text()
            if txt in self.ref_alleles:
                self.dialog.allele_table.table.selectRow(i)
                self.dialog.allele_table.table.setCurrentCell(i, 0)
                self.dialog.allele_table.table.remember_chosen()

        self.assertEqual(self.dialog.allele_table.count_field.text(), str(len(self.ref_alleles)))
        self.assertEqual(self.dialog.chosen_alleles, self.ref_alleles)

        self.dialog.proceed_btn3.click()
        self.dialog.proceed_btn4.click()

        # continue with NewAlleleForm:
        self.form.allele1_sec.checkbox.setChecked(True)  # choose first allele

        self.form.ok_btn.click()
        self.form.save_btn.click()

        # check results:
        new_ena_file_path = os.path.join(curr_settings["projects_dir"], self.project_name,
                                         self.sample_id_int,
                                         f"DKMS-LSL_{self.sample_id_int}_B_1.ena.txt")
        reference_file_path = os.path.join(self.mydir, "result.ena.txt")

        diff_ena_files = compare_2_files(new_ena_file_path, reference_file_path)
        self.assertEqual(len(diff_ena_files["added_sings"]), 0)
        self.assertEqual(len(diff_ena_files["deleted_sings"]), 0)

    def test6_test_restricted_db_moved(self):
        """test that restricted db files were moved to the sample's dir
        """
        sample_dir = os.path.join(curr_settings["projects_dir"], self.project_name,
                                  self.sample_id_int)
        local_name = f"{curr_settings['cell_line_token']}_{self.sample_id_int}_B_1"

        restricted_db_path_new = os.path.join(sample_dir, f"{local_name}_restricted_db")

        self.assertTrue(os.path.exists(restricted_db_path_new))
        for myfile in ["curr_version_hla.txt", "parsedhla.dump", "parsedhla.fa",
                       "parsedhla.fa.nhr", "parsedhla.fa.nin", "parsedhla.fa.nsq"]:
            self.assertTrue(os.path.exists(os.path.join(restricted_db_path_new, myfile)))

        temp_restricted_db_dir = os.path.join(curr_settings["temp_dir"], "restricted_db")
        self.assertFalse(os.listdir(temp_restricted_db_dir))  # should be empty


class TestHomozygousXML(unittest.TestCase):
    """test whether TypeLoader can handle an XML input file with only 1 allele
    """

    @classmethod
    def setUpClass(self):
        if skip_other_tests:
            self.skipTest(self, "Skipping TestHomozygousXML because skip_other_tests is set to True")
        else:
            self.project_name = project_name
            self.mydir = os.path.join(curr_settings["login_dir"], curr_settings["data_unittest"],
                                      "homozygous_xml")
            self.myfile = os.path.join(self.mydir, "input_file.xml")
            self.filetype = "XML"
            self.sample_id_int = "accept_me"
            self.sample_id_ext = "please"
            self.customer = "me"

    @classmethod
    def tearDownClass(self):
        pass

    def test_xml_file_handled(self):
        """test creating ENA flatfile from homozygous XML file
        """
        self.form = ALLELE.NewAlleleForm(log, mydb, self.project_name, curr_settings, None,
                                         self.sample_id_int, self.sample_id_ext,
                                         testing=True, incomplete_ok=True)
        log.info(f"XML raw file: {self.myfile}")
        self.form.file_widget.field.setText(self.myfile)
        self.form.upload_btn.setEnabled(True)
        self.form.upload_btn.click()

        self.form.allele1_sec.checkbox.setChecked(True)  # choose second allele

        self.assertEqual(self.form.allele1_sec.gene_field.text(), "HLA-B")
        self.assertEqual(self.form.allele1_sec.GenDX_result, "B*51:01:42-Novel-")
        self.assertEqual(self.form.allele1_sec.name_field.text(), "HLA-B*51:new")
        self.assertEqual(self.form.allele1_sec.product_field.text(), "MHC class I antigen")

        self.form.ok_btn.click()
        self.form.save_btn.click()

        new_ena_file_path = os.path.join(curr_settings["projects_dir"], self.project_name,
                                         self.sample_id_int,
                                         f"{curr_settings['cell_line_token']}_{self.sample_id_int}_B_1.ena.txt")
        reference_file_path = os.path.join(self.mydir, "result.ena.txt")

        diff_ena_files = compare_2_files(new_ena_file_path, reference_file_path)
        self.assertEqual(len(diff_ena_files["added_sings"]), 0)
        self.assertEqual(len(diff_ena_files["deleted_sings"]), 0)


class TestCleanStuff(unittest.TestCase):
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


# class TestPlaySound(unittest.TestCase):
#     """Play sound when finished
#     """
#
#     @classmethod
#     def setUpClass(self):
#         pass
#
#     @classmethod
#     def tearDownClass(self):
#         pass
#
#     def test_play_sound(self):
#         """test that sound was played
#          """
#         played = general.play_sound(log)
#         self.assertTrue(played)
#         # make sound distinguishable from normal running sounds, so tester knows tests are through
#         general.play_sound(log)
#         general.play_sound(log)


# ===========================================================
# functions:

def execute_db_query(query, num_columns, log, main_log, success_log, fail_log, format_string=""):
    """
    Executes statements on sqlite database
    returns success, data
    """
    log.info(main_log.format(format_string))
    return db_internal.execute_query(query, num_columns, log,
                                     success_log.format(format_string),
                                     fail_log.format(format_string))


def delete_written_samples(clear_every_row, table, log, column="", value=""):
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


def compare_2_files(query_path="", reference_path="", filetype="", query_var="", reference_var=""):
    result = {}

    if query_path == "":
        query_text = query_var
    else:
        with open(query_path, 'r') as myfile:
            query_text = myfile.read()

    if reference_path == "":
        reference_text = reference_var
    else:
        with open(reference_path, 'r') as myfile:
            reference_text = myfile.read()

    if filetype == "IPD":
        # change the date in order to compare both ipd files
        now = today.strftime("%d/%m/%Y")
        reference_text = re.sub('DT.*Submitted\)\n.*Release\)',
                                'DT   {} (Submitted)\nDT   {} (Release)'.format(now, now), reference_text)

    reference_text = reference_text.replace("{TL-VERSION}", __version__
                                            ).replace("{KIR-VERSION}", curr_settings["db_versions"]["KIR"]
                                                      ).replace("{HLA-VERSION}", curr_settings["db_versions"]["HLA"])

    diffInstance = difflib.Differ()
    diffList = list(diffInstance.compare(query_text.strip(), reference_text.strip()))
    result["added_sings"] = ' '.join(x[2:] for x in diffList if x.startswith('+ '))
    result["deleted_sings"] = ' '.join(x[2:] for x in diffList if x.startswith('- '))

    if result["added_sings"] != "" or result["deleted_sings"] != "":
        log.error("Differences found!")
        log.debug("New file: {}".format(query_path))
        log.debug("Reference file: {}".format(reference_path))
        log.warning("Differences (added): " + result["added_sings"])
        log.warning("Differences (deleted): " + result["deleted_sings"])

        print("Changed lines:")
        query = query_text.split("\n")
        ref = reference_text.split("\n")
        for i in range(len(query)):
            if query[i] != ref[i]:
                print("Query:\t\t", query[i])
                print("Reference:\t", ref[i])
    return result


# functions to order Tests:
# from: https://gist.github.com/catb0t/304ececa6c55f6e3788d
# case_factory: gets all Test_Classes and tests
# suiteFactory: orders the testcases

def suiteFactory(*testcases, testSorter=None, suiteMaker=unittest.makeSuite,
                 newTestSuite=unittest.TestSuite):
    """
    make a test suite from test cases, or generate test suites from test cases.
    *testcases     = TestCase subclasses to work on
    testSorter     = sort tests using this function over sorting by line number
    suiteMaker     = should quack like unittest.makeSuite.
    newTestSuite   = should quack like unittest.TestSuite.
    """

    if testSorter is None:
        log.info("Calling Testcases in order of writing...")
        ln = lambda f: getattr(tc, f).__code__.co_firstlineno
        testSorter = lambda a, b: ln(a) - ln(b)

    test_suite = newTestSuite()
    for tc in testcases:
        test_suite.addTest(suiteMaker(tc, sortUsing=testSorter))

    return test_suite


def caseFactory(
        scope=globals().copy(),
        caseSorter=lambda f: __import__("inspect").findsource(f)[1],
        caseSuperCls=unittest.TestCase,
        caseMatches=__import__("re").compile("^Test")
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
    log.exception(msg="Uncaught Exception", exc_info=(cls, exception, tb))
    # TODO: (future) maybe find a way to display the traceback only once, both in console and logfile?
    sys.__excepthook__(cls, exception, traceback)
    QCoreApplication.exit(1)


pass
# ===========================================================
# main:


sys.excepthook = log_uncaught_exceptions
cases = suiteFactory(*caseFactory())
runner = unittest.TextTestRunner(verbosity=2, failfast=True)
runner.run(cases)
# general.play_sound(log)
# general.play_sound(log)
# general.play_sound(log)
