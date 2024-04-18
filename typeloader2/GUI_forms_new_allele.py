#!/usr/bin/env python3
# -*- coding: cp1252 -*-
'''
Created on 13.03.2018

GUI_forms.py

widgits for adding new sequences or new projects to TypeLoader

@author: Bianca Schoene
'''

# import modules:

import sys, os
from collections import defaultdict

from PyQt5.QtWidgets import (QApplication, QGroupBox, QMessageBox,
                             QGridLayout, QFormLayout, QVBoxLayout,
                             QTextEdit, QLabel, QLineEdit, QCheckBox, QHBoxLayout, QFrame, QComboBox)
from PyQt5.Qt import QWidget, pyqtSlot, pyqtSignal, QDialog, QPushButton
from PyQt5.QtGui import QIcon
from pickle import load

from typeloader2 import general, typeloader_functions as typeloader

try:
    from .typeloader_core import errors
except ImportError:
    from typeloader2.typeloader_core import errors

from typeloader2.GUI_forms import (CollapsibleDialog, ChoiceSection, ChoiceButton, ChoiceTableWidget,
                       FileButton, ProceedButton, QueryButton, NewProjectButton,
                       check_project_open)
from typeloader2.GUI_misc import settings_ok


# ===========================================================
# parameters:

# ===========================================================
# classes:

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


class QueryBox(QDialog):
    """requests data from user that is not given via the file
    """
    sample_data = pyqtSignal(str, str, str, str)

    def __init__(self, log, settings, header_data=None, parent=None):
        super().__init__()
        self.settings = settings
        self.log = log
        if header_data:
            self.header_data = header_data
        else:
            self.header_data = defaultdict(None)
        self.init_UI()
        self.setWindowIcon(QIcon(general.favicon))
        if self.settings["modus"] == "debugging":
            self.fill_with_random_values()

    def init_UI(self):
        self.log.debug("Opening QueryBox...")
        layout = QFormLayout()
        self.setLayout(layout)
        self.title = "Add sample information"

        self.sample_int_field = QLineEdit(self.header_data["SAMPLE_ID_INT"], self)
        layout.addRow(QLabel("Internal Sample-ID:"), self.sample_int_field)

        self.sample_ext_field = QLineEdit(self.header_data["Spendernummer"], self)
        layout.addRow(QLabel("External Sample-ID:"), self.sample_ext_field)

        layout.addRow(QLabel(""))

        separator = QLabel("The following can also be provided later, but are necessary for ENA:")
        separator.setStyleSheet("font-weight: bold")
        layout.addRow(separator)

        self.sample_provenance_field = QComboBox(self)
        countries = typeloader.assemble_country_list(self.settings, self.log)
        self.sample_provenance_field.addItems(countries)
        layout.addRow(QLabel("Sample Provenance:"), self.sample_provenance_field)
        self.sample_provenance_field.setCurrentText(self.header_data["provenance"])

        if self.header_data["collection_date"]:
            self.sample_date_field = QLineEdit(self.header_data["collection_date"], self)
        else:
            self.sample_date_field = QLineEdit(self)
            self.sample_date_field.setPlaceholderText("at least the year; format YYYY-MM-DD")
        layout.addRow(QLabel("Date of Sample Collection:"), self.sample_date_field)

        self.ok_btn = QPushButton("Done", self)
        layout.addRow(self.ok_btn)
        self.ok_btn.clicked.connect(self.on_clicked)

    def on_clicked(self):
        """when ok_btn is clicked, get content of fields and emit it
        """
        self.log.debug("Getting info from query_box")
        sample_int = self.sample_int_field.text().strip()
        sample_ext = self.sample_ext_field.text().strip()
        sample_country = self.sample_provenance_field.currentText().strip()
        sample_date = self.sample_date_field.text().strip()

        if sample_date:
            ok, msg = typeloader.check_date(sample_date)
            if not ok:
                QMessageBox.warning(self, "Bad date", msg)
                return

        if sample_int and sample_ext:
            # TODO: (future) add sanity checks for sample names?
            self.log.debug(f"QueryBox emits ('{sample_int}', '{sample_ext}', '{sample_country}', '{sample_date}')")
            self.sample_data.emit(sample_int, sample_ext, sample_country, sample_date)
            self.close()

    def fill_with_random_values(self):
        """for debugging & development: generate random IDs & put them in QueryBox fields
        """
        import string
        if not self.sample_int_field.text():
            sample_ID = "ID1" + typeloader.id_generator(7, string.digits)
            self.sample_int_field.setText(sample_ID)
        if not self.sample_ext_field.text():
            spender = "DEDKM" + typeloader.id_generator(7, string.digits)
            self.sample_ext_field.setText(spender)


class AlleleSection(QGroupBox):
    """lists details about one allele, derived from input file & editable by user
    """
    selection_changed = pyqtSignal()

    def __init__(self, lbl, parent=None):
        super().__init__(parent)
        if parent:
            self.settings = parent.settings
        else:
            import GUI_login
            self.settings = GUI_login.get_settings("admin", self.log)
        self.lbl = lbl
        self.init_UI()
        self.unselect()
        self.checkbox.toggled.connect(self.toggle_selection)

    def init_UI(self):
        """setup the UI
        """
        self.fields = []
        layout = QGridLayout()
        self.setLayout(layout)
        self.lbl1 = QLabel(self.lbl)
        self.lbl1.setStyleSheet(general.label_style_main)
        self.checkbox = QCheckBox(self)
        layout.addWidget(self.lbl1, 0, 0, 1, 2)
        layout.addWidget(self.checkbox, 0, 2)

        lbl2 = QLabel("Allele details:")
        lbl2.setStyleSheet(general.label_style_2nd)
        layout.addWidget(lbl2, 1, 0, 1, 3)

        self.gene_field = QLineEdit("              ", self)
        self.gene_field.setWhatsThis("Gene of this allele")
        self.fields.append(self.gene_field)
        layout.addWidget(QLabel("\tGene:"), 2, 0)
        layout.addWidget(self.gene_field, 2, 1, 1, 2)

        self.name_field = QLineEdit(self)
        self.fields.append(self.name_field)
        self.name_field.setWhatsThis("Suggested internal name for this allele")
        layout.addWidget(QLabel("\tAllele name:"), 3, 0)
        layout.addWidget(self.name_field, 3, 1, 1, 2)

        self.product_field = QLineEdit(self)
        self.fields.append(self.product_field)
        self.product_field.setWhatsThis("Protein made by this allele, required for ENA")
        layout.addWidget(QLabel("\tProduct:"), 5, 0)
        layout.addWidget(self.product_field, 5, 1, 1, 2)

    def select(self):
        """select the whole section of this allele
        """
        if self.sender != self.checkbox:
            self.checkbox.setChecked(True)
        self.isSelected = True
        self.setStyleSheet(general.groupbox_style_normal)
        for field in self.fields:
            field.setDisabled(False)
        self.selection_changed.emit()

    def unselect(self):
        """unselect the whole section of this allele
        """
        if self.sender != self.checkbox:
            self.checkbox.setChecked(False)
        self.isSelected = False
        self.setStyleSheet(general.groupbox_style_inactive)
        for field in self.fields:
            field.setDisabled(True)
        self.selection_changed.emit()

    def toggle_selection(self):
        """toggle between selected and unselected state
        """
        if self.isSelected:
            self.unselect()
        else:
            self.select()


class NewAlleleForm(CollapsibleDialog):
    """a popup widget to create a new Typeloader Target Allele
    """
    new_allele = pyqtSignal(str)
    refresh_alleles = pyqtSignal(str, str)

    def __init__(self, log, mydb, current_project, settings, parent=None, sample_ID_int=None,
                 sample_ID_ext=None, testing=False, incomplete_ok=False,
                 startover=False):
        self.log = log
        self.mydb = mydb
        self.settings = settings
        if check_project_open(current_project, log, parent=parent):
            self.current_project = current_project
        else:
            self.current_project = ""
        self.startover = startover
        super().__init__(parent)
        log.debug("Opening 'New Allele' Dialog...")
        self.raw_path = None
        self.project = None
        self.parent = parent
        self.testing = testing
        self.incomplete_ok = incomplete_ok
        self.resize(1000, 800)
        if startover:
            self.setWindowTitle("Restart existing allele")
        else:
            self.setWindowTitle("Add new target allele")
        self.setWindowIcon(QIcon(general.favicon))
        self.show()
        self.blastXmlFilename = None
        self.myallele = None
        self.sample_name = sample_ID_int
        self.sample_id_ext = sample_ID_ext
        self.homozygous = False
        self.unsaved_changes = False
        self.upload_btn.check_ready()

        self.dialog = None
        self.restricted_db_path = None
        self.chosen_alleles = None

        ok, msg = settings_ok("new", self.settings, self.log)
        if not ok:
            QMessageBox.warning(self, "Missing settings", msg)
            self.close()

    def define_sections(self):
        """defining the dialog's sections
        """
        self.define_section1()
        self.define_section2()
        self.define_section3()

    def define_section1(self):
        """defining section 1: choose file to upload and project
        """
        mywidget = QWidget(self)
        layout = QHBoxLayout()
        mywidget.setLayout(layout)

        mypath = self.settings["raw_files_path"]
        file_btn = FileButton("Choose XML or Fasta file", mypath, self)
        self.file_widget = ChoiceSection("Raw File:", [file_btn], self.tree)
        self.file_widget.choice.connect(self.get_file)
        mypath = r"C:/Daten/local_data/TL_issue_data/ID14278154.xml"
        if self.settings["modus"] == "debugging":
            self.file_widget.field.setText(mypath)
        layout.addWidget(self.file_widget)

        proj_btn = QueryButton("Choose a (different) existing project",
                               "select project_name from projects where project_status = 'Open' order by project_name desc")
        new_proj_btn = NewProjectButton("Start a new project", self.log, self.mydb, self.settings)
        self.proj_widget = ChoiceSection("Project:", [proj_btn, new_proj_btn], self.tree)
        self.proj_widget.field.setText(self.current_project)
        proj_btn.change_to_normal(None)
        new_proj_btn.change_to_normal(None)
        if self.startover:
            proj_btn.setEnabled(False)
            new_proj_btn.setEnabled(False)

        self.proj_widget.choice.connect(self.get_project)
        layout.addWidget(self.proj_widget)

        self.upload_btn = ProceedButton("Load", [self.file_widget.field, self.proj_widget.field], self.log, 0)
        layout.addWidget(self.upload_btn)
        self.file_widget.choice.connect(self.upload_btn.check_ready)
        self.proj_widget.choice.connect(self.upload_btn.check_ready)
        self.upload_btn.proceed.connect(self.upload_file)

        self.sections.append(("(1) Upload raw file:", mywidget))

    @pyqtSlot(str)
    def get_file(self, file_path):
        """catches name of the file chosen in section1
        """
        self.raw_path = file_path
        self.log.debug("Chose file {}...".format(self.raw_path))

    @pyqtSlot(str)
    def get_project(self, project):
        """catches name of the project chosen in section1
        """
        self.project = project.strip()
        self.log.debug("Chose project {}...".format(self.project))

    @pyqtSlot(str, str, str, str)
    def get_sample_data_from_queryBox(self, sample_ID_int, sample_ID_ext, sample_country, sample_date):
        """accepts the data entered via QueryBox and stores it in self.header_data
        """
        self.header_data["LIMS_DONOR_ID"] = sample_ID_int
        self.sample_name = sample_ID_int
        self.header_data["Spendernummer"] = sample_ID_ext
        self.sample_id_ext = sample_ID_ext
        self.header_data["provenance"] = sample_country
        self.header_data["collection_date"] = sample_date

    @pyqtSlot(int)
    def upload_file(self, _=None):
        """uploads & parses chosen file
        """
        try:
            self.project = self.proj_widget.field.text().strip()
            self.upload_btn.setChecked(True)
            self.upload_btn.setEnabled(False)

            proj_open = check_project_open(self.project, self.log, self)

            if not proj_open:
                msg = f"Project {self.project} is currently closed! You cannot add alleles to closed projects.\n"
                msg += "To add alleles to this project, please open its ProjectView and click the 'Reopen Project' button!"
                msg += "\nAlternatively, please choose a different project."
                self.log.warning(msg)
                QMessageBox.warning(self, "This project is closed!", msg)
                return False

            raw_path = self.file_widget.field.text()

            # upload file to temp dir & parse it:
            self.log.debug("Uploading '{}' to temp dir...".format(os.path.basename(raw_path)))
            success, results = typeloader.handle_new_allele_parsing(self.project, None, None,
                                                                    raw_path, None, self.settings,
                                                                    self.log,
                                                                    self.restricted_db_path)
            if not success:
                msg = results
                self.log.warning("Could not upload new allele")
                self.log.warning(msg)
                if ":" in msg:
                    s = msg.split(":")
                    QMessageBox.warning(self, s[0], ":".join(s[1:]))
                    if s[0] == "Unknown file format!":
                        self.file_widget.reactivate()
                if "Did you maybe specify a wrong gene family?" in msg:
                    self.restricted_db_path = None
                return

            self.log.debug("\t=> success")
            self.success_upload = True
            (self.header_data, self.filetype, sample_name, self.targetFamily,
             self.temp_raw_file, self.blastXmlFile, self.fasta_filename,
             self.allelesFilename) = results

            if not self.sample_name:
                self.sample_name = sample_name
            self.header_data["sample_id_int"] = self.sample_name
            if self.sample_id_ext:
                self.header_data["Spendernummer"] = self.sample_id_ext

            if not self.sample_name:
                self.log.debug("Asking for sample info...")
                self.qbox = QueryBox(self.log, self.settings, self.header_data, self)
                self.qbox.sample_data.connect(self.get_sample_data_from_queryBox)
                self.qbox.exec_()
            if not self.sample_name:
                QMessageBox.warning(self, "No sample name",
                                    "Cannot proceed without sample IDs. Please retry!")
                return

            # process file & create Allele objects:
            self.header_data["sample_id_int"] = self.sample_name
            results = typeloader.process_sequence_file(self.project, self.filetype,
                                                       self.blastXmlFile, self.targetFamily,
                                                       self.fasta_filename, self.allelesFilename,
                                                       self.header_data, self.settings, self.log,
                                                       startover=self.startover)
            if not results[0]:  # something went wrong
                if results[1] == "Incomplete sequence":
                    reply = QMessageBox.question(self, results[1], results[2], QMessageBox.Yes |
                                                 QMessageBox.No, QMessageBox.No)
                    if reply == QMessageBox.Yes:
                        results = typeloader.process_sequence_file(self.project, self.filetype,
                                                                   self.blastXmlFile,
                                                                   self.targetFamily,
                                                                   self.fasta_filename,
                                                                   self.allelesFilename,
                                                                   self.header_data,
                                                                   self.settings,
                                                                   self.log, incomplete_ok=True,
                                                                   startover=self.startover)
                        if not results[0]:
                            QMessageBox.warning(self, results[1], results[2])
                            return
                    else:
                        return
                elif results[1] in ["Allele too divergent", "Too many possible alignments"]:
                    if self.restricted_db_path:
                        msg = "The file you're trying to upload could still not be handled "
                        msg += "with the given allele(s) as reference.\n\n"
                        msg += "Do you want to try again?"
                        reply = QMessageBox.question(self, "Allele still too divergent", msg,
                                                     QMessageBox.Yes | QMessageBox.No,
                                                     QMessageBox.No)
                        if reply == QMessageBox.No:
                            QMessageBox.information(self, "Allele too divergent",
                                                    results[2])
                            self.dialog.close()
                            return
                        else:
                            self.restricted_db_path = None
                            self.chosen_alleles = None
                            if self.dialog:
                                self.dialog.close()

                    self.dialog = ChooseReferenceAllelesDialog(self.log, self.settings,
                                                               self.filetype, self)
                    self.dialog.restricted_db_path.connect(self.catch_restricted_db_path)
                    self.dialog.restricted_alleles.connect(self.catch_chosen_alleles_for_restricted_db)
                    if not self.testing:
                        self.dialog.exec_()
                    return
                else:
                    QMessageBox.warning(self, results[1], results[2])
                    return

            self.success_parsing, self.myalleles, self.ENA_text = results
            if self.success_parsing:
                if self.startover:
                    if self.myalleles[0].gene != self.startover["gene"]:
                        msg = f"{self.startover['local_name']} is a(n) {self.startover['gene']} allele! The uploaded file contains a(n) " \
                              f"{self.myalleles[0].gene} sequence!\nRestarting an allele is only allowed with the " \
                              f"same locus, otherwise the allele name would not match the sequence."
                        QMessageBox.warning(self, "Locus does not match!", msg)
                        return
                if self.filetype == "XML":
                    self.allele1 = self.myalleles[0]
                    self.allele2 = self.myalleles[1]
                    self.upload_btn.setChecked(False)
                    self.upload_btn.setEnabled(True)
                    self.proceed_sections(0, 1)
                    self.fill_section2()

                else:  # Fasta File: move instantly to section 3:
                    self.myallele = self.myalleles[0]
                    self.ENA_widget.setText(self.ENA_text)
                    self.name_lbl.setText(self.myallele.newAlleleName)
                    self.upload_btn.setEnabled(True)
                    self.upload_btn.setChecked(False)
                    for warning in self.myalleles[0].warnings:
                        QMessageBox.warning(self, "Potential problem with annotating this allele", warning)
                    self.proceed_sections(0, 2)
            else:
                QMessageBox.warning(self, results[1], results[2])
        except Exception as E:
            self.log.error(E)
            self.log.exception(E)
            self.close()

    @pyqtSlot(str)
    def catch_restricted_db_path(self, restricted_db_path):
        """if a restricted db is created, store its path, then re-try allele upload

        :param restricted_db_path: path to the directory containing the restricted db
        :return: nothing
        """
        self.restricted_db_path = restricted_db_path
        if not restricted_db_path:
            self.log.debug("No restricted reference database was created")
            self.dialog.close()
            self.dialog = None
        else:
            self.log.debug(f"A restricted reference database was created under {restricted_db_path}")
            self.log.debug("Re-attempting file upload...")
            self.upload_file()

    @pyqtSlot(list)
    def catch_chosen_alleles_for_restricted_db(self, chosen_alleles):
        """if a restricted db is created, store the alleles chosen,
        so they can later be stored as comment in the database

        :param chosen_alleles: list of allele names chosen for the restricted db
        """
        if chosen_alleles:
            self.log.debug(f'Storing chosen reference alleles: {" & ".join(chosen_alleles)}')
            self.chosen_alleles = chosen_alleles
        else:
            self.log.debug("No reference alleles chosen")

    def define_section2(self):
        """defining section 2: Specify allele details
        """
        mywidget = QWidget(self)
        layout = QGridLayout()
        mywidget.setLayout(layout)

        a1_new = False
        self.allele1_sec = AlleleSection("Allele 1:", self)
        layout.addWidget(self.allele1_sec, 0, 0)

        a2_new = False
        self.allele2_sec = AlleleSection("Allele 2:", self)
        layout.addWidget(self.allele2_sec, 0, 1)

        # ToDo: add closest alleles!

        button_widget = QFrame(self)  # contains both-checkbox & proceed-button
        layout2 = QFormLayout()
        button_widget.setLayout(layout2)

        layout2.addRow(QLabel("\n\n"))

        if a1_new:
            self.allele1_sec.checkbox.setChecked(True)
            if a2_new:
                self.allele2_sec.checkbox.setChecked(True)
        elif a2_new:
            self.allele2_sec.checkbox.setChecked(True)
        self.allele1_sec.checkbox.clicked.connect(self.unselect_other_box)
        self.allele2_sec.checkbox.clicked.connect(self.unselect_other_box)

        self.ok_btn = ProceedButton("Proceed", [self.allele1_sec.checkbox, self.allele2_sec.checkbox], self.log,
                                    only1=True)
        self.ok_btn.check_ready()
        self.ok_btn.clicked.connect(self.make_ENA_file)
        self.allele1_sec.selection_changed.connect(self.ok_btn.check_ready)
        self.allele2_sec.selection_changed.connect(self.ok_btn.check_ready)

        layout2.addRow(self.ok_btn)
        layout.addWidget(button_widget, 0, 3)
        layout.setColumnStretch(0, 1)
        layout.setColumnStretch(1, 1)
        layout.setColumnStretch(2, 0)
        self.sections.append(("(2) Specify allele details:", mywidget))

    @pyqtSlot()
    def fill_section2(self):
        """fill fields in section 2 from results of parsing the file of section 1
        """
        self.allele1_sec.lbl1.setText("GenDX: " + self.allele1.gendx_result)
        self.allele1_sec.GenDX_result = self.allele1.gendx_result
        self.allele1_sec.gene_field.setText(self.allele1.gene)
        self.allele1_sec.name_field.setText(self.allele1.name)
        self.allele1_sec.product_field.setText(self.allele1.product)
        if "Novel" in self.allele1.gendx_result:
            self.allele1_sec.checkbox.setChecked(True)

        if self.allele2.gendx_result:
            gendx_result2 = "GenDX: " + self.allele2.gendx_result
            self.homozygous = False
        else:  # homozygous sample has no second allele
            gendx_result2 = ""
            self.homozygous = True
        self.allele2_sec.lbl1.setText(gendx_result2)
        self.allele2_sec.GenDX_result = self.allele2.gendx_result
        self.allele2_sec.gene_field.setText(self.allele2.gene)
        self.allele2_sec.name_field.setText(self.allele2.name)
        self.allele2_sec.product_field.setText(self.allele2.product)
        if "Novel" in self.allele2.gendx_result:
            self.allele2_sec.checkbox.setChecked(True)
        if self.homozygous:
            self.allele2_sec.checkbox.setChecked(False)
            self.allele2_sec.checkbox.setEnabled(False)

        if self.homozygous:
            msg = "This XML file contains only one allele. This is probably because of an allelic "
            msg += "dropout in your sample, which might lead to erroneous sequences!\n\n"
            msg += "Are you really sure you want to upload this file anyway?"
            self.log.warning(msg.replace("\n\n", " "))
            if not self.testing:
                reply = QMessageBox.question(self, "Homozygous sample", msg,
                                             QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
                if reply == QMessageBox.No:
                    self.log.info("User chose to abort.")
                    self.proceed_sections(1, 0)
                else:
                    self.log.info("User chose to continue.")

    @pyqtSlot()
    def select_both(self):
        """select or deselect both alleles at once using the "both" checkbox
        """
        if self.both_cbx.checkState():
            self.log.debug("Selecting both alleles..")
            self.allele1_sec.select()
            self.allele2_sec.select()
            self.msg.setStyleSheet(general.label_style_attention)
            self.ok_btn.check_ready()
        else:
            self.log.debug("Unselecting both alleles..")
            self.allele1_sec.unselect()
            self.allele2_sec.unselect()
            self.msg.setStyleSheet(general.label_style_normal)
            self.ok_btn.check_ready()

    @pyqtSlot()
    def unselect_other_box(self):
        """enforce that only one allele can be accepted at once
        """
        if self.homozygous:
            return
        if self.sender() == self.allele1_sec.checkbox:
            if self.allele1_sec.isSelected:
                self.allele2_sec.checkbox.setChecked(False)
            else:
                self.allele2_sec.checkbox.setChecked(True)
        elif self.sender() == self.allele2_sec.checkbox:
            if self.allele2_sec.isSelected:
                self.allele1_sec.checkbox.setChecked(False)
            else:
                self.allele1_sec.checkbox.setChecked(True)

    @pyqtSlot()
    def unselect_both_cbx(self):
        """if either allele is unselected manually but the both-checkbox is selected,
        unselect the both-checkbox
        """
        if self.both_cbx.checkState():
            if not (self.allele1_sec.isSelected and self.allele2_sec.isSelected):
                self.both_cbx.setChecked(False)
                self.msg.setStyleSheet(general.label_style_normal)
                self.ok_btn.check_ready()

    @pyqtSlot()
    def make_ENA_file(self):
        """creates the file for ENA out of an XML file
        """
        self.log.info("Creating EMBL file...")
        try:
            # get GUI data:
            self.allele1.geneName = self.allele1_sec.gene_field.text().strip()
            self.allele1.alleleName = self.allele1_sec.GenDX_result
            self.allele1.newAlleleName = self.allele1_sec.name_field.text().strip()
            self.allele1.productName_DE = self.allele1_sec.product_field.text().strip()
            self.allele1.productName_FT = self.allele1.productName_DE
            self.allele1.partner_allele = self.allele2_sec.name_field.text().strip()

            self.allele2.geneName = self.allele2_sec.gene_field.text().strip()
            self.allele2.alleleName = self.allele2_sec.GenDX_result
            self.allele2.newAlleleName = self.allele2_sec.name_field.text().strip()
            self.allele2.productName_DE = self.allele2_sec.product_field.text().strip()
            self.allele2.productName_FT = self.allele2.productName_DE
            self.allele2.partner_allele = self.allele1_sec.name_field.text().strip()

            if self.allele1_sec.checkbox.checkState():
                self.myallele = self.allele1
                other_allele_name = self.allele2.alleleName
                self.log.debug("Choosing allele 1...")
            elif self.allele2_sec.checkbox.checkState():
                self.myallele = self.allele2
                other_allele_name = self.allele1.alleleName
                self.log.debug("Choosing allele 2...")
                # TODO: (future) implement possibility to add both alleles
            else:
                QMessageBox.warning(self, "No allele chosen", "Please choose an allele to continue!")
                return

            if self.startover:
                self.myallele.local_name = self.startover["local_name"]
            typeloader.remove_other_allele(self.blastXmlFile, self.fasta_filename, other_allele_name, self.log)
            try:
                self.ENA_text = typeloader.make_ENA_file(self.blastXmlFile, self.targetFamily, self.myallele,
                                                         self.settings, self.log)
            except errors.IncompleteSequenceWarning as E:
                if self.incomplete_ok and self.testing:
                    proceed = True
                else:
                    proceed = False
                    reply = QMessageBox.question(self, "Incomplete Sequence", E.msg, QMessageBox.Yes | QMessageBox.No,
                                                 QMessageBox.No)
                    if reply == QMessageBox.Yes:
                        proceed = True
                if proceed:
                    try:
                        self.ENA_text = typeloader.make_ENA_file(self.blastXmlFile, self.targetFamily, self.myallele,
                                                                 self.settings, self.log, incomplete_ok=True)
                    except errors.MissingUTRError as E:
                        QMessageBox.warning(self, "Missing UTR", E.msg)
                        return
                    except Exception as E:
                        self.log.error(E)
                        self.log.exception(E)
                        QMessageBox.warning(self, "Error during ENA file creation", repr(E))
                        return
                else:
                    return

            except errors.MissingUTRError as E:
                QMessageBox.warning(self, "Missing UTR", E.msg)
                return
            except Exception as E:
                self.log.error(E)
                self.log.exception(E)
                QMessageBox.warning(self, "Error during ENA file creation", repr(E))
                return

            self.ENA_widget.setText(self.ENA_text)
            self.name_lbl.setText(self.myallele.newAlleleName)
            for warning in self.myalleles[0].warnings:
                QMessageBox.warning(self, "Potential problem with annotating this allele", warning)
            self.proceed_sections(1, 2)

        except Exception as E:
            self.log.error(E)
            self.log.exception(E)
            self.close()

    def define_section3(self):
        """defining section 3: check ENA-file & save allele
        """
        mywidget = QWidget(self)
        layout = QGridLayout()
        mywidget.setLayout(layout)

        # TODO: (future) implement option to display both alleles
        self.name_lbl = QLabel()
        self.name_lbl.setStyleSheet(general.label_style_2nd)
        layout.addWidget(self.name_lbl, 0, 0)

        self.ENA_widget = QTextEdit(self)
        self.ENA_widget.textChanged.connect(self.on_text_changed)
        layout.addWidget(self.ENA_widget, 1, 0, 1, 6)
        self.ENA_widget.setMinimumHeight(500)

        save_txt = "Save new target allele"
        if self.startover:
            save_txt = "Replace target allele"
        self.save_btn = ProceedButton(save_txt, [self.ENA_widget], self.log, 2, self)
        layout.addWidget(self.save_btn, 0, 5)
        self.save_btn.proceed.connect(self.save_allele)

        self.save_changes_btn = ProceedButton("Save changes!", [self.ENA_widget], self.log, parent=self)
        layout.addWidget(self.save_changes_btn, 3, 0, 1, 3)
        self.save_changes_btn.clicked.connect(self.save_changes)

        self.discard_btn = ProceedButton("Discard changes!", [self.ENA_widget], self.log, parent=self)
        layout.addWidget(self.discard_btn, 3, 3, 1, 3)
        self.discard_btn.clicked.connect(self.discard_changes)

        self.sections.append(("(3) Check ENA file and save allele:", mywidget))

    @pyqtSlot()
    def on_text_changed(self):
        """handle text edits in ENA text window
        """
        try:
            self.save_btn.check_ready()
            if self.ENA_widget.toPlainText() != self.ENA_text:
                self.save_changes_btn.setEnabled(True)
                self.save_changes_btn.setStyleSheet(general.btn_style_clickme)
                self.discard_btn.setEnabled(True)
                self.discard_btn.setStyleSheet(general.btn_style_clickme)
                self.unsaved_changes = True
        except Exception as E:
            self.log.error(E)
            self.log.exception(E)

    @pyqtSlot()
    def save_changes(self):
        """saves the edited file
        """
        self.log.debug("'Save changes' was clicked")
        try:
            txt = self.ENA_widget.toPlainText()
            if txt:
                self.log.debug("Saving changes?")
                reply = QMessageBox.question(self, 'Message',
                                             "Save changes to ENA-file for {}?".format(self.allele_name),
                                             QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)

                if reply == QMessageBox.Yes:
                    self.log.debug("Saving changes...")
                    self.ENA_text = txt
                    self.unsaved_changes = False
                    self.discard_btn.change_to_normal()
                    self.save_changes_btn.change_to_normal()
                    self.log.debug("=> Success")
                else:
                    self.log.debug("Not saving")
                    self.save_changes_btn.setChecked(False)
                    self.save_changes_btn.setStyleSheet(general.btn_style_clickme)
            else:
                self.log.debug("No text to save...")
        except Exception as E:
            self.log.error(E)
            self.log.exception(E)

    @pyqtSlot()
    def discard_changes(self):
        """reverts changes made to the displayed ENA text
        """
        try:
            self.log.debug("Discarding changes...")
            self.ENA_widget.setText(self.ENA_text)
            self.unsaved_changes = False
            self.discard_btn.change_to_normal()
            self.save_changes_btn.change_to_normal()
            self.log.debug("=> Success")
        except Exception as E:
            self.log.error(E)
            self.log.exception(E)

    @pyqtSlot(int)
    def save_allele(self, _):
        """saves the allele & closes the dialog
        """
        try:
            self.log.debug("Asking for confirmation before saving...")
            self.project = self.proj_widget.field.text()
            if self.project:
                if self.settings["modus"] == "staging":
                    reply = QMessageBox.Yes  # for automatic runthrough
                else:
                    reply = QMessageBox.question(self, 'Message',
                                                 "Save allele {} to project {}?".format(self.myallele.local_name,
                                                                                        self.project), QMessageBox.Yes |
                                                 QMessageBox.No, QMessageBox.No)
                if reply == QMessageBox.Yes:

                    if self.startover:
                        self.log.info("Renaming old files...")
                        for (old_file, new_file) in self.startover["rename_files"]:
                            try:
                                os.rename(old_file, new_file)
                            except FileNotFoundError:
                                self.log.debug(f"File {old_file} not found, probably already renamed")

                    # save sample files:

                    results = typeloader.save_new_allele(self.project, self.sample_name,
                                                         self.myallele.local_name, self.ENA_text,
                                                         self.filetype, self.temp_raw_file,
                                                         self.blastXmlFile, self.fasta_filename,
                                                         self.restricted_db_path,
                                                         self.settings, self.log)
                    (success, err_type, msg, files) = results
                    if not success:
                        QMessageBox.warning(self, err_type, msg)
                        self.save_btn.setStyleSheet(general.btn_style_normal)
                        self.save_btn.setChecked(False)
                        return False

                    # save to db & emit signals:
                    [self.raw_file, self.fasta_filename, self.blastXmlFile, self.ena_path] = files
                    (success, err_type, msg) = typeloader.save_new_allele_to_db(self.myallele,
                                                                                self.project,
                                                                                self.filetype,
                                                                                self.raw_file,
                                                                                self.fasta_filename,
                                                                                self.blastXmlFile,
                                                                                self.header_data,
                                                                                self.targetFamily,
                                                                                self.ena_path,
                                                                                self.chosen_alleles,
                                                                                self.settings,
                                                                                self.mydb, self.log,
                                                                                self.startover)
                    if success:
                        self.new_allele.emit(self.sample_name)
                        self.refresh_alleles.emit(self.project, self.sample_name)
                        self.close()
                    else:
                        self.log.info("Not successful!")
                        self.log.warning(err_type)
                        self.log.info(msg)
                        if err_type:  # if QMessageBox has already been shown, err_type is false => do nothing
                            self.save_btn.setStyleSheet(general.btn_style_normal)
                            self.save_btn.setChecked(False)
                            QMessageBox.warning(self, err_type, msg)
                            self.log.warning(msg)
                else:
                    self.log.debug("Not saving.")
                    self.save_btn.setChecked(False)
            else:
                QMessageBox.warning(self, "Project missing!",
                                    "Please specify project first!")
                self.proceed_sections(2, 0)
        except Exception as E:
            self.log.error(E)
            self.log.exception(E)


class ChooseReferenceAllelesDialog(CollapsibleDialog):
    """a dialog offering to create a resticted reference db containing only the chosen alleles
    """
    restricted_db_path = pyqtSignal(str)
    restricted_alleles = pyqtSignal(list)

    def __init__(self, log, settings, filetype, parent=None):
        self.settings = settings
        self.log = log
        self.filetype = filetype
        super().__init__(parent)

        self.setMinimumWidth(600)
        self.setMinimumHeight(600)
        self.setWindowTitle("Manual choice of reference alleles necessary")
        self.setWindowIcon(QIcon(general.favicon))

        self.chosen_alleles = []

    def define_sections(self):
        self.log.debug("Opening ChooseReferenceAllelesDialog...")
        self.define_section1()
        self.define_section2()
        self.define_section3()
        self.define_section4()

    def define_section1(self):
        """defines section 1, where the user is informed about the situation and asked whether to
        proceed with a restricted database, or abort
        """
        mywidget = QWidget(self)
        layout = QVBoxLayout()
        mywidget.setLayout(layout)

        intro_txt1 = "TypeLoader cannot find a suitable reference gene in the reference database.\n"
        intro_txt1 += "This may be because this allele is too dissimilar to all known "
        intro_txt1 += "reference alleles."
        intro_txt2 = "Do you know the correct reference alleles for the alleles contained in this "
        intro_txt2 += "file? \nIf yes, you can proceed to select them and retry with a reference "
        intro_txt2 += "database\nrestricted to the chosen alleles."
        layout.addWidget(QLabel(intro_txt1))
        layout.addWidget(QLabel(intro_txt2))

        self.proceed_btn1 = ProceedButton("Proceed", log=self.log, section=0)
        self.proceed_btn1.proceed.connect(self.proceed_sections)
        abort_btn = ChoiceButton("Abort")
        abort_btn.clicked.connect(self.abort)
        choice_row = QWidget()
        choice_layout = QHBoxLayout()
        choice_row.setLayout(choice_layout)
        choice_layout.addWidget(self.proceed_btn1)
        choice_layout.addWidget(abort_btn)
        layout.addWidget(choice_row)

        self.sections.append(("Try again with restricted reference database?", mywidget))

    def define_section2(self):
        """defines section 2, where the user can choose which gene system the allele belongs to
        """
        mywidget = QWidget(self)
        layout = QVBoxLayout()
        mywidget.setLayout(layout)

        intro_lbl = QLabel("Which gene family does this alle belong to?")
        intro_lbl.setStyleSheet(general.label_style_2nd)
        layout.addWidget(intro_lbl)
        self.hla_btn = ChoiceButton("HLA", self)
        kir_btn = ChoiceButton("KIR", self)
        mic_btn = ChoiceButton("MIC", self)
        mysec = ChoiceSection("Please choose a gene family", [self.hla_btn, kir_btn, mic_btn], self,
                              log=self.log)
        self.target_field = mysec.field
        layout.addWidget(mysec)

        self.proceed_btn2 = ProceedButton("Proceed", [self.target_field], self.log, 1, self)
        layout.addWidget(self.proceed_btn2)
        self.proceed_btn2.proceed.connect(self.proceed_sections)
        self.proceed_btn2.clicked.connect(self.fill_allele_table)
        self.target_field.textChanged.connect(self.proceed_btn2.check_ready)

        self.sections.append(("Choose gene family", mywidget))

    def define_section3(self):
        """defines section 3, where the user can select which alleles to add to the restricted
        reference
        """
        mywidget = QWidget(self)
        layout = QVBoxLayout()
        mywidget.setLayout(layout)

        intro1 = "Please select all alleles that should be used as a reference.\n"
        intro_lbl1 = QLabel(intro1)
        intro_lbl1.setStyleSheet(general.label_style_2nd)
        layout.addWidget(intro_lbl1)

        intro2 = "Make sure there is a suitable reference allele for each allele in your file!"
        intro2 += "\n\n(Note that only full-length alleles are available as reference alleles.)"
        layout.addWidget(QLabel(intro2))

        self.allele_table = ChoiceTableWidget(["Alleles"], [], self.log, None)
        self.allele_table.chosen_items.connect(self.catch_chosen_alleles)
        layout.addWidget(self.allele_table)

        self.proceed_btn3 = ProceedButton("Proceed", [self.allele_table.count_field], self.log, 2,
                                          self)
        self.allele_table.table.clicked.connect(self.proceed_btn3.check_ready)
        self.proceed_btn3.proceed.connect(self.proceed_sections)
        self.proceed_btn3.clicked.connect(self.fill_list_with_chosen_alleles)

        layout.addWidget(self.proceed_btn3)

        self.sections.append(("Choose reference alleles", mywidget))

    def define_section4(self):
        mywidget = QWidget(self)
        layout = QVBoxLayout()
        mywidget.setLayout(layout)

        intro = "You have chosen the following reference alleles:"
        layout.addWidget(QLabel(intro))

        self.chosen_alleles_widget = QTextEdit()
        self.chosen_alleles_widget.setFixedHeight(100)
        layout.addWidget(self.chosen_alleles_widget)

        q_text = "Do you want to create a restricted reference using only these alleles "
        q_text2 = "and reattempt uploading your allele using this reference?"
        layout.addWidget(QLabel(q_text))
        layout.addWidget(QLabel(q_text2))

        btn_widget = QWidget()
        btn_layout = QHBoxLayout()
        btn_widget.setLayout(btn_layout)
        layout.addWidget(btn_widget)

        self.proceed_btn4 = ProceedButton("Proceed", log=self.log, section=3)
        btn_layout.addWidget(self.proceed_btn4)
        self.proceed_btn4.clicked.connect(self.create_restricted_ref_and_retry)

        return_btn = QPushButton("Return to allele choice section")
        return_btn.clicked.connect(self.return_to_allele_choice)
        btn_layout.addWidget(return_btn)

        self.sections.append(("Create and run restricted reference", mywidget))

    def fill_allele_table(self):
        """populates the allele table in section 3 with allele names according to the target chosen in
        section 2
        """
        file_dic = {"HLA": "hla_allelenames.dump",
                    "MIC": "hla_allelenames.dump",
                    "KIR": "KIR_allelenames.dump"}
        target_dic = {"HLA": "hla",
                      "MIC": "hla",
                      "KIR": "kir"}

        target = self.target_field.text()
        if target not in file_dic:
            msg = f"'{target}' is not a gene system known to TypeLoader!\n"
            msg += "Please return to the previous section and choose one of the given gene systems."
            QMessageBox.warning(self, "Unknown gene system", msg)
            return None

        self.target_family = target_dic[target]
        allele_name_file = file_dic[target]
        allele_name_file = os.path.join(self.settings["root_path"],
                                        self.settings["general_dir"],
                                        self.settings["reference_dir"],
                                        allele_name_file)
        with open(allele_name_file, "rb") as f:
            allele_names = load(f)

        allele_names_restricted = [a for a in allele_names if a.startswith(target)]
        self.allele_table.table.data = allele_names_restricted
        self.allele_table.table.fill_UI()

    @pyqtSlot()
    def fill_list_with_chosen_alleles(self):
        self.log.debug(f"You have chosen {len(self.chosen_alleles)} alleles. Proceed?")
        txt = ""
        for allele in sorted(self.chosen_alleles):
            txt += f" - {allele}\n"
        self.chosen_alleles_widget.setText(txt[:-1])

    @pyqtSlot(list)
    def catch_chosen_alleles(self, items):
        self.log.debug(f"Caught {len(items)} chosen alleles: {' & '.join(items)}")
        self.chosen_alleles = items

    @pyqtSlot()
    def create_restricted_ref_and_retry(self):
        # check whether length of chosen alleles makes sense:
        if self.filetype == "XML":
            n = len(self.chosen_alleles)
            if n < 2:
                msg = f"You have chosen only {n} allele although you are trying to upload "
                msg += "an XML file.\n"
                msg += "GenDX XML files usually contain 2 alleles. "
                msg += "You need to provide a reference allele for each of those.\n\n"
                msg += "Are you sure you have selected all necessary reference alleles?"
                reply = QMessageBox.question(self, "Only 1 reference allele for XML file", msg,
                                             QMessageBox.Yes | QMessageBox.No, QMessageBox.No
                                             )
                if reply == QMessageBox.No:
                    self.proceed_btn4.setEnabled(True)
                    self.proceed_btn4.check_ready()
                    self.return_to_allele_choice()
                    return

        # create restricted db:
        from typeloader2.typeloader_core import update_reference
        self.log.debug("Deleting old restricted database, if necessary...")
        for (dirpath, subdirs, files) in os.walk(os.path.join(self.settings["temp_dir"],
                                                              "restricted_db")):
            for f in files:
                myfile = os.path.join(dirpath, f)
                os.remove(myfile)

        self.log.debug("Creating restricted database...")
        ref_path_official = os.path.join(self.settings["root_path"], self.settings["general_dir"],
                                         self.settings["reference_dir"])
        restricted_db_dir = os.path.join(self.settings["temp_dir"], "restricted_db")
        blast_path = os.path.dirname(self.settings["blast_path"])
        success, msg = update_reference.make_restricted_db(self.target_family,
                                                           ref_path_official,
                                                           self.chosen_alleles,
                                                           restricted_db_dir,
                                                           blast_path,
                                                           self.log)
        if not success:
            QMessageBox.warning(self, "Error while creating restricted db",
                                "Could not create the restricted reference database\n\n" + msg)
            self.abort()

        self.log.debug(f"Successfully created restricted db under {restricted_db_dir}")
        self.restricted_db_path.emit(restricted_db_dir)
        self.restricted_alleles.emit(self.chosen_alleles)
        self.log.debug("Closing ChooseReferenceAllelesDialog")
        self.close()

    @pyqtSlot()
    def return_to_allele_choice(self):
        self.log.debug("Returning to allele choice section...")
        self.proceed_sections(3, 2)

    @pyqtSlot()
    def abort(self):
        """if user chooses to abort, display a message, then close the dialog
        """
        self.log.info("User chose to abort ChooseReferenceAllelesDialog")
        msg = "Sorry, TypeLoader apparently can't currently handle this allele."
        msg += "\nAborting upload..."
        QMessageBox.information(self, "Can't handle divergent allele", msg)
        self.restricted_db_path.emit(None)
        self.restricted_alleles.emit([])
        self.log.info(msg)
        self.close()


# ===========================================================
# main:

if __name__ == '__main__':
    from typeloader_GUI import create_connection, close_connection
    import GUI_login

    log = general.start_log(level="DEBUG")
    log.info("<Start {}>".format(os.path.basename(__file__)))
    sys.excepthook = log_uncaught_exceptions
    mysettings = GUI_login.get_settings("admin", log)
    mydb = create_connection(log, mysettings["db_file"])

    app = QApplication(sys.argv)

    ex = NewAlleleForm(log, mydb, "20230421_ADMIN_mix_1", mysettings)
    #ex = QueryBox(log, mysettings, None)
    ex.show()

    result = app.exec_()
    close_connection(log, mydb)
    log.info("<End>")
    sys.exit(result)
