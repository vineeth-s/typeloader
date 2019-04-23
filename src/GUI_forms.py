#!/usr/bin/env python3
# -*- coding: cp1252 -*-
'''
Created on 13.03.2018

GUI_forms.py

components for forms and dialogs

@author: Bianca Schoene
'''

# import modules:

import sys, os
from PyQt5.QtSql import QSqlQueryModel, QSqlQuery
from PyQt5.QtWidgets import (QApplication, QListView, QGroupBox, QGridLayout, 
                             QFileDialog, QLabel, QPushButton, QLineEdit, QDialog, 
                             QTreeWidgetItem, QVBoxLayout, QHBoxLayout, QFrame,
                             QTableWidget, QWidget, QCheckBox, QTableWidgetItem)
from PyQt5.Qt import pyqtSlot, pyqtSignal, QTreeWidget
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor

import general, db_internal
from GUI_forms_new_project import NewProjectForm

#===========================================================
# parameters:

from __init__ import __version__
#===========================================================
# classes:
 
class QueryDialog(QDialog):
    """a dialog to choose an item from a query
    """
    choice = pyqtSignal(str)
     
    def __init__(self, query):
        super().__init__()
        self.query = query
        self.create_model()
        self.init_UI()
         
    def create_model(self):
        """creates the model as QSqlQueryModel,
        using the given query
        """
        self.model = QSqlQueryModel()
        q = QSqlQuery()
        q.exec_(self.query)
        self.model.setQuery(q)
        
         
    def init_UI(self):
        """setup the UI
        """
        layout = QVBoxLayout()
        self.setLayout(layout)
        self.resize(200,200)
        self.title = "Choose an existing project"
         
        self.list = QListView(self)
        layout.addWidget(self.list)
        self.list.setModel(self.model)
        self.list.setWhatsThis("Choose a project by clicking on it")
         
        self.btn = QPushButton("Accept", self)
        layout.addWidget(self.btn)
        self.btn.clicked.connect(self.on_btn_clicked)
        self.btn.setWhatsThis("Click here to accept your selection (works only if a project has been selected)")
     
    def on_btn_clicked(self):
        """when self.btn is clicked, accept the choice and emit it as self.choice
        """
        selected = self.list.selectedIndexes()
        if selected:
            index = selected[0]
            chosen = self.model.data(index, Qt.DisplayRole)
            self.choice.emit(chosen)
            self.close()
        
        self.choice.emit("")
        self.close()
 
 
class SectionExpandButton(QPushButton):
    """a QPushbutton that can expand or collapse its section
    """
    def __init__(self, item, text = "", parent = None):
        super().__init__(text, parent)
        self.section = item
        self.clicked.connect(self.on_clicked)
         
    def on_clicked(self):
        """toggle expand/collapse of section by clicking
        """
        if self.section.isExpanded():
            self.section.setExpanded(False)
        else:
            self.section.setExpanded(True)
             
 
class ProceedButton(QPushButton):
    """a QPushButton that checks whether all necessary data has been given to a section
    & emits number of next section when it is ok to proceed;
    all fields that need to be checked should be given as list <items>;
    => they need to contain text to be accepted;
    section is the section number this Button belongs to
    """
    proceed = pyqtSignal(int)
     
    def __init__(self, text = "", items = [], log = None, section = 0, parent = None,
                 only1 = False):
        """constructor
        """
        super().__init__(text, parent)
        self.text = text
        self.log = log
        self.items = items
        self.section = section
        self.setEnabled(False)
        self.setCheckable(True)
        self.clicked.connect(self.on_clicked)
        self.only1 = only1
        self.check_ready()
         
    @pyqtSlot(str)
    def enable(self, sig):
        """enables the button if a non-False signal is passed
        """
        self.log.debug("Enabling button {}".format(self.text))
        if sig:
            self.setEnabled(True)
         
    @pyqtSlot()
    def check_ready(self, debugging = False):
        """checks all items: if all evaluate to True, enables proceeding 
        """
        self.log.debug("Ready for proceeding?")
        ready = True
        active_fields = 0
        for item in self.items:
            item_type = str(type(item))
            if "QTableWidgetItem" in item_type:
                active_fields += 1
                text = item.text()
                try:
                    text = int(text)
                except ValueError:
                    pass
                if not text:
                    ready = False
                if debugging:
                    print("text:", text)
            elif "QCheckBox" in item_type:
                if item.isChecked():
                    active_fields += 1
                if debugging:
                    print("checkbox:", item.isChecked())
            else:
                if item.isEnabled():
                    active_fields += 1
                    if "QTextEdit" in item_type:
                        text = item.toPlainText()
                    else:
                        text = item.text()
                    if not text:
                        ready = False
                    if debugging:
                        print("text:", text)
        
        if active_fields == 0: # if nothing selected
            ready = False
        elif self.only1: # if only one should be selected
            if active_fields != 1:
                ready = False
        
        if ready:
            self.log.debug("\tReady!")
            self.setEnabled(True)
            self.setStyleSheet(general.btn_style_ready)
        else:
            self.log.debug("\tNot ready!")
            self.setDisabled(True)
            self.setStyleSheet(general.btn_style_normal)
             
    @pyqtSlot(str)
    def change_to_normal(self, _ = None):
        """sets button to normal look once it has been used
        (accepts an optional unused string so it can be connected to signals emitting a string)
        """
        self.setStyleSheet(general.btn_style_normal)
        self.setChecked(False)
         
    def on_clicked(self):
        """when button is clicked,
        emit whether section is ready for proceeding
        """
        self.log.debug("Proceed to section {}!".format(self.section + 1))
        self.check_ready()
        self.setChecked(True)
        self.proceed.emit(self.section + 1)
         
 
class ChoiceButton(QPushButton):
    """a QPushButton to open a dialog, select something with it, 
    and emit the chosen item as signal done;
    reimplement the open_dialog() method for individual kinds of ChoiceButtons
    """
    done = pyqtSignal(str)
     
    def __init__(self, text, parent=None):
        """constructor
        """
        self.text = text
        super().__init__(text, parent)
        self.setStyleSheet(general.btn_style_clickme)
        self.clicked.connect(self.open_dialog)
         
    @pyqtSlot()
    def open_dialog(self):
        """reimplement for individual buttons
        """
        pass
         
    @pyqtSlot(str)
    def change_to_normal(self, _ = None):
        """sets button to normal look once it has been used
        (accepts an optional unused argument so it can be connected to signals emitting a single object)
        """
        self.setStyleSheet(general.btn_style_normal)
         
    @pyqtSlot(str)
    def emit_choice(self, choice):
        """emits choice made with dialog
        """
        if choice:
            self.done.emit(choice)
            self.change_to_normal("")
         
         
class FileButton(ChoiceButton):
    """a QPushButton to select and a file;
    triggers a FileDialog and emits the name of the chosen file
    """
    def __init__(self, text, default_path = "/home", parent=None, log = None):
        """constructor
        """
        self.default_path = default_path
        super().__init__(text, parent)
        self.log = log
         
    @pyqtSlot()
    def open_dialog(self):
        """opens a file selection dialog to choose a file,
        emits chosen file as signal file_chosen
        """
        try:
            myfiles = QFileDialog.getOpenFileName(self, self.text, self.default_path) # TODO: (future) filter to file type
            myfile = myfiles[0] # only return first file
            if myfile:
                self.change_to_normal("")
                self.emit_choice(myfile)
        except Exception as E:
            if self.log:
                self.log.exception(E)
                
     
     
class QueryButton(ChoiceButton):
    """a QPushButton to select something from a query;
    triggers a Dialog and emits the name of the chosen item
    """
    def __init__(self, text, query = "", parent=None):
        """constructor
        """
        self.query = query
        super().__init__(text, parent)
         
    @pyqtSlot()
    def open_dialog(self):
        """opens a QueryDialog to choose an item from a QListView,
        passes chosen item to self.emit_choice()
        """
        self.dialog = QueryDialog(self.query)
        self.dialog.show()
        self.dialog.choice.connect(self.emit_choice)
         
 
class NewProjectButton(ChoiceButton):
    """a QPushButton that opens a NewProjectDialog to create a new project
    and emits the chosen project as signal <done>
    """
    def __init__(self, text, log, mydb, settings, parent=None):
        """constructor
        """
        super().__init__(text, parent)
        self.log = log
        self.mydb = mydb
        self.settings = settings
         
    @pyqtSlot()
    def open_dialog(self):
        """opens a NewProjectDialog to create a new project;
        passes chosen project to self.emit_choice()
        """
        self.dialog = NewProjectForm(self.log, self.mydb, self.settings)
        self.dialog.show()
        self.dialog.project_changed.connect(self.emit_choice)
         
     
class CollapsibleDialog(QDialog):
    """a dialog to which collapsible sections can be added;
    reimplement define_sections() to define sections and
        add them as (title, widget) tuples to self.sections
     
    reimplemented from http://www.fancyaddress.com/blog/qt-2/create-something-like-the-widget-box-as-in-the-qt-designer/
    """
    def __init__(self, parent = None):
        super().__init__(parent)
        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        layout = QVBoxLayout()
        layout.addWidget(self.tree)
        self.setLayout(layout)
        self.tree.setIndentation(0)
        self.sections = []
        self.section_dic = {}
        self.define_sections()
        self.add_sections()
         
    def add_sections(self):
        """adds a collapsible sections for every 
        (title, widget) tuple in self.sections
        """
        for (i, (title, widget)) in enumerate(self.sections):
            button = self.add_button(title)
            section = self.add_widget(button, widget)
            button.addChild(section)
            if i == 0:
                button.setExpanded(True)
            self.section_dic[i] = (button, section)
 
    def define_sections(self):
        """reimplement this to define all your sections
        and add them as (title, widget) tuples to self.sections
        """
        widget = QFrame(self.tree)
        layout = QHBoxLayout(widget)
        layout.addWidget(QLabel("Bla"))
        layout.addWidget(QLabel("Blubb"))
        title = "Section 1"
        self.sections.append((title, widget))
 
    def add_button(self, title):
        """creates a QTreeWidgetItem containing a button 
        to expand or collapse its section
        """
        item = QTreeWidgetItem()
        self.tree.addTopLevelItem(item)
        self.tree.setItemWidget(item, 0, SectionExpandButton(item, text = title))
        return item
 
    def add_widget(self, button, widget):
        """creates a QWidgetItem containing the widget,
        as child of the button-QWidgetItem
        """
        section = QTreeWidgetItem(button)
        section.setDisabled(True)
        self.tree.setItemWidget(section, 0, widget)
        return section
     
    @pyqtSlot(int, int)
    def proceed_sections(self, old_section, new_section):
        """collapses old section and expands next section
        """
        (button_new, _) = self.section_dic[new_section]
        (button_old, _) = self.section_dic[old_section]
        button_new.setExpanded(True)
        button_old.setExpanded(False)
        try:
            self.sender().setChecked(False) # if sent from a button, un-press it
        except:
            pass
         
     
class ChoiceSection(QGroupBox):
    """a widget containing a read-only-field and several buttons,
    any of which can be used to make the choice, 
    which is displayed in the field & emitted as signal choice 
    """
    choice = pyqtSignal(str)
    def __init__(self, field_text, btns, parent=None, label_width = None):
        super().__init__(parent)
        self.field_text = field_text
        self.buttons = btns
        self.label_width = label_width
        self.init_UI()
        #TODO: add log to log choice
         
    def init_UI(self):
        grid = QGridLayout()
        self.setLayout(grid)
        
        label = QLabel(self.field_text, self)
        if self.label_width:
            label.setMinimumWidth(self.label_width)
            label.setMaximumWidth(self.label_width)
        grid.addWidget(label, 0, 0)
        self.field = QLineEdit(self)
        self.field.setReadOnly(True)
        self.field.setStyleSheet(general.label_style_entry)
        grid.addWidget(self.field, 0, 1)
         
        self.button_dic = {}
        row = 1
         
        for (i, button) in enumerate(self.buttons):
            row += 1
            grid.addWidget(button, row, 0, 1, 2)
            self.button_dic[i] = button
         
            if i < len(self.buttons) - 1:
                row += 1
                or_lbl = QLabel("or", self)
                or_lbl.setAlignment(Qt.AlignCenter)
                grid.addWidget(or_lbl, row, 0, 1, 2)
         
        for i in self.button_dic:
            btn = self.button_dic[i]
            btn.done.connect(self.field.setText)
            btn.done.connect(self.choice.emit)
            #TODO: log choice
            for j in self.button_dic:
                if i != j:
                    btn2 = self.button_dic[j]
                    btn.done.connect(btn2.change_to_normal)
     
    @pyqtSlot()
    def reactivate(self):
        """returns buttons to clickme-style if they need to be re-used
        """
        for i in self.button_dic:
            btn = self.button_dic[i]
            btn.setStyleSheet(general.btn_style_clickme)


class FileChoiceTable(QTableWidget):
    """displays samples so user can choose some of them
    """
    files = pyqtSignal(int)
    files_chosen = pyqtSignal(int)
    
    def __init__(self, project, log, header, query, num_columns, myfilter, 
                 allele_status_column = None, instant_accept_status = None, 
                 parent = None):
        super().__init__()
        self.log = log
        if parent:
            self.settings = parent.settings
        else:
            import GUI_login
            self.settings = GUI_login.get_settings("admin", self.log)
        self.header = header
        self.allele_status_column = allele_status_column
        self.instant_accept_status = instant_accept_status
        self.query = query
        self.num_columns = num_columns
        self.myfilter = myfilter
        self.keep_choices = False
        self.check_dic = {}
        self.init_UI()
    
    def init_UI(self):
        self.setColumnCount(len(self.header))
        self.setHorizontalHeaderLabels(self.header)
        self.horizontalHeader().setStretchLastSection(True)
        self.verticalHeader().hide()
    
    def reset_filter(self, myfilter = ""):
        """sets a new filter for use in self.get_data()
        """
        self.myfilter = myfilter
        
    def get_data(self):
        """get alleles from database
        """
        success, data = db_internal.execute_query(self.query + self.myfilter, self.num_columns, 
                                                  self.log, "retrieving data for FileChoiceTable from database", 
                                                  "Database error", self)
        if success:
            self.data = data
        self.log.debug("Emitting 'files = {}'".format(len(self.data)))
        self.files.emit(len(self.data))
    
    def fill_UI(self):
        """fills table with data
        """
        self.get_data()
        
        if self.keep_choices: # if table is refreshed just to update data before re-attempt, don't repopulate
            self.count_chosen()
            return
        
        rows = len(self.data) + 1
        self.setRowCount(rows)
        
        all_checked = True
        for (i, row) in enumerate(self.data):
            cell_widget = QWidget()
            mini_layout = QHBoxLayout(cell_widget)
            cell_widget.setLayout(mini_layout)
            self.check_dic[i] = QCheckBox(self)
            self.check_dic[i].clicked.connect(self.count_chosen)
            self.check_dic[i].clicked.connect(self.unselect_select_all)
            mini_layout.addWidget(self.check_dic[i])
            mini_layout.setAlignment(Qt.AlignCenter)
            self.setCellWidget(i, 0, cell_widget)
            for (k, item) in enumerate(row):
                self.setItem(i, k + 1, QTableWidgetItem(str(item)))
            if self.allele_status_column:
                status = row[self.allele_status_column]
                color = general.color_dic[general.allele_status_dic[status.lower()]]
                status_item = QTableWidgetItem(status)
                status_item.setBackground(QColor(color))
                if self.instant_accept_status:
                    if status == self.instant_accept_status:
                        self.check_dic[i].setChecked(True)
                        status_item.setFont(general.font_bold)
                    else:
                        all_checked = False
            self.setItem(i, self.allele_status_column + 1, status_item)
        
        # add select-all row:
        cell_widget = QWidget()
        mini_layout = QHBoxLayout(cell_widget)
        cell_widget.setLayout(mini_layout)
        self.check_all = QCheckBox(self)
        mini_layout.addWidget(self.check_all)
        if all_checked:
            self.check_all.setChecked(True)
        mini_layout.setAlignment(Qt.AlignCenter)
        self.check_all.clicked.connect(self.toggle_select_all)
        self.setCellWidget(rows-1, 0, cell_widget)
        self.setItem(rows-1, 1, QTableWidgetItem(""))
        self.setItem(rows-1, 2, QTableWidgetItem("Select All"))
        self.setItem(rows-1, 3, QTableWidgetItem(""))       
        self.setItem(rows-1, 4, QTableWidgetItem(""))
        
        self.resizeColumnsToContents()
        self.count_chosen()
    
    def count_chosen(self):
        """counts and emits the number of currently chosen files
        """
        self.log.debug("Recounting chosen files...")
        n = 0
        for i in self.check_dic:
            box = self.check_dic[i]
            if box.checkState():
                n += 1
        self.log.debug("\t=> Currently {} files chosen".format(n))
        self.files_chosen.emit(n)
    
    def unselect_select_all(self):
        """unchecks the 'select all' checkbox when a single file is unchecked manually
        """
        if not self.sender().checkState():
            self.check_all.setChecked(False)
        
    def toggle_select_all(self):
        """select or deselect all alleles at once
        """
        if self.check_all.checkState():
            self.log.debug("Selecting all alleles")
            for i in self.check_dic:
                box = self.check_dic[i]
                box.setChecked(True)
                self.files_chosen.emit(len(self.data))
        else:
            self.log.debug("Deselecting all alleles")
            for i in self.check_dic:
                box = self.check_dic[i]
                box.setChecked(False)
                self.files_chosen.emit(0)

pass
#===========================================================
# main:
        
if __name__ == '__main__':
    from typeloader_GUI import create_connection, close_connection
    import GUI_login
    log = general.start_log(level="DEBUG")
    log.info("<Start {} V{}>".format(os.path.basename(__file__), __version__))
    settings_dic = GUI_login.get_settings("admin", log)
    mydb = create_connection(log, settings_dic["db_file"])
    
    app = QApplication(sys.argv)
    ex = QueryDialog("select project_name from projects")
    ex.show()
    
    result = app.exec_()
    close_connection(log, mydb)
    log.info("<End>")
    sys.exit(result)

