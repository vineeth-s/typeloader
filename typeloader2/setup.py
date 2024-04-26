#!/usr/bin/env python3
# -*- coding: cp1252 -*-
"""
Created on 03.09.2018

builds the TypeLoader executable for Windows using cx_Freeze

@author: schoene
"""
import os
from cx_Freeze import Executable
from setuptools import setup
import sys
from pathlib import Path
import distutils.command.build

from typeloader_installer_updater import BUILD_DIR

# remove dummy table files
for myfile in os.listdir(os.path.join(os.path.dirname(__file__), "tables")):
    if "dummy" in myfile:
        os.remove(os.path.join("tables", myfile))

build_exe_options = {
    "include_files": [
        "config_raw.ini",
        "config_base.ini",
        "config_company.ini",
        "LICENSE.txt",
        "sound_done.mp3",
        "blastn/",
        "icons/",
        "sample_files/",
        "tables/",
    ],
    "excludes": [
        "tkinter",
        "unittest",
    ]}

base = None
if sys.platform == "win32":
    base = "Win32GUI"

icon_path = os.path.join("icons", "TypeLoader.ico")


# Override build command
class BuildCommand(distutils.command.build.build):
    def initialize_options(self):
        distutils.command.build.build.initialize_options(self)
        build_dir = Path(BUILD_DIR).parent
        build_dir.mkdir(parents=True, exist_ok=True)
        self.build_base = str(build_dir)


# stdout = sys.stdout
# stderr = sys.stderr

# with open('typeloader_build.log', 'w') as log_file:
#     sys.stdout = log_file
#     sys.stderr = log_file

setup(name="TypeLoader",
      version="2.15.0",
      description="TypeLoader",
      package_dir={"": "."},
      options={"build_exe": build_exe_options},
      executables=[Executable("typeloader_GUI.pyw",
                              base=base,
                              icon=icon_path,
                              target_name="TypeLoader.exe"
                              )],
      cmdclass={"build": BuildCommand},
      )

# sys.stdout = stdout
# sys.stderr = stderr
