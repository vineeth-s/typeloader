#!/usr/bin/env python3
# -*- coding: cp1252 -*-
"""
Created on 2020-08-31

This script re-builds TypeLoader.exe and updates `typeloader_installer.nsi`,
the NIS script to create the TypeLoader Windows installer.

@author: Bianca Schoene
"""

# import modules:
import os
from collections import defaultdict
import pathlib
import subprocess
import shutil
import general

# ===========================================================
# parameters:

NEW_VERSION = "2.15.0"
BUILD_DIR = r"C:\Users\schoene\WorkFolders\Code\TL_build\exe.win-amd64-3.10"
# BUILD_DIR = r"C:\Daten\Arbeit\DKMS\TypeLoader\TL_build\exe.win-amd64-3.10"

INSTALLER_SCRIPT = "typeloader_installer.nsi"
INSTALLER_SCRIPT_NEW = "typeloader_installer_new.nsi"
NSIS_PATH = r"C:\Program Files (x86)\NSIS\makensis.exe"
IGNORE_FILES = [pathlib.Path(BUILD_DIR, "config_base.ini"),
                pathlib.Path(BUILD_DIR, "config_company.ini"),
                ]
UPDATE_INSTALLER = True
BUILD_INSTALLER = True

# ===========================================================
# functions:


def remove_subdir(subdir_name: str, parent_path: pathlib.Path, log):
    """Remove subdir and all files in it from parent_path.
    """
    dir_path = parent_path / subdir_name
    if dir_path.exists():
        log.debug(f"\t- removing subfolder {subdir_name} from {parent_path}...")
        shutil.rmtree(dir_path)
    assert not dir_path.exists()


def cleanup_build_dir(log):
    """Remove redundant directories (artifacts of typeloader being a package)
    """
    log.info("Cleaning up build dir...")
    global BUILD_DIR

    BUILD_DIR = pathlib.Path(BUILD_DIR)
    cleanup_dir = {BUILD_DIR / "lib/typeloader2": [
        ".idlerc",
        ".idea",
        "flowcharts",
        "icons",
        "installer",
        "logs",
        "reference_data2",
        "sample_files",
        "tables",
        "temp",
        "tests",
    ],
        BUILD_DIR / "lib": [
            "authuser",
            "typeloader_core"
        ],
        BUILD_DIR: [
            "ENA_Webin_CLI",
        ],
    }
    for parent_path in cleanup_dir:
        for sub_dir in cleanup_dir[parent_path]:
            remove_subdir(sub_dir, parent_path, log)

    log.info("\t=> done!")


def rebuild_exe(log):
    """renames old BUILD_DIR, if it exists,
    then calls setup.py via poetry to create a fresh build
    """
    log.info("Re-building TypeLoader.exe...")

    # rename old build dir:
    build_dir = pathlib.Path(BUILD_DIR)
    build_dir_old = pathlib.Path(f"{BUILD_DIR}_old")
    if build_dir.exists():
        if build_dir_old.exists():
            shutil.rmtree(build_dir_old)
        log.info("\tRenaming old build-dir...")
        build_dir.rename(build_dir_old)

    # re-build:
    log.info("\tRe-building...")
    ret_code = subprocess.call(["poetry", "run", "python", "setup.py", "build"], shell=True)
    if ret_code == 0:
        log.info("\t\t=> Success!")

        cleanup_build_dir(log)

        return True
    else:
        log.error("\t\tBuild failed!")
        return False


def gather_files(build_dir, log):
    """collects all files in the current BUILD_DIR,
    returns them as a defaultdict[filepath] = True (False as default value)
    """
    log.info(f"Searching for files in {build_dir}...")
    found_files = defaultdict(lambda: False)
    for root, dirs, files in os.walk(build_dir):
        for filename in files:
            myfile = pathlib.Path(root, filename)
            found_files[myfile] = True
    log.info(f"\t=> found {len(found_files)} files")
    return found_files


def read_script(nsi_file, log):
    """reads the existing .nsi script and retrieves which files are already covered
    returns them as a defaultdict[filepath] = True (False as default value)
    """
    log.info(f"Reading {nsi_file}...")
    with open(nsi_file, "r") as f:
        files = defaultdict(lambda: False)
        for line in f:
            if line.startswith("  File"):
                myfile = pathlib.Path(line.split('"')[1])
                files[myfile] = True

    log.info(f"\t=> found {len(files)} files")
    return files


def check_all_files_contained(files_found, files_listed, log):
    """checks files contained in the current BUILD_DIR, whether they are already contained in the installer;
    returns a list of all files that are not. (These are new in this build and were not contained in the previous one.)
    """
    log.info(f"Checking for files missing in {INSTALLER_SCRIPT}...")
    missing_files = defaultdict(list)
    nr_missing_files = 0
    for file in files_found:
        if not files_listed[file]:
            if file not in IGNORE_FILES:
                mydir = file.parent
                filename = file.name
                missing_files[mydir].append(filename)
                nr_missing_files += 1
                log.debug(f"File {file} is not contained in {INSTALLER_SCRIPT}!")
    log.info(f"\t=> found {nr_missing_files} missing files")
    return missing_files


def check_all_files_found(files_found, files_listed, log):
    """checks files listed in the installer, whether they are contained in the current BUILD_DIR;
    returns a list of all files that are not. (These are leftovers from previous builds.)
    """
    log.info(f"Checking for files missing in {BUILD_DIR}...")
    missing_files = defaultdict(lambda: False)
    for file in files_listed:
        if not files_found[file]:
            if file not in IGNORE_FILES:
                missing_files[file] = True
                log.debug(f"File {file} is not contained in {BUILD_DIR} => probably deprecated!")
    log.info(f"\t=> found {len(missing_files)} deprecated files")
    return missing_files


def check_old_script_for_consistency_with_new_files(log):
    """checks consistency between old installer-script and new build-dir;
    returns list of files missing from the installer and list of files in the installer that are not contained
    in the build-dir (each as a list of pathlib.Path objects)
    """
    log.info(f"Checking {INSTALLER_SCRIPT} for consistency with {BUILD_DIR}...")
    files_found = gather_files(BUILD_DIR, log)
    files_listed = read_script(INSTALLER_SCRIPT, log)

    missing_files = check_all_files_contained(files_found, files_listed, log)
    deprecated_files = check_all_files_found(files_found, files_listed, log)
    return missing_files, deprecated_files


def adjust_installer(missing_files, deprecated_files, log):
    """writes a new installer file with changes as needed
        - add missing files in appropriate places
        - remove adding of deprecated files
        - DOES NOT remove delete-lines for deprecated files (these can be leftovers from previous versions!)
        - adjusts version number
    """
    log.info(f"Writing new installer to {INSTALLER_SCRIPT_NEW}...")

    build_dir = str(pathlib.Path(BUILD_DIR))
    ena_webin_dir = pathlib.Path(BUILD_DIR, "ENA_Webin_CLI")
    removed_files = 0
    added_files = 0
    delete_me = []

    with open(INSTALLER_SCRIPT, "r") as f, open(INSTALLER_SCRIPT_NEW, "w") as g:
        section = "header"
        subsection = None
        new_webin_cli = False

        for line in f:
            myline = line
            # parse sections:
            if line.startswith('Section'):
                if line.startswith('Section "MainSection" SEC01'):
                    section = "SEC01"  # copy files
                elif line.startswith('Section -AdditionalIcons'):
                    section = "rest"
                elif line.startswith("Section Uninstall"):
                    section = "Uninstall"

            # parse individual lines:
            if section == "header":
                if line.startswith("!define PRODUCT_VERSION"):
                    myline = f'!define PRODUCT_VERSION "{NEW_VERSION}"\n'

            elif section == "SEC01":
                # add missing files to appropriate section if it exists:
                if line.startswith("  SetOutPath"):
                    mydir = line.split('"')[1]  # output dir
                    subsection = mydir
                    target_dir = pathlib.Path(mydir.replace("$INSTDIR", build_dir))

                    if target_dir in missing_files:
                        if target_dir == ena_webin_dir:
                            new_webin_cli = True
                        for file in missing_files[target_dir]:
                            myfile = pathlib.Path(target_dir, file)
                            newline = f'  File "{myfile}"\n'
                            myline += newline
                            added_files += 1
                            delete_me.append(myfile)

                        missing_files.pop(target_dir)

                elif line.startswith("  File"):
                    myfile = pathlib.Path(line.split('"')[1])
                    if deprecated_files[myfile]:
                        myline = ""
                        removed_files += 1

                    if subsection == r"$INSTDIR\ENA_Webin_CLI":
                        if new_webin_cli:  # delete old webin CLI version
                            myline = line.replace(f'File "{build_dir}', 'Delete "$INSTDIR')

                # add remaining missing_files, if necessary:
                elif line.strip() == r'CreateDirectory "$SMPROGRAMS\typeloader"':
                    if missing_files:
                        myline = ""
                        for mydir in missing_files:
                            target_dir = f"{mydir}".replace(build_dir, "$INSTDIR")
                            target_dir = pathlib.Path(target_dir)
                            myline += f'  SetOutPath "{target_dir}"\n'
                            for myfile in missing_files[mydir]:
                                myline += f'  File "{pathlib.Path(mydir, myfile)}"\n'
                        myline += line

            elif section == "Uninstall":
                if line.startswith('  RMDir'):
                    if delete_me:
                        myline = ""
                        for file in delete_me:
                            myfile = f"{file}".replace(build_dir, "$INSTDIR")
                            myline += f'  Delete "{myfile}"\n'
                            delete_me.pop(0)
                        myline += "\n" + line

            g.write(myline)

    log.info(f"\t=> added {added_files} new file(s)")
    log.info(f"\t=> removed {removed_files} deprecated file(s)")


def compile_installer(log) -> bool:
    """compiles the new installer script into a new Setup.exe using NSIS.

    return success (bool)
    """
    log.info("Compiling installer with NSIS...")
    if not os.path.isfile(NSIS_PATH):
        log.warning(f"\t=> Could not find NSIS under {NSIS_PATH}. Please adjust path and run again.")
        return

    cmd = [NSIS_PATH, "/NOTIFYHWND", "134742", os.path.abspath(INSTALLER_SCRIPT)]
    result_code = subprocess.call(cmd, shell=True)
    if result_code == 0:
        log.info("\t\t=> Success!")
        return True
    else:
        log.info("\t\t=> Failed!")
        return False


def wrap_up_and_compile_installer(changes, log):
    """ties up all loose ends and compiles the new installer into a Setup.exe
    """
    log.info("Wrapping up...")
    if pathlib.Path(INSTALLER_SCRIPT_NEW).exists():
        log.info("\tReplacing installer-script with updated version...")
        pathlib.Path(INSTALLER_SCRIPT).unlink()
        pathlib.Path(INSTALLER_SCRIPT_NEW).rename(INSTALLER_SCRIPT)

    build_dir_old = pathlib.Path(f"{BUILD_DIR}_old")
    if build_dir_old.exists():
        log.info("\tRemoving old build-dir...")
        shutil.rmtree(build_dir_old, ignore_errors=True)

    if changes:
        log.info("Successfully updated!")
    else:
        log.info("No changes were necessary! The old installer script is still good to go. :-)")

    success = compile_installer(log)

    if success:
        new_name = f"TypeLoader_Setup_V{NEW_VERSION}.exe"
        try:
            pathlib.Path(new_name).unlink(missing_ok=True)
            os.rename("TypeLoader_Setup.exe", new_name)
        except FileNotFoundError:
            pass

        log.info(f"The new installer can be found in this directory as {new_name}. \n"
                 "Please test it thoroughly before deploying!")

    else:
        log.info("No installer could be built")


# ===========================================================
# main:

def main(log):
    rebuild = general.confirm("Should I completely rebuild the TypeLoader.exe?")
    if rebuild:
        success = rebuild_exe(log)
    else:
        log.info("\t=> User chose not to re-build")
        success = True

    if not success:
        return

    changes = False
    if UPDATE_INSTALLER:
        missing_files, deprecated_files, = check_old_script_for_consistency_with_new_files(log)
        if missing_files or deprecated_files:
            changes = True
        adjust_installer(missing_files, deprecated_files, log)

    if BUILD_INSTALLER:
        wrap_up_and_compile_installer(changes, log)

    else:
        log.info("Parameter said not to update the installer script => leaving it as it was")

    general.play_sound(log)


if __name__ == "__main__":
    logger = general.start_log(level="DEBUG")
    logger.info("<Start>")
    main(logger)
    logger.info("<End>")
