#!/usr/bin/env python

"""
update_referce.py

handles updating of the reference data for TypeLoader
"""
import os
import zipfile
from pathlib import Path
import datetime
import re, subprocess, shutil, socket
import urllib.request
from urllib.error import URLError
import hashlib
from typing import Tuple
import logging

if __name__ == "__main__":
    import hla_embl_parser
else:
    from . import hla_embl_parser

remote_db_path = {
    "hla_path": "https://github.com/DKMS-LSL/IMGTHLA2/raw/Latest/hla.dat.zip",
    "kir_path": "https://github.com/DKMS-LSL/IPDKIR/raw/Latest/KIR.dat"
}

remote_db_path_old_files = {
    "hla_path": "https://github.com/ANHIG/IMGTHLA/raw/VERSION/hla.dat.zip",
    "kir_path": "https://github.com/ANHIG/IPDKIR/raw/VERSION/KIR.dat"
}

remote_checksumfile_index = {
    "hla_checksums_file": "https://raw.githubusercontent.com/DKMS-LSL/IMGTHLA2/Latest/md5checksum.txt",
    "kir_checksums_file": "https://raw.githubusercontent.com/DKMS-LSL/IPDKIR/Latest/md5checksum.txt"
}

COUNTRY_URL = "https://raw.githubusercontent.com/DKMS-LSL/typeloader_reference_parser/master/data/countries.csv"
COUNTRY_FILE = "collection_country_options.csv"


def read_remote_file(myurl, proxy, timeout, log, return_binary=False, to_file=None):
    """reads a remote file from a given URL, either using the given proxy or not if none is given,
    returns the data as string
    """
    if proxy:
        log.debug("Using proxy...")
        proxy_handler = urllib.request.ProxyHandler({"https": proxy})
        opener = urllib.request.build_opener(proxy_handler)
    else:
        log.debug("Not using proxy...")
        opener = urllib.request.build_opener()

    with opener.open(myurl, timeout=timeout) as request:
        if return_binary:
            if to_file:
                with open(to_file, 'wb') as g:
                    shutil.copyfileobj(request, g)
                return True
            else:
                return request.read()
        else:
            html = request.read()
            data = html.decode("UTF-8", "ignore")
            return data


def read_local_md5_checkfile(ref_path: str, db_name: str, log) -> Tuple[str | None, bool]:
    """Get md5 checksum and last modified date of local reference file.

    Returns Tuple of (md5, modified_today):
        - md5 = md5 of current .dat file (None if checksum-file not found)
        - modified_today = True if modified today, else False
    """
    log.info("\tGetting MD5 checksum of local file...")
    curr_md5_file = os.path.join(ref_path, f"curr_md5_{db_name}.txt")
    if not os.path.exists(curr_md5_file):
        return None, False

    with open(curr_md5_file, "r") as f:
        line = f.read().strip()
        [md5, timestamp] = line.split()
        log.debug(f"\t\t=> MD5 of local file: {md5}")
        log.debug(f"\t\t=> last updated: {timestamp}")

    today = datetime.datetime.now().strftime('%d.%m.%y')
    if timestamp == today:
        log.debug("\t\t=> file already updated today")
        return md5, True

    return md5, False


def get_remote_md5checksum(db_name, IPD_db_name, proxy, log):
    """retrieves MD5 checksum from IPD's checksum_file
    """
    log.debug("\tGetting checksum of current remote file...")
    remote_checksumfile = remote_checksumfile_index["%s_checksums_file" % db_name]
    log.info(str(remote_checksumfile))

    checksum_data = read_remote_file(remote_checksumfile, proxy, 10, log)

    # The checksum_data are in lines of the form MD5 (hla.dat) = 2dde3a26abf52c11a70aae7fa8f14666\n
    for pattern in [f".*{IPD_db_name}.dat.zip\) = (.*)\n", f".*{IPD_db_name}.dat\) = (.*)\n"]:
        datFile_regex = re.compile(pattern)
        match = datFile_regex.search(checksum_data)
        if match:
            md5 = match.groups()[0]
            log.debug("\t=> {}".format(md5))
            return md5

    # try different pattern because IPD is stupidly inconsistent:
    for pattern2 in [f"\w+\s+ *{IPD_db_name}.dat.zip\n", f"\w+\s+ *{IPD_db_name}.dat\n"]:
        datFile_regex2 = re.compile(pattern2)
        match = datFile_regex2.search(checksum_data)
        if match:
            md5 = match.group().split()[0]
            log.debug("\t=> {}".format(md5))
            return md5

    md5 = None
    msg = "Could not find MD5 checksum of remote {} file!".format(IPD_db_name)
    log.error(msg)
    return md5


def get_local_md5checksum(local_reference_file, log):
    """gets the MD5 checksum of the current local database version
    """
    log.debug("\tGetting checksum of local file {}...".format(local_reference_file))

    md5 = hashlib.md5(open(local_reference_file, "rb").read()).hexdigest()
    log.debug("\t=> {}".format(md5))
    return md5


def make_blast_db(target, ref_dir, blast_path, log):
    """calls blast to create the local blast database for TypeLoader;
    returns success (=BOOL), msg (=String if error; None if not)
    """
    log.debug("\tCreating local blast database...")
    makeblastdb = os.path.join(blast_path, "makeblastdb")
    cmd_list = [makeblastdb, "-dbtype", "nucl", "-in", "{}".format(os.path.join(ref_dir, "parsed{}.fa".format(target)))]
    try:
        subprocess.run(cmd_list, check=True, shell=False)
        return True, None

    except Exception as E:
        log.exception(E)
        msg = "Could not create local blast database:\n{}\n".format(repr(E))
        cmd = " ".join(cmd_list)
        msg += cmd
        log.debug(cmd)
        return False, msg


def move_files(ref_path_temp, ref_path, target, log):
    """moves all files from ref_path_temp to ref_path, replacing existing files
    """
    log.debug("\tReplacing old files with new files...")
    for myfile in os.listdir(ref_path_temp):
        if target.lower() in myfile.lower():
            src_path = os.path.join(ref_path_temp, myfile)
            target_path = os.path.join(ref_path, myfile)
            log.debug("\t\t- {}".format(myfile))
            if os.path.exists(target_path):
                os.remove(target_path)
            shutil.move(src_path, target_path)


def check_database(db_name, reference_local_path, proxy, log, skip_if_updated_today=True):
    """checks IPD GitHub account for new releases
    """
    log.info("Checking {} for IPD update...".format(db_name.upper()))
    if db_name == "kir":
        use_dbname = "KIR"  # biological databases and consistency in naming are arch enemies
    else:
        use_dbname = db_name

    local_reference_file = os.path.join(reference_local_path, "{}.dat".format(use_dbname))
    local_md5, updated_today = read_local_md5_checkfile(reference_local_path, use_dbname, log)
    if os.path.isfile(local_reference_file):
        if skip_if_updated_today:
            if updated_today:
                return False, "Reference file was already updated today"

        remote_md5 = get_remote_md5checksum(db_name, use_dbname, proxy, log)
        if not remote_md5:
            log.info("Aborting attempt to update the {} reference".format(db_name))
            return False, f"Could not reach or parse IPD's checksum file for {db_name}!"

        if remote_md5 == local_md5:
            msg = "\t=> {} is up to date".format(db_name)
            log.info(msg)
            return False, msg

    else:
        log.info("=> No local reference file found")
        local_md5 = None  # if reference file does not exist, download it

    msg = "=> Found new reference version for {}".format(db_name)
    log.info(msg)
    return True, msg


def update_database(db_name, reference_local_path, blast_path, proxy, log, version=None):
    """updates a reference database
    """
    log.info("Retrieving new database version for {}...".format(db_name))
    if db_name == "kir":
        use_dbname = "KIR"  # biological databases and consistency in naming are arch enemies
        is_zipped = False
    else:
        use_dbname = db_name
        is_zipped = True

    ref_path_temp = os.path.join(reference_local_path, "temp")
    os.makedirs(ref_path_temp, exist_ok=True)

    if version:
        remote_db_file = remote_db_path_old_files["%s_path" % db_name].replace("VERSION", version)
    else:
        remote_db_file = remote_db_path["%s_path" % db_name]
    log.debug(f"\tdownloading new file from {remote_db_file}...")
    local_db_file = os.path.join(ref_path_temp, "%s.dat" % use_dbname)

    try:
        if is_zipped:
            local_db_file_zipped = local_db_file + ".zip"
            read_remote_file(remote_db_file, proxy, 60, log, return_binary=True, to_file=local_db_file_zipped)
            local_md5 = get_local_md5checksum(local_db_file_zipped, log)
            with zipfile.ZipFile(local_db_file_zipped) as zf:
                zf.extractall(ref_path_temp)
        else:
            db_response = read_remote_file(remote_db_file, proxy, 60, log, return_binary=True)
            with open(local_db_file, "wb") as db_local:
                db_local.write(db_response)
            local_md5 = get_local_md5checksum(local_db_file, log)

        log.debug("\t => successfully downloaded new {} file".format(db_name))
        log.debug(f"\t => MD5 of downloaded file: {local_md5}")
    except urllib.error.HTTPError:
        msg = f"Sorry, could not find file {remote_db_file}!\n\n" \
              f"Possibly, version {version} of {db_name.upper()} does not exist?"
        return False, msg
    except socket.timeout:
        msg = "Reference file took too long to download. :-( Maybe the connection is slow or you need a proxy?"
        return False, msg

    log.debug("\tCreating parsed files...")
    version = hla_embl_parser.make_parsed_files(use_dbname, ref_path_temp, log)

    success, msg = make_blast_db(use_dbname, ref_path_temp, blast_path, log)

    if success:
        curr_md5_file = os.path.join(ref_path_temp, f"curr_md5_{use_dbname}.txt")
        log.debug(f"Writing md5 checksum to local file {curr_md5_file}...")
        with open(curr_md5_file, "w") as g:
            g.write(f"{local_md5} {datetime.datetime.now().strftime('%d.%m.%y')}")

        update_msg = f"Updated the reference data for {db_name.upper()} to version {version}."
        move_files(ref_path_temp, reference_local_path, db_name, log)
    else:
        log.error(msg)
        update_msg = f"Tried to update the reference data for {db_name.upper()} to version {version}, "
        update_msg += f"but could not process it (see below). We'll continue to use the old files for now.\n\n{msg}"

    log.info(update_msg)
    return success, update_msg


def make_restricted_db(db_name, ref_path, restricted_to, target_dir, blast_path, log):
    """creates a limited version of the given database, restricted to the given alleles
    """
    log.info(f"Create local reference version of {db_name} restricted to {', '.join(restricted_to)}...")
    if db_name == "kir":
        use_dbname = "KIR"  # biological databases and consistency in naming are arch enemies
    else:
        use_dbname = db_name

    os.makedirs(target_dir, exist_ok=True)

    log.debug("\tCreating parsed files...")
    hla_embl_parser.make_parsed_files(use_dbname, ref_path, log,
                                      restricted_to=restricted_to,
                                      target_dir=target_dir)

    success, msg = make_blast_db(use_dbname, target_dir, blast_path, log)

    if success:
        log.info("Success!")
    else:
        log.error(msg)
    return success, msg


def update_country_data(target_dir: str, proxy: str, log: logging.Logger) -> Tuple[bool, str | None]:
    """Get current list of possible countries from GitHub and store them locally.

    Returns tuple:
        - success (bool)
        - error message (if any encountered, else None)
    """
    log.info("Retrieving current country list...")
    try:
        data = read_remote_file(COUNTRY_URL, proxy=proxy, timeout=20, log=log)
    except Exception as E:
        return False, f"Could not retrieve list of possible countries from {COUNTRY_URL}:\n\n{E}"

    log.info("Updating local file...")
    country_file = Path(target_dir) / COUNTRY_FILE

    with open(country_file, "w") as g:
        g.write(data)
    log.info("Success!")
    return True, None


# ========================================================

def start_log(include_lines=False, error_to_email=False, info_to_file=False,
              debug_to_file=False,
              elaborate=False, level="DEBUG"):
    """starts a logger and returns it for logging,
    if log file is wanted, info_to_file should be the log file destination (path/filename);
    default logging level is set to DEBUG
    """
    import socket

    script_name = os.path.basename(__file__)
    log = logging.getLogger(script_name)

    # set level:
    level = level.upper()
    if level == "DEBUG":
        log.setLevel(logging.DEBUG)
    elif level == "INFO":
        log.setLevel(logging.INFO)
    elif level == "WARNING":
        log.setLevel(logging.WARNING)
    elif level == "ERROR":
        log.setLevel(logging.ERROR)
    elif level == "CRITICAL":
        log.setLevel(logging.CRITICAL)
    else:
        log.setLevel(logging.DEBUG)
        log.warning("Unknown loglevel '%s'! Using DEBUG instead..." % level)

    # define logging format:
    if include_lines:
        formatter = logging.Formatter('%(levelname)s [%(asctime)s] - [ln.%(lineno)d] %(message)s')
    elif elaborate:
        formatter = logging.Formatter(
            "[%(asctime)s] %(levelname)s [%(filename)s:%(lineno)s - %(funcName)s] %(message)s")
    else:
        formatter = logging.Formatter('%(levelname)s [%(asctime)s] - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    # establish stream handler (instead of print):
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    log.addHandler(stream_handler)

    # establish email handler, if needed: (level ERROR)
    if error_to_email:
        myport = 25
        myserver = '192.168.2.23'
        nodename = socket.gethostname()
        email_adresses = ['schoene@dkms-lab.de']
        email_handler = logging.handlers.SMTPHandler(mailhost=(myserver, myport),
                                                     fromaddr='%s@%s' % (script_name, nodename),
                                                     toaddrs=email_adresses,
                                                     subject='Error in %s' % script_name)
        email_handler.setFormatter(formatter)
        email_handler.setLevel(logging.ERROR)
        log.addHandler(email_handler)

    # establish file handler, if needed: (level INFO)
    if info_to_file:
        file_handler = logging.FileHandler(info_to_file)
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.INFO)
        log.addHandler(file_handler)

    if debug_to_file:
        log.debug_handler = logging.FileHandler(debug_to_file)
        log.debug_handler.setFormatter(formatter)
        log.debug_handler.setLevel(logging.DEBUG)
        log.addHandler(log.debug_handler)

    return log


pass

# ===========================================================
# main:
# def main():
#     log = start_log(level="DEBUG")
#     log.info("<Start>")
#     db_list = ["hla", "kir"]
#     blast_path = r"Y:\Projects\typeloader\blast-2.7.1+\bin"
#     reference_local_path = r"Y:\Projects\typeloader\temp\_general\reference_data"
#
#     for db_name in db_list:
#         new_version, msg = check_database(db_name, reference_local_path, log, skip_if_updated_today=False)
#         if new_version:
#             success, update_msg = update_database(db_name, reference_local_path, blast_path, log)
#
#     log.info("<End>")
#
#
if __name__ == '__main__':
    log = start_log(level="DEBUG")
    log.info("<Start>")
    # proxy = "10.78.205.144:3128"

    db_name = "hla"
    reference_local_path = r"\\labor.local\data\Typeloader\_general\reference_data"
    blast_path = r"C:\Users\schoene\WorkFolders\Code\typeloader2\typeloader2\blastn"
    proxy=None

    update_database(db_name, reference_local_path, blast_path, proxy, log, version=None)
