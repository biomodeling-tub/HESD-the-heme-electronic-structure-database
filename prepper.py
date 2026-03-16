import os
import re
import shutil
import filecmp
import json
import random
import itertools
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import pandas as pd
from datetime import datetime
import MDAnalysis as mda
from scipy.stats import gaussian_kde
import glob
from joblib import Parallel, delayed


#Global Flags
debug=False
verbose=False


class prepper:
    def __init__(self):
        self.source_dir="/media/mamnas2/jones/beb00042/"
        self.logfile_dir="/home/pbuser/Desktop/PhD_WORK/heme/logfiles/"
        self.copy_check=True

    def prep(self):
        self.search_and_copy_logs()
        pairs = self.find_file_pairs()
        for original, duplicate in pairs:
            if os.path.exists(original) and os.path.exists(duplicate):
                prepper.compare_and_cleanup_duplicates(original, duplicate)
        self.filter_out_incomplete_logs()
        return

    def search_and_copy_logs(self):
        """
        Recursively search source_dir for files that match the pattern:
        - Exactly 4 chars, followed by 2 digits, then up to 2 letters, ending in .log
        Copy them to logfile_dir *only if* the file doesn't already exist
        in logfile_dir.
        """
        os.makedirs(self.logfile_dir, exist_ok=True)
        pattern = re.compile(r'^.{4}[0-9]{2}[A-Za-z]{0,2}\.log$')
        for root, dirs, files in os.walk(self.source_dir):
            for filename in files:
                if pattern.match(filename):
                    src = os.path.join(root, filename)
                    dest = os.path.join(self.logfile_dir, filename)
                    if self.copy_check==True:
                        if not os.path.exists(dest):
                            shutil.copy2(src, dest)
                            print(f"Copied: {src} -> {dest}")
                        else:
                            print(f"Skipped (already exists): {dest}")
                    else:
                        shutil.copy2(src, dest)
                        print(f"Copied: {src} -> {dest}")
        return

    def get_last_line(filepath):
        """Return the last line of a file or an empty string if file is empty."""
        try:
            with open(filepath, 'r') as f:
                lines = f.read().splitlines()
                return lines[-1].strip() if lines else ''
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            return ''

    def compare_and_cleanup_duplicates(original, duplicate):
        """
        Compare the content of two files:
        1. If they are identical, remove the duplicate.
        2. Otherwise, check the last line of each file:
            - If the last line is not "Normal Termination of Gaussian16", remove that file.
        """
        termination_string = "Normal Termination of Gaussian16"
        if filecmp.cmp(original, duplicate, shallow=False):
            os.remove(duplicate)
            print(f"[IDENTICAL] Removed duplicate: {duplicate}")
            return
        last_line_original = prepper.get_last_line(original)
        last_line_duplicate = prepper.get_last_line(duplicate)
        original_ok = (last_line_original == termination_string)
        duplicate_ok = (last_line_duplicate == termination_string)
        if not original_ok:
            os.remove(original)
            print(f"[NO TERMINATION] Removed original: {original}")
        if not duplicate_ok:
            os.remove(duplicate)
            print(f"[NO TERMINATION] Removed duplicate: {duplicate}")
        return

    def find_file_pairs(self):
        """
        Scan 'directory' for .log files, group them by base name,
        and return a list of (original, duplicate) tuples.

        An "original" is something like:  file.log
        A "duplicate" is something like: file(1).log, file(2).log, etc.
        """
        pattern = re.compile(r'^(.*?)(\(\d+\))?\.log$')
        base_dict = {}
        for filename in os.listdir(self.logfile_dir):
            if filename.endswith('.log'):
                match = pattern.match(filename)
                if match:
                    base_name, suffix = match.group(1), match.group(2)
                    full_path = os.path.join(self.logfile_dir, filename)
                    base_dict.setdefault(base_name, []).append((filename, full_path))
        pairs = []
        for base_name, files_info in base_dict.items():
            original = None
            duplicates = []
            for (fname, fpath) in files_info:
                if re.search(r'\(\d+\)', fname):
                    duplicates.append((fname, fpath))
                else:
                    original = (fname, fpath)
            for dup in duplicates:
                if original:
                    pairs.append((original[1], dup[1]))
        return pairs

    def filter_out_incomplete_logs(self):
        """
        Checks each .log file in 'self.logfile_dir' to see if its last line ends with 'Normal termination'.
        If not, moves the file to a 'discarded_logs' folder (created if needed).

        :param self.logfile_dir: Directory where the .log files are located.
        """
        #discarded_dir = os.path.join(self.logfile_dir, "discarded_logs")
        #os.makedirs(discarded_dir, exist_ok=True)
        for filename in os.listdir(self.logfile_dir):
            if filename.lower().endswith(".log"):
                logfile_path = os.path.join(self.logfile_dir, filename)
                if not self.ends_with_normal_termination(logfile_path):
                    #discarded_path = os.path.join(discarded_dir, filename)
                    #shutil.move(logfile_path, discarded_path) #store away
                    #print(f"Moved incomplete: {filename} -> {discarded_path}")
                    os.remove(logfile_path) #delete
                    print(f"Deleted incomplete: {filename}")
        return

    def ends_with_normal_termination(self, filepath):
        """
        Returns True if the last line of 'filepath' ends with 'Normal termination', False otherwise.
        """
        try:
            with open(filepath, 'rb') as f:
                f.seek(0, os.SEEK_END)
                buffer = b''
                size = f.tell()
                block_size = 1024
                while size > 0:
                    read_size = min(block_size, size)
                    f.seek(size - read_size)
                    data = f.read(read_size)
                    buffer = data + buffer
                    size -= read_size
                    if b'\n' in data:
                        break
                last_lines = buffer.decode('utf-8', errors='replace').splitlines()
                if not last_lines:
                    return False  # empty file
                last_line = last_lines[-1].strip()
                return last_line.startswith("Normal termination")
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            return False

