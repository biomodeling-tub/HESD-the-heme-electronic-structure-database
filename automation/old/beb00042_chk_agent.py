import os, sys, json
import subprocess as proc
from datetime import datetime
from datetime import date
import hashlib
from beb00042_agent import *

if __name__ == "__main__":
    chk_files_in_folder(pdb_id=pdb_id, chk_files_dict=chk_files_dict)