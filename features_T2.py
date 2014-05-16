# -*- coding: utf-8 -*-
"""
Created on Thu May 15 16:15:44 2014

@ author (C) Cristina Gallego, University of Toronto
----------------------------------------------------------------------
 """

import os, os.path
import sys
import string

from sys import argv, stderr, exit
import shlex, subprocess
import re

import numpy as np
import dicom
import psycopg2
import sqlalchemy as al
import sqlalchemy.orm

import pandas as pd
from query_database import *

from inputs_init import *
from display import *
import pylab    



if __name__ == '__main__':
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__))
    print path_rootFolder
   
    # Open filename list
    file_ids = open(sys.argv[1],"r")
    init_flag=1
    
    for fileline in file_ids:
        # Get the line: StudyNumber    DicomExamNumber    MRN    chosen_lesions_id    StudyDate    SeriesID    image_pos_pat    image_ori_pat
        fileline = fileline.split()
        cond = fileline[0] 
        StudyID = fileline[1]  
        DicomExamNumber = fileline[2]
        Lesions_id = fileline[3]
        dateID = fileline[4]
        SeriesID = fileline[5] # corresponds to dynamic sequence;
        

