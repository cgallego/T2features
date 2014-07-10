# -*- coding: utf-8 -*-
"""
Procees T2 signal from batch list 

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
import processDicoms

from inputs_init import *
from display import *
from features_dynamic import *
from features_morphology import *
from features_texture import *
import pylab      

from features_T2 import *
import annot

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
        BenignNMaligNAnt = fileline[9]
        Diagnosis = fileline[10]
        										

        #############################
        ###### 2) Querying Research database for clinical, pathology, radiology data
        ############################# 
        print "Executing SQL connection..."
        # Format query StudyID
        if (len(StudyID) >= 4 ): fStudyID=StudyID
        if (len(StudyID) == 3 ): fStudyID='0'+StudyID
        if (len(StudyID) == 2 ): fStudyID='00'+StudyID
        if (len(StudyID) == 1 ): fStudyID='000'+StudyID
           
        # Format query redateID
        redateID = dateID[0:4]+'-'+dateID[4:6]+'-'+dateID[6:8]

        # perform query        
        try:
            queryData = Query()
            is_mass, colLabelsmass, is_nonmass, colLabelsnonmass = queryData.queryDatabase4T2(fStudyID, redateID)
            print is_mass
            print is_nonmass
            
        except Exception: 
            queryData = Query()
            queryData.queryDatabasewNoGuiNodate(fStudyID)
                    
        # ask for info about lesion row data from query
        if "mass" == cond[0:-1]:
            rowCase = raw_input('pick row for MASS (0-n) or x: ')
            if (rowCase == 'x'):
                ## append collection of cases 
                massframe = pd.DataFrame(data=array( ['NA'] * 10))
                massframe = massframe.transpose()
                massframe.columns = ['finding.side_int','finding.size_x_double','finding.size_y_double','finding.size_z_double','finding.mri_dce_init_enh_int','finding.mri_dce_delay_enh_int','finding.curve_int','finding.mri_nonmass_dist_int', 'finding.mri_nonmass_int_enh_int','finding.t2_signal_int']
            else:
                massframe = pd.DataFrame(data=array( is_mass[int(rowCase)] ))
                massframe = massframe.transpose()
                massframe.columns = list(colLabelsmass)
            
        if "nonmass" == cond[0:-1]:
            rowCase = raw_input('pick row for NONMASS (0-n) or x: ')
            if (rowCase == 'x'):
                ## append collection of cases 
                nonmassframe = pd.DataFrame(data=array( ['NA'] * 10))
                nonmassframe = nonmassframe.transpose()
                nonmassframe.columns = ['finding.side_int','finding.size_x_double','finding.size_y_double','finding.size_z_double','finding.mri_dce_init_enh_int','finding.mri_dce_delay_enh_int','finding.curve_int','finding.mri_nonmass_dist_int', 'finding.mri_nonmass_int_enh_int','finding.t2_signal_int']
            else:
                nonmassframe = pd.DataFrame(data=array( is_nonmass[int(rowCase)] ))
                nonmassframe = nonmassframe.transpose()
                nonmassframe.columns = list(colLabelsnonmass)
                    
        #slice data, get only 1 record    
        if "mass" == cond[0:-1]:
            dataCase = massframe 
        if "nonmass" == cond[0:-1]: 
            dataCase = nonmassframe
            
        casesFrame = pd.DataFrame(columns=dataCase.columns)
        casesFrame = casesFrame.append(dataCase) # 20             
               
        # Choose series T2
        # Get study image folder
        if "mass" in cond:
            img_folder = 'Z:/Cristina/MassNonmass/mass/'
        if "nonmass" in cond:  
            img_folder = 'Z:/Cristina/MassNonmass/nonmass/'
            
        [abspath_ExamID, eID, SeriesIDall, studyFolder, dicomInfo] = processDicoms.get_series(StudyID, img_folder)
                        
        # ask for which series to load
        print "\n----------------------------------------------------------"
        choseSerie = raw_input('Enter n T2 Series to load (0-n), or x if NO T2w sequence and pass: ')
        
        #############################
        # Reveal annotations                                      
        #############################
        annotflag = False 
        annotationsfound = [] 
        for iSer in SeriesIDall:
            exam_loc = img_folder+StudyID+'/'+eID+'/'+iSer
            print "Path Series annotation inspection: %s" % iSer
            annotationsfound, annotflag = annot.list_ann(exam_loc, annotflag, annotationsfound) 
            
        if (choseSerie == 'x'):
            ## append collection of cases 
            casesFrame['T2SeriesID']='NONE'
            casesFrame['cond']=cond
            casesFrame['id']=fStudyID
            casesFrame['DicomExamNumber']=DicomExamNumber
            casesFrame['Lesions_id']=Lesions_id
            casesFrame['BenignNMaligNAnt']=BenignNMaligNAnt
            casesFrame['Diagnosis']=Diagnosis
            casesFrame['T2dims']="[]"
            casesFrame['T2spacing']="[]"
            casesFrame['T2fatsat']="[]"
            # append to frame        
            casesFrame['T2_muscleSI']="[]"
            casesFrame['T2_muscleSIstd']="[]"
            casesFrame['muscle_scalar_range']="[]"
            casesFrame['bounds_muscleSI']="[]"
            casesFrame['T2_lesionSI']="[]"
            casesFrame['T2_lesionSIstd']="[]"
            casesFrame['lesion_scalar_range']="[]"
            casesFrame['LMSIR']="[]"
            if (annotflag):
                casesFrame['annotations'] = str( annotationsfound )
            else:
                casesFrame['annotations'] = "NA"
            casesFrame['T2min_F_r_i']="[]"
            casesFrame['T2max_F_r_i']="[]"
            casesFrame['T2mean_F_r_i']="[]"
            casesFrame['T2var_F_r_i']="[]"
            casesFrame['T2skew_F_r_i']="[]"
            casesFrame['T2kurt_F_r_i']="[]"
            casesFrame['T2grad_margin']="[]"
            casesFrame['T2grad_margin_var']="[]"
            casesFrame['T2RGH_mean']="[]"
            casesFrame['T2RGH_var']="[]"
            casesFrame.set_index('id',inplace=False)     
            print casesFrame
            
        else:    
            T2SeriesID = SeriesIDall[int(choseSerie)]
            
            ## append collection of cases 
            casesFrame['T2SeriesID']=T2SeriesID
            casesFrame['cond']=cond
            casesFrame['id']=fStudyID
            casesFrame['DicomExamNumber']=DicomExamNumber
            casesFrame['Lesions_id']=Lesions_id
            casesFrame['BenignNMaligNAnt']=BenignNMaligNAnt
            casesFrame['Diagnosis']=Diagnosis
                    
            path_T2Series = img_folder+StudyID+'/'+eID+'/'+T2SeriesID
            print "\nPath to T2 location: %s" %  path_T2Series
           
            #############################                  
            ###### 3) Process T2 and visualize
            #############################
            ###### Start by Loading 
            print "Start by loading T2 volume..."
            load = Inputs_init()        
            load.readT2(path_T2Series)
            
            casesFrame['T2dims']= str(load.T2dims)
            casesFrame['T2spacing']=str(load.T2spacing)
            casesFrame['T2fatsat']=str(load.T2fatsat)
            print casesFrame
            
            print "\n Load Segmentation..."
            data_loc='Z:\Cristina\MassNonmass'+os.sep+cond[:-1]
            print "\n Load DCE-MRI series..."
            [series_path, phases_series, lesionID_path] = load.readVolumes(data_loc, StudyID, DicomExamNumber, SeriesID, Lesions_id)
            print "Path to lesion segmentation: %s" % lesionID_path
            lesion3D = load.loadSegmentation(lesionID_path)
            print "Data Structure: %s" % lesion3D.GetClassName()
            print "Number of points: %d" % int(lesion3D.GetNumberOfPoints())
            print "Number of cells: %d" % int(lesion3D.GetNumberOfCells())
            
            print "\n Visualize volumes..."
            # Create only 1 display
            loadDisplay = Display()
            lesion3D_mesh = loadDisplay.addSegment(lesion3D, (0,1,1), interact=False)
            loadDisplay.visualize(load.DICOMImages, load.image_pos_pat, load.image_ori_pat, sub=True, postS=4, interact=False)
            print "\n Visualize segmentation ..."
            loadDisplay.addT2visualize(load.T2Images, load.T2image_pos_pat, load.T2image_ori_pat, load.T2dims, load.T2spacing, interact=True)
            
            #############################
            # Extract Lesion and Muscle Major pectoralies signal                                   
            #############################             
            T2 = features_T2()  
            [T2_muscleSI, muscle_scalar_range]  = T2.extract_muscleSI(load.T2Images, load.T2image_pos_pat, load.T2image_ori_pat, loadDisplay.origin, loadDisplay.iren1, loadDisplay.renderer1, loadDisplay.picker, loadDisplay.xImagePlaneWidget, loadDisplay.yImagePlaneWidget, loadDisplay.zImagePlaneWidget)
            print "ave. T2_muscleSI: %d" % mean(T2_muscleSI)
                    
            [T2_lesionSI, lesion_scalar_range]  = T2.extract_lesionSI(load.T2Images, lesion3D, load.T2image_pos_pat, load.T2image_ori_pat)
            print "ave. T2_lesionSI: %d" % mean(T2_lesionSI)
            
            LMSIR = mean(T2_lesionSI)/mean(T2_muscleSI)
            print "LMSIR: %d" % LMSIR
            
            # append to frame        
            casesFrame['T2_muscleSI']=str(mean(T2_muscleSI))
            casesFrame['T2_muscleSIstd']=str(std(T2_muscleSI))
            casesFrame['muscle_scalar_range']=str(muscle_scalar_range)
            casesFrame['bounds_muscleSI']=str(T2.bounds_muscleSI)
            
            casesFrame['T2_lesionSI']=str(mean(T2_lesionSI))
            casesFrame['T2_lesionSIstd']=str(std(T2_lesionSI))
            casesFrame['lesion_scalar_range']=str(lesion_scalar_range)
                    
            casesFrame['LMSIR']= str( LMSIR ) 
            print "\n Is there anotations..."
            print annotflag
            if (annotflag):
                casesFrame['annotations'] = str( annotationsfound )
            else:
                casesFrame['annotations'] = "NA"
                
            casesFrame.set_index('id',inplace=False)       
            print casesFrame
        
            #############################
            # Extract morphological and margin features from T2                                   
            #############################
            print "\n Extract T2 Morphology features..."
            morphoT2features = T2.extractT2morphology(load.T2Images, lesion3D, load.T2image_pos_pat, load.T2image_ori_pat)
            print "\n=========================================="
            print morphoT2features
            print "\n=========================================="
            
            print "\n Extract T2 Texture features..."
            textureT2features = T2.extractT2texture(load.T2Images, lesion3D, load.T2image_pos_pat, load.T2image_ori_pat)
            print "\n=========================================="
            print textureT2features
            print "\n=========================================="
            
            #############################
            ###### Finish tidying up and save to file
            ## append collection of cases
            #############################            
            morphoT2features['id']=fStudyID
            morphoT2features.set_index('id',inplace=False)
            casesFrame = pd.merge(casesFrame, morphoT2features, on='id', how='inner')
            
            # deal with closing windows, plots, renders, actors
            pylab.close('all')
            loadDisplay.iren1.TerminateApp()
            loadDisplay.renWin1.Finalize() 
        
        os.chdir(path_rootFolder)
        if(init_flag): 
            allcasesFrame = pd.DataFrame(columns=casesFrame.columns)
            init_flag=False  
            
        allcasesFrame = allcasesFrame.append(casesFrame, ignore_index=True)
        allcasesFrame.to_csv('casesFrames_T2features1.csv')  
         
            
    file_ids.close()