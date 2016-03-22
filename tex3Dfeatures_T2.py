# -*- coding: utf-8 -*-
"""
Extract pool of descriptors for T2 features

USAGE: 
=============
from features_T2 import *
T2 = features_T2()  
[T2_muscleSI, muscle_scalar_range]  = T2.extract_muscleSI(load.T2Images, load.T2image_pos_pat, load.T2image_ori_pat, loadDisplay.origin, loadDisplay.iren1, loadDisplay.renderer1, loadDisplay.picker, loadDisplay.xImagePlaneWidget, loadDisplay.yImagePlaneWidget, loadDisplay.zImagePlaneWidget)
[T2_lesionSI, lesion_scalar_range]  = T2.extract_lesionSI(load.T2Images, lesion3D, load.T2image_pos_pat, load.T2image_ori_pat)

Created on Wed Jun 04 13:50:08 2014

@ author (C) Cristina Gallego, University of Toronto, 2013
----------------------------------------------------------------------
"""
import os, os.path
import sys
import string
from sys import argv, stderr, exit
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import processDicoms
import dicom

import pandas as pd
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

from display import *
from features_morphology import *

from scipy import stats
from skimage.feature import greycomatrix, greycoprops

#!/usr/bin/env python
class tex3Dfeatures_T2(object):
    """
    USAGE:
    =============
    T2 = tex3Dfeatures_T2()
    Based on features_T2, this version inplements 3D texture and features based on nondirectional GLCM,
    as proposed by Chen et al (http://www.ncbi.nlm.nih.gov/pubmed/17763361) MRM 2007
    """
    def __init__(self): 
        """ initialize visualization with standard vtk actors, renders, windowsn, interactors """           
        # use cell picker for interacting with the image orthogonal views.
        self.morphologyT2Features = []    
        self.textureT2Features = []
        self.T2series = " "
        
        # Initialize features
        self.i_var = []
        self.alln_F_r_i=[]
        self.allmin_F_r_i=[]
        self.allmax_F_r_i=[]
        self.allmean_F_r_i=[]
        self.allvar_F_r_i=[]
        self.allskew_F_r_i=[]
        self.allkurt_F_r_i=[]
        self.grad_margin = 0
        self.grad_margin_var = 0   
        self.RGH_mean = 0
        self.RGH_var = 0
        
    def __call__(self):       
        """ Turn Class into a callable object """
        features_T2()
    
    def SelectPolygons(self, widget, event_string):
        # As can be seen the callback takes two arguments.  The first being the object that generates the event 
        # and the second argument the event name (which is a string).        
        self.planes = vtk.vtkPlanes()
        self.boxWidget.GetPlanes(self.planes)
        self.boundsPlane_presel = self.planes.GetPoints().GetBounds()
        
        return

    # Create a Python function to create the text for the text mapper used
    # to display the results of picking.
    def annotatePick(self, object, event):
        if(self.picker.GetCellId() < 0):
            self.textActor.VisibilityOff()     
        else:
    		print "pick"
    		selPt = self.picker.GetSelectionPoint()
    		self.pickPos = self.picker.GetPickPosition()
    		print self.pickPos
    				
    		self.textMapper.SetInput("(%.6f, %.6f, %.6f)"% self.pickPos)
    		self.textActor.SetPosition(selPt[:2])
    		self.textActor.VisibilityOn()
      
        return         
    
    def histeq(self, im, nbr_bins=256):
        #get image histogram  
        imhist,bins = histogram(im.flatten(),nbr_bins,normed=True)
        cdf = imhist.cumsum() #cumulative distribution function
        cdf = 255 * cdf / cdf[-1] #normalize
        
        #use linear interpolation of cdf to find new pixel values
        im2 = interp(im.flatten(),bins[:-1],cdf)
        
        return im2.reshape(im.shape), cdf 
        
    def extract_muscleSI(self, T2image, image_pos_pat, image_ori_pat, iren, ren, picker, xplane, yplane, zplane):
        """ extract_muscleSI: Place manually a widget over muscle, extract reference SI
        
        INPUTS:
        =======        
        images: (vtkImageData)   list of Input image to Transform
        OUTPUTS:
        =======
        muscleSI (float)    Signal intensity from muscle
        muscleSIcoords (float[3])            cords where Signal intensity from muscle is measured
        
        """
        ## Transform T2 img
        loadDisplay = Display()
        
        # Proceed to build reference frame for display objects based on DICOM coords   
        [t_T2image, transform_cube] = loadDisplay.dicomTransform(T2image, image_pos_pat, image_ori_pat)
            
        # Calculate the center of the volume
        t_T2image.UpdateInformation() 
               
        print "\nBoxwidget placed..."
        #################################################################
        # The box widget observes the events invoked by the render window
        # interactor.  These events come from user interaction in the render
        # window.
        # Place the interactor initially. The output of the reader is used to
        # place the box widget.
        self.boxWidget = vtk.vtkBoxWidget()
        self.boxWidget.SetInteractor(iren)
        self.boxWidget.SetPlaceFactor(1)
        self.boxWidget.SetInput(t_T2image)
        self.boxWidget.ScalingEnabledOff()
        self.boxWidget.OutlineCursorWiresOff()
                
        # custom interaction
        self.picker = picker
        self.textMapper = vtk.vtkTextMapper()
        tprop = self.textMapper.GetTextProperty()
        tprop.SetFontFamilyToArial()
        tprop.SetFontSize(10)
        tprop.BoldOn()
        tprop.ShadowOn()
        tprop.SetColor(1, 0, 0)
                
        self.textActor = vtk.vtkActor2D()
        self.textActor.VisibilityOff() 
        self.textActor.SetMapper(self.textMapper)
        ren.AddActor2D(self.textActor)
        picker.AddObserver("EndPickEvent", self.annotatePick)
        #iren.Start()        
                
        # Construct a bounding box
        bwidg = [0,0,0,0,0,0]     
        bwidg[0] = self.pickPos[0]; bwidg[1] = self.pickPos[0]+5; 
        bwidg[2] = self.pickPos[1]; bwidg[3] = self.pickPos[1]+5;
        bwidg[4] = self.pickPos[2]; bwidg[5] = self.pickPos[2]+5;
        self.bounds_muscleSI = bwidg
        print "\nbounds_muscleSI "
        print self.bounds_muscleSI
        
        self.boxWidget.PlaceWidget( self.bounds_muscleSI )
        # turn off planes
        xplane.Off()
        yplane.Off()
        self.boxWidget.On()
        
        ##########
        ### Set image stencil for muscle
        # create a simple box VOI mask shape using previously found boundsPlane_preselected
        VOIStencil = vtk.vtkROIStencilSource()
        VOIStencil.SetShapeToBox()
        VOIStencil.SetInformationInput(t_T2image)
        VOIStencil.SetBounds( self.bounds_muscleSI )    
        VOIStencil.Update()
                                
        # cut the corresponding VOI region and set the background:
        extractVOI_imgstenc = vtk.vtkImageStencil()
        extractVOI_imgstenc.SetInput(t_T2image)
        extractVOI_imgstenc.SetStencil(VOIStencil.GetOutput())
        extractVOI_imgstenc.ReverseStencilOff()
        extractVOI_imgstenc.SetBackgroundValue(5000)
        extractVOI_imgstenc.SetInformation(t_T2image.GetInformation())
        extractVOI_imgstenc.Update()
            
        # take out average image
        finalmuscleSIIm = vtk.vtkImageData()
        finalmuscleSIIm = extractVOI_imgstenc.GetOutput()
                
        ## Display histogram 
        dims = finalmuscleSIIm .GetDimensions()
        scalars = finalmuscleSIIm.GetPointData().GetScalars()
        np_scalars = vtk_to_numpy(scalars)      
        np_scalars = np_scalars.reshape(dims[2], dims[1], dims[0]) 
        np_scalars = np_scalars.transpose(2,1,0)
        muscleSI = np_scalars[np_scalars<5000]
                
        muscle_scalar_range = [muscleSI.min(), muscleSI.max()]
        print "\nMuscle scalar Range:"
        print muscle_scalar_range[0], muscle_scalar_range[1]
        
        return muscleSI, muscle_scalar_range, self.bounds_muscleSI


    def extract_muscleS_mhaI(self, T2image, image_pos_pat, image_ori_pat, iren, ren, picker, xplane, yplane, zplane):
        """ extract_muscleSI: Place manually a widget over muscle, extract reference SI
        
        INPUTS:
        =======        
        images: (vtkImageData)   list of Input image to Transform
        OUTPUTS:
        =======
        muscleSI (float)    Signal intensity from muscle
        muscleSIcoords (float[3])            cords where Signal intensity from muscle is measured
        
        """
        ## Transform T2 img
        loadDisplay = Display()
        
        # Proceed to build reference frame for display objects based on DICOM coords   
        t_T2image = loadDisplay.mhaTransform(T2image, image_pos_pat, image_ori_pat)
            
        # Calculate the center of the volume
        t_T2image.UpdateInformation() 
               
        print "\nBoxwidget placed..."
        #################################################################
        # The box widget observes the events invoked by the render window
        # interactor.  These events come from user interaction in the render
        # window.
        # Place the interactor initially. The output of the reader is used to
        # place the box widget.
        self.boxWidget = vtk.vtkBoxWidget()
        self.boxWidget.SetInteractor(iren)
        self.boxWidget.SetPlaceFactor(1)
        self.boxWidget.SetInput(t_T2image)
        self.boxWidget.ScalingEnabledOff()
        self.boxWidget.OutlineCursorWiresOff()
                
        # custom interaction
        self.picker = picker
        self.textMapper = vtk.vtkTextMapper()
        tprop = self.textMapper.GetTextProperty()
        tprop.SetFontFamilyToArial()
        tprop.SetFontSize(10)
        tprop.BoldOn()
        tprop.ShadowOn()
        tprop.SetColor(1, 0, 0)
                
        self.textActor = vtk.vtkActor2D()
        self.textActor.VisibilityOff() 
        self.textActor.SetMapper(self.textMapper)
        ren.AddActor2D(self.textActor)
        picker.AddObserver("EndPickEvent", self.annotatePick)
        iren.Start()        
                
        # Construct a bounding box
        bwidg = [0,0,0,0,0,0]     
        bwidg[0] = self.pickPos[0]; bwidg[1] = self.pickPos[0]+5; 
        bwidg[2] = self.pickPos[1]; bwidg[3] = self.pickPos[1]+5;
        bwidg[4] = self.pickPos[2]; bwidg[5] = self.pickPos[2]+5;
        self.bounds_muscleSI = bwidg
        print "\nbounds_muscleSI "
        print self.bounds_muscleSI
        
        self.boxWidget.PlaceWidget( self.bounds_muscleSI )
        # turn off planes
        xplane.Off()
        yplane.Off()
        self.boxWidget.On()
        
        ##########
        ### Set image stencil for muscle
        # create a simple box VOI mask shape using previously found boundsPlane_preselected
        VOIStencil = vtk.vtkROIStencilSource()
        VOIStencil.SetShapeToBox()
        VOIStencil.SetInformationInput(t_T2image)
        VOIStencil.SetBounds( self.bounds_muscleSI )    
        VOIStencil.Update()
                                
        # cut the corresponding VOI region and set the background:
        extractVOI_imgstenc = vtk.vtkImageStencil()
        extractVOI_imgstenc.SetInput(t_T2image)
        extractVOI_imgstenc.SetStencil(VOIStencil.GetOutput())
        extractVOI_imgstenc.ReverseStencilOff()
        extractVOI_imgstenc.SetBackgroundValue(5000)
        extractVOI_imgstenc.SetInformation(t_T2image.GetInformation())
        extractVOI_imgstenc.Update()
            
        # take out average image
        finalmuscleSIIm = vtk.vtkImageData()
        finalmuscleSIIm = extractVOI_imgstenc.GetOutput()
                
        ## Display histogram 
        dims = finalmuscleSIIm .GetDimensions()
        scalars = finalmuscleSIIm.GetPointData().GetScalars()
        np_scalars = vtk_to_numpy(scalars)      
        np_scalars = np_scalars.reshape(dims[2], dims[1], dims[0]) 
        np_scalars = np_scalars.transpose(2,1,0)
        muscleSI = np_scalars[np_scalars<5000]
                
        muscle_scalar_range = [muscleSI.min(), muscleSI.max()]
        print "\nMuscle scalar Range:"
        print muscle_scalar_range[0], muscle_scalar_range[1]
        
        return muscleSI, muscle_scalar_range, self.bounds_muscleSI
        
        
    def load_muscleSI(self, t_T2image, m_bounds, iren, ren, picker, xplane, yplane, zplane):
        """ load_muscleSI: Place automatically a widget over muscle location from file
        
        INPUTS:
        =======        
        images: (vtkImageData)   list of Input image to Transform
        OUTPUTS:
        =======
        muscleSI (float)    Signal intensity from muscle
        muscleSIcoords (float[3])            cords where Signal intensity from muscle is measured
        
        """
        ## Transform T2 img
#        loadDisplay = Display()
#        
#        # Proceed to build reference frame for display objects based on DICOM coords   
#        [t_T2image, transform_cube] = loadDisplay.dicomTransform(T2image, image_pos_pat, image_ori_pat)
#            
#        # Calculate the center of the volume
#        t_T2image.UpdateInformation() 
               
        print "\nBoxwidget placed..."
        #################################################################
        # The box widget observes the events invoked by the render window
        # interactor.  These events come from user interaction in the render
        # window.
        # Place the interactor initially. The output of the reader is used to
        # place the box widget.
        self.boxWidget = vtk.vtkBoxWidget()
        self.boxWidget.SetInteractor(iren)
        self.boxWidget.SetPlaceFactor(1)
        self.boxWidget.SetInput(t_T2image)
        self.boxWidget.ScalingEnabledOff()
        self.boxWidget.OutlineCursorWiresOn()                
                
        # Construct a bounding box
        bwidg = [0,0,0,0,0,0]     
        bwidg[0] = m_bounds[0]; bwidg[1] = m_bounds[1]; 
        bwidg[2] = m_bounds[2]; bwidg[3] = m_bounds[3];
        bwidg[4] = m_bounds[4]; bwidg[5] = m_bounds[5];
        self.bounds_muscleSI = bwidg
        print "\nbounds_muscleSI "
        print self.bounds_muscleSI

        # add to visualize        
        self.boxWidget.PlaceWidget( self.bounds_muscleSI )
        self.boxWidget.On()
        #iren.Start() 
        
        ############################# Do extract_muscleSI 
        print "\n Re-extract muscle VOI? "
        rexorNot=0
        rexorNot = int(raw_input('type 1 to Re-extract or anykey: '))
        
        if rexorNot == 1:
            # custom interaction
            self.picker = picker
            self.textMapper = vtk.vtkTextMapper()
            tprop = self.textMapper.GetTextProperty()
            tprop.SetFontFamilyToArial()
            tprop.SetFontSize(10)
            tprop.BoldOn()
            tprop.ShadowOn()
            tprop.SetColor(1, 0, 0)
                    
            self.textActor = vtk.vtkActor2D()
            self.textActor.VisibilityOff() 
            self.textActor.SetMapper(self.textMapper)
            ren.AddActor2D(self.textActor)
            picker.AddObserver("EndPickEvent", self.annotatePick)
            iren.Start()        
                    
            # Construct a bounding box
            bwidg = [0,0,0,0,0,0]     
            bwidg[0] = self.pickPos[0]; bwidg[1] = self.pickPos[0]+5; 
            bwidg[2] = self.pickPos[1]; bwidg[3] = self.pickPos[1]+5;
            bwidg[4] = self.pickPos[2]; bwidg[5] = self.pickPos[2]+5;
            self.bounds_muscleSI = bwidg
            print "\nbounds_muscleSI "
            print self.bounds_muscleSI
            
            self.boxWidget.PlaceWidget( self.bounds_muscleSI )
            # turn off planes
            xplane.Off()
            yplane.Off()
            self.boxWidget.On()            
            
        ##########
        ### Set image stencil for muscle
        # create a simple box VOI mask shape using previously found boundsPlane_preselected
        VOIStencil = vtk.vtkROIStencilSource()
        VOIStencil.SetShapeToBox()
        VOIStencil.SetBounds( self.bounds_muscleSI )    
        VOIStencil.SetInformationInput(t_T2image)
        VOIStencil.Update()
                                
        # cut the corresponding VOI region and set the background:
        extractVOI_imgstenc = vtk.vtkImageStencil()
        extractVOI_imgstenc.SetInput(t_T2image)
        extractVOI_imgstenc.SetStencil(VOIStencil.GetOutput())
        extractVOI_imgstenc.ReverseStencilOff()
        extractVOI_imgstenc.SetBackgroundValue(5000)
        extractVOI_imgstenc.Update()
            
        # take out average image
        finalmuscleSIIm = vtk.vtkImageData()
        finalmuscleSIIm = extractVOI_imgstenc.GetOutput()
        finalmuscleSIIm.Update()
                
        ## Display histogram 
        dims = finalmuscleSIIm .GetDimensions()
        scalars = finalmuscleSIIm.GetPointData().GetScalars()
        np_scalars = vtk_to_numpy(scalars)      
        np_scalars = np_scalars.reshape(dims[2], dims[1], dims[0]) 
        np_scalars = np_scalars.transpose(2,1,0)
        muscleSI = np_scalars[np_scalars<5000]
        print "ave. T2_muscleSI: %d" % mean(muscleSI)
        
        muscle_scalar_range = [muscleSI.min(), muscleSI.max()]
        print "\nMuscle scalar Range:"
        print muscle_scalar_range[0], muscle_scalar_range[1]
        
        return muscleSI, muscle_scalar_range, self.bounds_muscleSI
        
        
    def load_muscleSI_mha(self, T2image, image_pos_pat, image_ori_pat, m_bounds, iren):
        """ load_muscleSI: Place automatically a widget over muscle location from file
        
        INPUTS:
        =======        
        images: (vtkImageData)   list of Input image to Transform
        OUTPUTS:
        =======
        muscleSI (float)    Signal intensity from muscle
        muscleSIcoords (float[3])            cords where Signal intensity from muscle is measured
        
        """
        ## Transform T2 img
        loadDisplay = Display()
        
        # Proceed to build reference frame for display objects based on DICOM coords   
        t_T2image = loadDisplay.mhaTransform(T2image, image_pos_pat, image_ori_pat)
            
        # Calculate the center of the volume
        t_T2image.UpdateInformation() 
               
        print "\nBoxwidget placed..."
        #################################################################
        # The box widget observes the events invoked by the render window
        # interactor.  These events come from user interaction in the render
        # window.
        # Place the interactor initially. The output of the reader is used to
        # place the box widget.
        self.boxWidget = vtk.vtkBoxWidget()
        self.boxWidget.SetInteractor(iren)
        self.boxWidget.SetPlaceFactor(1)
        self.boxWidget.SetInput(t_T2image)
        self.boxWidget.ScalingEnabledOff()
        self.boxWidget.OutlineCursorWiresOn()                
                
        # Construct a bounding box
        bwidg = [0,0,0,0,0,0]     
        bwidg[0] = m_bounds[0]; bwidg[1] = m_bounds[1]; 
        bwidg[2] = m_bounds[2]; bwidg[3] = m_bounds[3];
        bwidg[4] = m_bounds[4]; bwidg[5] = m_bounds[5];
        self.bounds_muscleSI = bwidg
        print "\nbounds_muscleSI "
        print self.bounds_muscleSI

        # add to visualize        
        self.boxWidget.PlaceWidget( self.bounds_muscleSI )
        self.boxWidget.On()
        
        ##########
        ### Set image stencil for muscle
        # create a simple box VOI mask shape using previously found boundsPlane_preselected
        VOIStencil = vtk.vtkROIStencilSource()
        VOIStencil.SetShapeToBox()
        VOIStencil.SetBounds( self.bounds_muscleSI )    
        VOIStencil.SetInformationInput(t_T2image)
        VOIStencil.Update()
                                
        # cut the corresponding VOI region and set the background:
        extractVOI_imgstenc = vtk.vtkImageStencil()
        extractVOI_imgstenc.SetInput(t_T2image)
        extractVOI_imgstenc.SetStencil(VOIStencil.GetOutput())
        extractVOI_imgstenc.ReverseStencilOff()
        extractVOI_imgstenc.SetBackgroundValue(5000)
        extractVOI_imgstenc.Update()
            
        # take out average image
        finalmuscleSIIm = vtk.vtkImageData()
        finalmuscleSIIm = extractVOI_imgstenc.GetOutput()
        finalmuscleSIIm.Update()
                
        ## Display histogram 
        dims = finalmuscleSIIm .GetDimensions()
        scalars = finalmuscleSIIm.GetPointData().GetScalars()
        np_scalars = vtk_to_numpy(scalars)      
        np_scalars = np_scalars.reshape(dims[2], dims[1], dims[0]) 
        np_scalars = np_scalars.transpose(2,1,0)
        muscleSI = np_scalars[np_scalars<5000]
                
        muscle_scalar_range = [muscleSI.min(), muscleSI.max()]
        print "\nMuscle scalar Range:"
        print muscle_scalar_range[0], muscle_scalar_range[1]
        
        return muscleSI, muscle_scalar_range, self.bounds_muscleSI
        
        
    def extract_lesionSI(self, t_T2image, lesion3D, loadDisplay):
        """ extract_lesionSI: Use lesion segmentation to extract lesion SI
        
        INPUTS:
        =======        
        images: (vtkImageData)   list of Input image to Transform
        lesion3D: (vtkPolyData)     segmentation
        OUTPUTS:
        =======
        lesionSI (float)    Signal intensity from lesion
        lesion_scalar_range (float[3])    SI range inside lesion
        
        """
        ## Transform T2 img
#        loadDisplay = Display()
#        
#        # Proceed to build reference frame for display objects based on DICOM coords   
#        [t_T2image, transform_cube] = loadDisplay.dicomTransform(T2image, image_pos_pat, image_ori_pat)
#            
#        # Update info
#        t_T2image.UpdateInformation() 
                
        # create a simple box VOI mask shape using previously found boundsPlane_preselected
        VOIStencil = vtk.vtkROIStencilSource()
        VOIStencil.SetShapeToBox()
        VOIStencil.SetBounds( lesion3D.GetBounds() )    
        VOIStencil.SetInformationInput(t_T2image)
        VOIStencil.Update()
        
        self.boxWidget.PlaceWidget( lesion3D.GetBounds() )
        self.boxWidget.On()
                                        
        # cut the corresponding VOI region and set the background:
        extractVOIlesion_imgstenc = vtk.vtkImageStencil()
        extractVOIlesion_imgstenc.SetInput(t_T2image)
        extractVOIlesion_imgstenc.SetStencil(VOIStencil.GetOutput())
        extractVOIlesion_imgstenc.ReverseStencilOff()
        extractVOIlesion_imgstenc.SetBackgroundValue(5000)
        extractVOIlesion_imgstenc.Update()
            
        # take out average image
        finallesionSIIm = vtk.vtkImageData()
        finallesionSIIm = extractVOIlesion_imgstenc.GetOutput()
                
        ## extract scalars 
        dims = finallesionSIIm.GetDimensions()
        scalars = finallesionSIIm.GetPointData().GetScalars()
        np_scalars = vtk_to_numpy(scalars)    
        np_scalars = np_scalars.reshape(dims[2], dims[1], dims[0]) 
        np_scalars = np_scalars.transpose(2,1,0)
        lesionSI = np_scalars[np_scalars<5000]
                
        print "\nLesion_scalar Range:"
        lesion_scalar_range = [lesionSI.min(), lesionSI.max()]
        print lesion_scalar_range[0], lesion_scalar_range[1]
        
#        fileT2name = pathSegment+'/'+nameSegment
#        T2imseg = vtk.vtkMetaImageWriter()
#        T2imseg.SetFileName(fileT2name+'_T2imgseg.mhd')
#        T2imseg.SetInput(extractVOIlesion_imgstenc.GetOutput())
#        T2imseg.Write()
        
        return lesionSI, lesion_scalar_range
        
    def extract_lesionSI_mha(self, T2image, lesion3D, image_pos_pat, image_ori_pat, loadDisplay):
        """ extract_lesionSI: Use lesion segmentation to extract lesion SI
        
        INPUTS:
        =======        
        images: (vtkImageData)   list of Input image to Transform
        lesion3D: (vtkPolyData)     segmentation
        OUTPUTS:
        =======
        lesionSI (float)    Signal intensity from lesion
        lesion_scalar_range (float[3])    SI range inside lesion
        
        """
        ## Transform T2 img
        loadDisplay = Display()
        
        # Proceed to build reference frame for display objects based on DICOM coords   
        t_T2image = loadDisplay.mhaTransform(T2image, image_pos_pat, image_ori_pat)
            
        # Update info
        t_T2image.UpdateInformation() 
        
        # create a simple box VOI mask shape using previously found boundsPlane_preselected
        VOIStencil = vtk.vtkROIStencilSource()
        VOIStencil.SetShapeToBox()
        VOIStencil.SetBounds( lesion3D.GetBounds() )    
        VOIStencil.SetInformationInput(t_T2image)
        VOIStencil.Update()
        
        self.boxWidget.PlaceWidget( lesion3D.GetBounds() )
        self.boxWidget.On()
                                        
        # cut the corresponding VOI region and set the background:
        extractVOIlesion_imgstenc = vtk.vtkImageStencil()
        extractVOIlesion_imgstenc.SetInput(t_T2image)
        extractVOIlesion_imgstenc.SetStencil(VOIStencil.GetOutput())
        extractVOIlesion_imgstenc.ReverseStencilOff()
        extractVOIlesion_imgstenc.SetBackgroundValue(5000)
        extractVOIlesion_imgstenc.Update()
            
        # take out average image
        finallesionSIIm = vtk.vtkImageData()
        finallesionSIIm = extractVOIlesion_imgstenc.GetOutput()
                
        ## extract scalars 
        dims = finallesionSIIm.GetDimensions()
        scalars = finallesionSIIm.GetPointData().GetScalars()
        np_scalars = vtk_to_numpy(scalars)    
        np_scalars = np_scalars.reshape(dims[2], dims[1], dims[0]) 
        np_scalars = np_scalars.transpose(2,1,0)
        lesionSI = np_scalars[np_scalars<5000]
                
        print "\nLesion_scalar Range:"
        lesion_scalar_range = [lesionSI.min(), lesionSI.max()]
        print lesion_scalar_range[0], lesion_scalar_range[1]
        
        return lesionSI, lesion_scalar_range
        
    
    def createMaskfromMesh(self, VOI_mesh, im):
        """ Takes an image and a VOI_mesh and returns a boolean image with only 1s inside the VOI_mesh """                
        
        
        VOIStencil = vtk.vtkROIStencilSource()
        VOIStencil.SetShapeToBox()
        VOIStencil.SetBounds( VOI_mesh.GetBounds() )    
        VOIStencil.SetInformationInput(im)
        VOIStencil.Update()
         
        # cut the corresponding white image and set the background:
        imgstenc = vtk.vtkImageStencil()
        imgstenc.SetInput(im)
        imgstenc.SetStencil(VOIStencil.GetOutput())
        imgstenc.ReverseStencilOff()
        imgstenc.SetBackgroundValue(5000)
        imgstenc.Update()
        
        # write to image        
        dims = im.GetDimensions()
        scalars = imgstenc.GetOutput().GetPointData().GetScalars()
        np_scalars = vtk_to_numpy(scalars)     
        np_scalars = np_scalars.reshape(dims[2], dims[1], dims[0]) 
        np_scalars = np_scalars.transpose(2,1,0)
                
        return np_scalars
        
        
    def extractT2morphology(self, t_T2image, VOI_mesh):
        """ Start pixVals for collection pixel values at VOI """
        pixVals_margin = []; pixVals = []
        Fmargin = {}; voxel_frameS = {}
        
        # necessary to read point coords
        VOIPnt = [0,0,0]
        ijk = [0,0,0]
        pco = [0,0,0]
        
        ## Transform T2 img
#        localloadDisplay = Display()
#        
#        # Proceed to build reference frame for display objects based on DICOM coords   
#        [t_T2image, transform_cube] = localloadDisplay.dicomTransform(T2image, image_pos_pat, image_ori_pat)
#        # Update info
#        t_T2image.UpdateInformation() 
        
        # create mask from segmenation
        np_VOI_mask = self.createMaskfromMesh(VOI_mesh, t_T2image)
        
        for j in range( VOI_mesh.GetNumberOfPoints() ):
            VOI_mesh.GetPoint(j, VOIPnt)      
            
            # extract pixID at location VOIPnt
            pixId = t_T2image.FindPoint(VOIPnt[0], VOIPnt[1], VOIPnt[2])
            im_pt = [0,0,0]
            
            t_T2image.GetPoint(pixId,im_pt)           
            inorout = t_T2image.ComputeStructuredCoordinates( im_pt, ijk, pco)
            if(inorout == 0):
                pass
            else:
                pixValx = t_T2image.GetScalarComponentAsFloat( ijk[0], ijk[1], ijk[2], 0)
                pixVals_margin.append(pixValx)
                
        # Now collect pixVals
        print "\n Saving %s" % 'Fmargin'
        Fmargin['Fmargin'] = pixVals_margin
        pixVals_margin = []
        
        # extract pixID at inside VOIPnt
        VOI_scalars = t_T2image.GetPointData().GetScalars()
        np_VOI_imagedata = vtk_to_numpy(VOI_scalars)     
        
        dims = t_T2image.GetDimensions()
        spacing = t_T2image.GetSpacing()
        np_VOI_imagedata = np_VOI_imagedata.reshape(dims[2], dims[1], dims[0]) 
        np_VOI_imagedata = np_VOI_imagedata.transpose(2,1,0)
        
        #################### HERE GET INTERNAL PIXELS IT AND MASK IT OUT   
        VOI_imagedata = np_VOI_mask[np_VOI_mask<5000]
        
        print "shape of VOI_imagedata  Clipped:"
        print VOI_imagedata.shape
            
        for j in range( len(VOI_imagedata) ):
            pixValx = VOI_imagedata[j]
            pixVals.append(pixValx)
                
        # Now collect pixVals
        print "\n Saving %s" % 'F'
        voxel_frameS['F'] = pixVals   
        pixVals = []
                
        ##############################################################
        # Collect to Compute inhomogeneity variance of uptake and other variables
        F_r_i =  array(voxel_frameS['F']).astype(float)
        n_F_r_i, min_max_F_r_i, mean_F_r_i, var_F_r_i, skew_F_r_i, kurt_F_r_i = stats.describe(F_r_i)
                
        print("Number of internal voxels: {0:d}".format(n_F_r_i))
        self.alln_F_r_i = n_F_r_i
        print("Minimum: {0:8.6f} Maximum: {1:8.6f}".format(min_max_F_r_i[0], min_max_F_r_i[1]))
        self.allmin_F_r_i = min_max_F_r_i[0]
        self.allmax_F_r_i = min_max_F_r_i[1]
        print("Mean: {0:8.6f}".format(mean_F_r_i))
        self.allmean_F_r_i = mean_F_r_i
        print("Variance F_r_i: {0:8.6f}".format(var_F_r_i))
        self.allvar_F_r_i = var_F_r_i
        print("Skew : {0:8.6f}".format(skew_F_r_i))
        self.allskew_F_r_i = (skew_F_r_i)
        print("Kurtosis: {0:8.6f}".format(kurt_F_r_i))
        self.allkurt_F_r_i = kurt_F_r_i
        
        # Extract features for sharpness of lesion margin, compute Margin gradient iii_var
        # The gradient is computed using convolution with a 3D sobel filter using scipy.ndimage.filters.sobel
        # The function generic_gradient_magnitude calculates a gradient magnitude using the function passed through derivative to calculate first derivatives. 
        F_rmargin =  array(Fmargin['Fmargin']).astype(float)

        # Collect to Compute variance of uptake and other variables
        sobel_grad_margin_delta = generic_gradient_magnitude(F_rmargin, sobel) 
    
        # compute feature Margin Gradient
        n, min_max, mean_sobel_grad_margin, var_sobel_grad_margin, skew, kurt = stats.describe(sobel_grad_margin_delta)
            
        print("Margin Gradient: {0:8.6f}".format(mean_sobel_grad_margin))
        self.grad_margin =  mean_sobel_grad_margin
        
        print("Variance of Margin Gradient: {0:8.6f}".format( var_sobel_grad_margin ))
        self.grad_margin_var = var_sobel_grad_margin
        
        ####################################
        # Radial gradient analysis ref[9] white paper
        ###################################
        # Radial gradient analysis is based on examination of the angles between voxel-value gradients
        # and lines intersecting a single point near the center of the suspect lesion, lines in radial directions. 
        # Radial gradient values are given by the dot product of the gradient direction and the radial direction.
        RGH_mean = []
        RGH_var = []        
        H_norm_p = []

        for j in range( VOI_mesh.GetNumberOfPoints() ):
            VOI_mesh.GetPoint(j, VOIPnt)
            
            r = array(VOIPnt)
            rc = array(VOI_mesh.GetCenter())
            norm_rdir = r-rc/linalg.norm(r-rc)
           
            # Find point for gradient vectors at the margin point
            pixId = t_T2image.FindPoint(VOIPnt[0], VOIPnt[1], VOIPnt[2])
            sub_pt = [0,0,0]            
            t_T2image.GetPoint(pixId, sub_pt)
            
            ijk = [0,0,0]
            pco = [0,0,0]
            
            grad_pt = [0,0,0];
            
            inorout = t_T2image.ComputeStructuredCoordinates( sub_pt, ijk, pco)
            if(inorout == 0):
                print "point outside data"
            else:
                t_T2image.GetPointGradient( ijk[0], ijk[1], ijk[2], t_T2image.GetPointData().GetScalars(), grad_pt)
                
            #############
            # Compute vector in the direction gradient at margin point
            grad_marginpt = array([grad_pt])
            norm_grad_marginpt = grad_marginpt/linalg.norm(grad_marginpt)
            
            # Compute dot product (unit vector for dot product)
            p_dot = dot(norm_grad_marginpt, norm_rdir)
            norm_p_dot = np.abs(p_dot)[0]
            
            H_norm_p.append(norm_p_dot)    
                
                
        # The histogram of radial gradient values quantifying the frequency of occurrence of the dot products in a given region of interest
        # radial gradient histogram. The hist() function now has a lot more options
        # first create a single histogram
                    
        # the histogram of the data with histtype='step'
#        plt.figure()
#        nsamples, bins, patches = plt.hist(array(H_norm_p), 50, normed=1, histtype='bar',facecolor='blue', alpha=0.75)
#        n, min_max, mean_bins, var_bins, skew, kurt = stats.describe(nsamples)
        
        mean_bins = np.mean(H_norm_p)
        var_bins = np.var(H_norm_p)
        
        print("\n mean RGB: {0:8.6f}".format( mean_bins ))
        print("variance RGB: {0:8.6f}".format( var_bins ))
        
        # Append data
        RGH_mean =  mean_bins 
        RGH_var = var_bins
        
        self.RGH_mean = RGH_mean
        self.RGH_var = RGH_var
        
        # add a line showing the expected distribution
        # create a histogram by providing the bin edges (unequally spaced)
        plt.xlabel('normalized dot product |R.G|')
        plt.ylabel('Probability')
        plt.title('radial gradient histogram')
        plt.grid(True)            
                
        ##################################################
        # orgamize into dataframe
        self.morphologyT2Features = DataFrame( data=array([[ self.allmin_F_r_i, self.allmax_F_r_i, self.allmean_F_r_i, self.allvar_F_r_i, self.allskew_F_r_i, self.allkurt_F_r_i,
            self.grad_margin, self.grad_margin_var, self.RGH_mean, self.RGH_var]]), 
            columns=['T2min_F_r_i', 'T2max_F_r_i', 'T2mean_F_r_i', 'T2var_F_r_i', 'T2skew_F_r_i', 'T2kurt_F_r_i', 
            'T2grad_margin', 'T2grad_margin_var', 'T2RGH_mean', 'T2RGH_var'])
        
        return self.morphologyT2Features
        
    
    def extractT2texture(self, t_T2image, image_pos_pat, image_ori_pat, VOI_mesh, patchdirname):
        ################################### 
        # Haralick et al. defined 10 grey-level co-occurrence matrix (GLCM) enhancement features 
        # 3D textrue features
#        self = funcD
#        T2image = self.load.T2Images
#        image_pos_pat = self.load.T2image_pos_pat
#        image_ori_pat = self.load.T2image_ori_pat
#        VOI_mesh = funcD.lesion3D
#        loadDisplay = funcD.loadDisplay
        ## Transform T2 img
        #localloadDisplay = Display()
        
        # Proceed to build reference frame for display objects based on DICOM coords   
        #[t_T2image, transform_cube] = loadDisplay.dicomTransform(T2image, image_pos_pat, image_ori_pat)
        # Update info
        t_T2image.UpdateInformation() 
        
        print "\nLoading VOI_mesh MASK... "
        VOIPnt = [0,0,0]; pixVals_margin = [];  pixVals = []
            
        #################### TO NUMPY
        VOI_scalars = t_T2image.GetPointData().GetScalars()
        np_VOI_imagedata = vtk_to_numpy(VOI_scalars)

        dims = t_T2image.GetDimensions()
        spacing = t_T2image.GetSpacing()
        np_VOI_imagedata = np_VOI_imagedata.reshape(dims[2], dims[1], dims[0])
        np_VOI_imagedata = array(np_VOI_imagedata.transpose(2,1,0))

        #****************************"
        print "Normalizing data..."
        np_VOI_imagedataflat = np_VOI_imagedata.flatten().astype(float)
        np_VOI_imagedata_num = (np_VOI_imagedataflat-min(np_VOI_imagedataflat))
        np_VOI_imagedata_den = (max(np_VOI_imagedataflat)-min(np_VOI_imagedataflat))
        np_VOI_imagedata_flatten = 255*(np_VOI_imagedata_num/np_VOI_imagedata_den)
        np_VOI_imagedata = np_VOI_imagedata_flatten.reshape(dims[0], dims[1], dims[2])
        
        # Prepare lesion localization and PATCH size for Texture analysis    
        # lesion_centroid        
        lesionBox = VOI_mesh.GetBounds()
        self.boxWidget.PlaceWidget( lesionBox )
        self.boxWidget.On()
        
        deltax = lesionBox[1] - lesionBox[0]
        print deltax
        deltay = lesionBox[3]-lesionBox[2]
        print  deltay
        deltaz = lesionBox[5]-lesionBox[4]
        print deltaz
        
        # find largest dimension and use that for lesion radius
        if deltax > deltay:
            lesion_dia = deltax/spacing[0] 
            print "deltax was largest %d... " % lesion_dia
        else: 
            lesion_dia = deltay/spacing[1]
            print "deltay was largest %d... " % lesion_dia
             
        lesion_radius = lesion_dia/2
        lesionthick = deltaz/spacing[2]/2
        print "VOI_efect_diameter %s... " % str(lesion_radius)
        
        lesion_location = VOI_mesh.GetCenter()
        print "lesion_location %s... " % str(lesion_location)
        lesion_patches = []
        
        ######### translate lesion location to ijk coordinates
        # sub_pre_transformed_image.FindPoint(lesion_location
        pixId = t_T2image.FindPoint(lesion_location[0], lesion_location[1], lesion_location[2])
        sub_pt = [0,0,0]            
        t_T2image.GetPoint(pixId, sub_pt)
        
        ijk = [0,0,0]
        pco = [0,0,0]
                
        inorout = t_T2image.ComputeStructuredCoordinates( sub_pt, ijk, pco)
        print "coresponding ijk_vtk indices:"
        print ijk
        ijk_vtk = ijk
        
        # Perform texture classification using grey level co-occurrence matrices (GLCMs).
        # A GLCM is a histogram of co-occurring greyscale values at a given offset over an image.
        # compute some GLCM properties each patch
        # p(i,j) is the (i,j)th entry in a normalized spatial gray-level dependence matrix;
        lesion_patches = np_VOI_imagedata[
                int(ijk_vtk[0] - lesion_radius-1):int(ijk_vtk[0] + lesion_radius+1),
                int(ijk_vtk[1] - lesion_radius-1):int(ijk_vtk[1] + lesion_radius),
                int(ijk_vtk[2] - lesionthick/2-1):int(ijk_vtk[2] + lesionthick/2+1) ]

        print '\n Lesion_patches:'
        print lesion_patches
        
        #################### RESAMPLE TO ISOTROPIC
        # pass to vtk
        lesion_patches_shape = lesion_patches.shape
        vtklesion_patches = lesion_patches.transpose(2,1,0)
        vtklesion_patches_data = numpy_to_vtk(num_array=vtklesion_patches.ravel(), deep=True, array_type=vtk.VTK_FLOAT)

        # probe into the unstructured grid using ImageData geometry
        vtklesion_patches = vtk.vtkImageData()
        vtklesion_patches.SetExtent(0,lesion_patches_shape[0]-1,0,lesion_patches_shape[1]-1,0,lesion_patches_shape[2]-1)
        vtklesion_patches.SetOrigin(0,0,0)
        vtklesion_patches.SetSpacing(spacing)
        vtklesion_patches.GetPointData().SetScalars(vtklesion_patches_data)

        # write ####################
        isowriter = vtk.vtkMetaImageWriter()
        isowriter.SetInput(vtklesion_patches)
        isowriter.SetFileName(patchdirname+"_T2w.mha")
        isowriter.Write()

        isopix = mean(vtklesion_patches.GetSpacing())
        resample = vtk.vtkImageResample ()
        resample.SetInput( vtklesion_patches )
        resample.SetAxisOutputSpacing( 0, isopix )
        resample.SetAxisOutputSpacing( 1, isopix )
        resample.SetAxisOutputSpacing( 2, isopix )
        resample.Update()
        isoImage = resample.GetOutput()

        #################### get isotropic patches
        ISO_scalars = isoImage.GetPointData().GetScalars()
        np_ISO_Image = vtk_to_numpy(ISO_scalars)

        isodims = isoImage.GetDimensions()
        isospacing = isoImage.GetSpacing()
        np_ISO_Imagedata = np_ISO_Image.reshape(isodims[2], isodims[1], isodims[0])
        np_ISO_Imagedata = array(np_ISO_Imagedata.transpose(2,1,0))

        #################### save isotropic patches
        fig = plt.figure()
        # patch histogram
        ax1 = fig.add_subplot(221)
        n, bins, patches = plt.hist(array(np_ISO_Imagedata.flatten()), 50, normed=1, facecolor='green', alpha=0.75)
        ax1.set_ylabel('histo')

        ax2 = fig.add_subplot(222)
        plt.imshow(np_ISO_Imagedata[:,:,np_ISO_Imagedata.shape[2]/2])
        plt.gray()
        ax2.set_ylabel('iso: '+str(isodims) )

        ax3 = fig.add_subplot(223)
        plt.imshow(lesion_patches[:,:,lesion_patches.shape[2]/2])
        plt.gray()
        ax3.set_ylabel('original: '+str(lesion_patches_shape))

        ax4 = fig.add_subplot(224)
        plt.imshow(np_VOI_imagedata[:,:,ijk_vtk[2]])
        plt.gray()
        ax4.set_ylabel('lesion centroid(ijk): '+str(ijk_vtk))

        # FInally display
        # plt.show()
        plt.savefig(patchdirname+'_T2w.png', format='png')
        ####################

        #################### 3D GLCM
        from graycomatrix3D import glcm3d

        patch = np_ISO_Imagedata.astype(np.uint8)
        lev = int(patch.max()+1) # levels
        
        # perfor glmc extraction in all 13 directions in 3D pixel neighbors
        g1 = glcm3d(lev, patch, offsets=[0,0,1])  # orientation 0 degrees (example same slices: equal to Nslices*0degree 2D case)
        g2 = glcm3d(lev, patch, offsets=[0,1,-1]) # orientation 45 degrees (example same slices: equal to Nslices*45degree 2D case)
        g3 = glcm3d(lev, patch, offsets=[0,1,0]) # orientation 90 degrees (example same slices: equal to Nslices*90degree 2D case)
        g4 = glcm3d(lev, patch, offsets=[0,1,1]) # orientation 135 degrees (example same slices: equal to Nslices*135degree 2D case)
        g5 = glcm3d(lev, patch, offsets=[1,0,-1]) # 0 degrees/45 degrees (example same slices: equal to (Nslices-1)*0degree 2D case)
        g6 = glcm3d(lev, patch, offsets=[1,0,0])  # straight up (example same slices: equal to np.unique())
        g7 = glcm3d(lev, patch, offsets=[1,0,1]) # 0 degree/135 degrees (example same slices: equal to (Nslices-1)*transpose of 0degree 2D case)
        g8 = glcm3d(lev, patch, offsets=[1,1,0]) # 90 degrees/45 degrees (example same slices: equal to (Nslices-1)*90 degree 2D case)
        g9 = glcm3d(lev, patch, offsets=[1,-1,0])    # 90 degrees/135 degrees (example same slices: equal to (Nslices-1)*transpose of 90 degree 2D case)
        g10 = glcm3d(lev, patch, offsets=[1,1,-1])    # 45 degrees/45 degrees (example same slices: equal to (Nslices-1)*45 degree 2D case)
        g11 = glcm3d(lev, patch, offsets=[1,-1,1])   # 45 degree/135 degrees (example same slices: equal to (Nslices-1)*transpose of 45 degree 2D case)
        g12 = glcm3d(lev, patch, offsets=[1,1,1])    # 135 degrees/45 degrees (example same slices: equal to (Nslices-1)*135 degree 2D case)
        g13 = glcm3d(lev, patch, offsets=[0,0,1])    # 135 degrees/135 degrees (example same slices: equal to (Nslices-1)*transpose of 135 degree 2D case)

        # plot                
        fig = plt.figure()
        fig.add_subplot(431);           plt.imshow(g1); plt.gray()
        fig.add_subplot(432);           plt.imshow(g2); plt.gray()
        fig.add_subplot(433);           plt.imshow(g3); plt.gray()
        fig.add_subplot(434);           plt.imshow(g4); plt.gray()
        fig.add_subplot(435);           plt.imshow(g5); plt.gray()
        fig.add_subplot(436);           plt.imshow(g6); plt.gray()
        fig.add_subplot(437);           plt.imshow(g7); plt.gray()
        fig.add_subplot(438);           plt.imshow(g8); plt.gray()
    
        # add all directions to make features non-directional
        g = g1+g2+g3+g4+g5+g6+g7+g8+g9+g10+g11+g12+g13
    
        fig.add_subplot(439);           plt.imshow(g); plt.gray()
        #plt.show()
        
        ### glcm normalization ###
        if g.sum() != 0:
            g = g.astype(float)/g.sum()
        
        ### compute auxiliary variables ###
        (num_level, num_level2) = g.shape
        I, J = np.ogrid[0:num_level, 0:num_level]
        I = 1+ np.array(range(num_level)).reshape((num_level, 1))
        J = 1+ np.array(range(num_level)).reshape((1, num_level))
        diff_i = I - np.apply_over_axes(np.sum, (I * g), axes=(0, 1))[0, 0]
        diff_j = J - np.apply_over_axes(np.sum, (J * g), axes=(0, 1))[0, 0]
        std_i = np.sqrt(np.apply_over_axes(np.sum, (g * (diff_i) ** 2),axes=(0, 1))[0, 0])
        std_j = np.sqrt(np.apply_over_axes(np.sum, (g * (diff_j) ** 2),axes=(0, 1))[0, 0])
        cov = np.apply_over_axes(np.sum, (g * (diff_i * diff_j)),axes=(0, 1))[0, 0]
        
        gxy = np.zeros(2*g.shape[0]-1)   ### g x+y
        gx_y = np.zeros(g.shape[0])  ### g x-y
        for i in xrange(g.shape[0]):
            for j in xrange(g.shape[0]):
                gxy[i+j] += g[i,j]
                gx_y[np.abs(i-j)] += g[i,j]
        
        mx_y = (gx_y*np.arange(len(gx_y))).sum()
        v = np.zeros(11)
        i,j = np.indices(g.shape)+1
        ii = np.arange(len(gxy))+2
        ii_ = np.arange(len(gx_y))
        
        ### compute descriptors ###
        v[0] = np.apply_over_axes(np.sum, (g ** 2), axes=(0, 1))[0, 0] # energy or Angular second moment
        v[1] = np.apply_over_axes(np.sum, (g * ((I - J) ** 2)), axes=(0, 1))[0, 0] # Contrast
        if std_i>1e-15 and std_j>1e-15: # handle the special case of standard deviations near zero
            v[2] = cov/(std_i*std_j)#v[2] = greycoprops(g,'correlation') # Correlation
        else:
            v[2] = 1
        v[3] = np.apply_over_axes(np.sum, (g* (diff_i) ** 2),axes=(0, 1))[0, 0]# Variance or Sum of squares
        v[4] = np.sum(g * (1. / (1. + (I - J) ** 2))) # Inverse difference moment
        v[5] = (gxy*ii).sum() # Sum average
        v[6] = ((ii-v[5])*(ii-v[5])*gxy).sum() # Sum variance
        v[7] = -1*(gxy*np.log10(gxy+ np.spacing(1))).sum() # Sum entropy
        v[8] = -1*(g*np.log10(g+np.spacing(1))).sum() # Entropy
        v[9] = ((ii_-mx_y)*(ii_-mx_y)*gx_y).sum() # Difference variance
        v[10] = -1*(gx_y*np.log10(gx_y++np.spacing(1))).sum() # Difference entropy
            
        # writing to file from row_lesionID Drow_PathRepID
        ##################################################
        # writing to file from row_lesionID Drow_PathRepID
        print "\n Append texture features for each post contrast"
        [self.T2texture_energy_nondir, self.T2texture_contrast_nondir, self.T2texture_correlation_nondir,
         self.T2texture_variance_nondir, self.T2texture_inversediffmoment_nondir, self.T2texture_sumaverage_nondir,
         self.T2texture_sumvariance_nondir, self.T2texture_sumentropy_nondir, self.T2texture_entropy_nondir,
         self.T2texture_diffvariance_nondir,self.T2texture_diffentropy_nondir] = v
              
        ##################################################
        # orgamize into dataframe
        self.textureT2Features = DataFrame( data=array([[ self.T2texture_energy_nondir, self.T2texture_contrast_nondir, self.T2texture_correlation_nondir,
         self.T2texture_variance_nondir, self.T2texture_inversediffmoment_nondir, self.T2texture_sumaverage_nondir,
         self.T2texture_sumvariance_nondir, self.T2texture_sumentropy_nondir, self.T2texture_entropy_nondir,
         self.T2texture_diffvariance_nondir,self.T2texture_diffentropy_nondir ]]), 
        columns=['texture_energy_nondir','texture_contrast_nondir','texture_correlation_nondir', 'texture_variance_nondir', 'texture_inversediffmoment_nondir', 'texture_sumaverage_nondir', 'texture_sumvariance_nondir', 
                     'texture_sumentropy_nondir', 'texture_entropy_nondir', 'texture_diffvariance_nondir', 'texture_diffentropy_nondir'])

        return self.textureT2Features