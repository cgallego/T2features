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
class features_T2(object):
    """
    USAGE:
    =============
    T2 = features_T2()
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
        
        
    def extract_muscleSI(self, T2image, image_pos_pat, image_ori_pat, origin, iren, ren, picker, xplane, yplane, zplane):
        """ extract_muscleSI: Place manually a widget over muscle, extract reference SI
        
        INPUTS:
        =======        
        images: (vtkImageData)   list of Input image to Transform
        origin: float[3]    list of origin coordinates
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
        
        return muscleSI, muscle_scalar_range


    def extract_lesionSI(self, T2image, lesion3D, image_pos_pat, image_ori_pat):
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
        [t_T2image, transform_cube] = loadDisplay.dicomTransform(T2image, image_pos_pat, image_ori_pat)
            
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
        
        
        lesion_scalar_range = [lesionSI.min(), lesionSI.max()]
        print "\nLesion_scalar Range:"
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
        
        
    def extractT2morphology(self, T2image, VOI_mesh, image_pos_pat, image_ori_pat):
        """ Start pixVals for collection pixel values at VOI """
        pixVals_margin = []; pixVals = []
        Fmargin = {}; voxel_frameS = {}
        
        # necessary to read point coords
        VOIPnt = [0,0,0]
        ijk = [0,0,0]
        pco = [0,0,0]
        
        ## Transform T2 img
        loadDisplay = Display()
        
        # Proceed to build reference frame for display objects based on DICOM coords   
        [t_T2image, transform_cube] = loadDisplay.dicomTransform(T2image, image_pos_pat, image_ori_pat)
        # Update info
        t_T2image.UpdateInformation() 
        
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
            
            # Compute dot product (unit vector for dot product)
            p_dot = dot(grad_marginpt, norm_rdir)
            norm_p_dot = linalg.norm(p_dot)
            
            H_norm_p.append(norm_p_dot)    
        
        # The histogram of radial gradient values quantifying the frequency of occurrence of the dot products in a given region of interest
        # radial gradient histogram. The hist() function now has a lot more options
        # first create a single histogram
                    
        # the histogram of the data with histtype='step'
        plt.figure()
        n, bins, patches = plt.hist(array(H_norm_p), 50, normed=1, histtype='bar',facecolor='blue', alpha=0.75)
        n, min_max, mean_bins, var_bins, skew, kurt = stats.describe(bins)
        
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
        
    
    def extractT2texture(self, T2image, VOI_mesh, image_pos_pat, image_ori_pat):
        ################################### 
        # Haralick et al. defined 10 grey-level co-occurrence matrix (GLCM) enhancement features 
        # (energy, maximum probability, contrast, homogeneity, entropy, correlation, sum average, sum variance, difference average and difference variance) to describe texture
        # For both mass and non-mass lesions
        averg_texture_contrast=array([0,0,0,0]).reshape(1,4)    
        averg_texture_homogeneity=array([0,0,0,0]).reshape(1,4)
        averg_texture_dissimilarity=array([0,0,0,0]).reshape(1,4)
        averg_texture_correlation=array([0,0,0,0]).reshape(1,4)
        averg_texture_ASM=array([0,0,0,0]).reshape(1,4)
        averg_texture_energy=array([0,0,0,0]).reshape(1,4)
        
        # N is the number of distinct gray levels in the histogram equalized image;
        # obtain vols of interest
        ## Transform T2 img
        loadDisplay = Display()
        
        # Proceed to build reference frame for display objects based on DICOM coords   
        [t_T2image, transform_cube] = loadDisplay.dicomTransform(T2image, image_pos_pat, image_ori_pat)
        # Update info
        t_T2image.UpdateInformation() 
        
        print "\nLoading VOI_mesh MASK... "
        VOIPnt = [0,0,0]; pixVals_margin = [];  pixVals = []
            
        #################### HERE GET INTERNAL PIXELS IT AND MASK IT OUT
        VOI_scalars = t_T2image.GetPointData().GetScalars()
        np_VOI_imagedata = vtk_to_numpy(VOI_scalars)     
        
        dims = t_T2image.GetDimensions()
        spacing = t_T2image.GetSpacing()
        np_VOI_imagedata = np_VOI_imagedata.reshape(dims[2], dims[1], dims[0]) 
        np_VOI_imagedata = array(np_VOI_imagedata.transpose(2,1,0))
        
        # Prepare lesion localization and PATCH size for Texture analysis    
        # lesion_centroid        
        lesionBox = VOI_mesh.GetBounds()
        deltax = lesionBox[1] - lesionBox[0]
        print deltax 
        deltay = lesionBox[3]-lesionBox[2]
        print  deltay
        deltaz = lesionBox[5]-lesionBox[4]
        print deltaz
        
        # find largest dimension and use that for lesion radius
        if deltax > deltay:
            if deltax > deltaz:
                lesion_dia = deltax 
                print "deltax was largest %d... " % lesion_dia
            else: 
                lesion_dia = deltaz
                print "deltaz was largest %d... " % lesion_dia
        else:
            if deltay > deltaz:
                lesion_dia = deltay
                print "deltay was largest %d... " % lesion_dia
            else:
                lesion_dia = deltaz
                print "deltaz was largest %d... " % lesion_dia
              
        lesion_radius = lesion_dia/2
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
    
        # compute texture cdf        
        eq_numpy_pre_VOI_imagedata, cdf = self.histeq(np_VOI_imagedata[:,:,int(ijk_vtk[2])])    
        plt.figure()
        plt.subplot(231)
        n, bins, patches = plt.hist(array(np_VOI_imagedata[:,:,int(ijk_vtk[2])].flatten()), 50, normed=1, facecolor='green', alpha=0.75)
        
        plt.subplot(233)
        n, bins, patches = plt.hist(array(eq_numpy_pre_VOI_imagedata.flatten()), 50, normed=1, facecolor='green', alpha=0.75)
         
        plt.subplot(234)
        plt.imshow(np_VOI_imagedata[:,:,int(ijk_vtk[2])])
        plt.gray()
        
        plt.subplot(236)
        plt.imshow(eq_numpy_pre_VOI_imagedata)
        plt.gray()
        
        # FInally display
        plt.show()
        
        #****************************"
        # Do histogram cast
        if( np_VOI_imagedata.max() > 256):
            # Do only VOI histo equ
            np_VOI_imagedata = np_VOI_imagedata.astype('uint8')
            eq_numpy_pre_VOI_imagedata, cdf = self.histeq(np_VOI_imagedata[:,:,int(ijk_vtk[2])])
        else:
            # Do only VOI histo equ
            np_VOI_imagedata = np_VOI_imagedata.astype('uint8')
            eq_numpy_pre_VOI_imagedata, cdf = self.histeq(np_VOI_imagedata[:,:,int(ijk_vtk[2])])
        
        
        # Get the shape of the numpy image data to confirm dimensions
        eq_numpy_pre_VOI_imagedata = eq_numpy_pre_VOI_imagedata.reshape(dims[0], dims[1], 1)
        
        # Perform texture classification using grey level co-occurrence matrices (GLCMs).
        # A GLCM is a histogram of co-occurring greyscale values at a given offset over an image.
        # compute some GLCM properties each patch
        # p(i,j) is the (i,j)th entry in a normalized spatial gray-level dependence matrix; 
        lesion_patches = []
        lesion_patches = np_VOI_imagedata[
                int(ijk_vtk[0] - lesion_radius):int(ijk_vtk[0] + lesion_radius),
                int(ijk_vtk[1] - lesion_radius):int(ijk_vtk[1] + lesion_radius),
                int(ijk_vtk[2])]
        
        patches_shape = lesion_patches.shape
                            
        for k in range(patches_shape[0]):
            for l in range(patches_shape[1]):
                if (lesion_patches[k,l] < 0):
                    lesion_patches[k,l] = 0
        
        print '\n Lesion_patches:'
        print lesion_patches

        #skimage.feature.greycomatrix(image, distances, angles, levels=256, symmetric=False, normed=False)
        glcm = greycomatrix(lesion_patches, [3], [0, 45*pi/180, 90*pi/180, 135*pi/180], 256, symmetric=True, normed=True)
        texture_contrast = greycoprops(glcm, 'contrast')
        texture_homogeneity = greycoprops(glcm, 'homogeneity')
        texture_dissimilarity = greycoprops(glcm, 'dissimilarity')
        texture_correlation = greycoprops(glcm, 'correlation')
        texture_ASM = greycoprops(glcm, 'ASM')
        texture_energy = greycoprops(glcm, 'energy')
    
        # display the image patches
        plt.subplot(3, len(lesion_patches), len(lesion_patches) * 1 )
        plt.imshow(lesion_patches, cmap=plt.cm.gray, interpolation='nearest',
               vmin=0, vmax=255)
        plt.xlabel('lesion_patches')
        
        # display original image with locations of patches
        plt.subplot(3, 2, 1)
        plt.imshow(np_VOI_imagedata[:,:,int(ijk_vtk[2])] , cmap=plt.cm.gray, interpolation='nearest',
               vmin=0, vmax=255)
        # Plot
        # create the figure
        plt.plot(ijk_vtk[0] - lesion_radius/2, ijk_vtk[1] - lesion_radius/2, 'gs')
        plt.xlabel('Original Image')
        plt.xticks([])
        plt.yticks([])
        plt.axis('image')
        
        # FInally display
        plt.show()
        
        #plt.close()     
            
        # writing to file from row_lesionID Drow_PathRepID
        print "\n Average texture features for each orientation"
        contrast = texture_contrast
        print contrast
        [self.contrast_zero, self.contrast_quarterRad, self.contrast_halfRad, self.contrast_halfRad] =  contrast[0,0], contrast[0,1], contrast[0,2], contrast[0,3]
        
        homogeneity = texture_homogeneity
        print homogeneity                
        [self.homogeneity_zero, self.homogeneity_quarterRad, self.homogeneity_halfRad, self.homogeneity_threeQuaRad] = homogeneity[0,0], homogeneity[0,1], homogeneity[0,2], homogeneity[0,3]

        dissimilarity =texture_dissimilarity
        print dissimilarity                
        [self.dissimilarity_zero, self.dissimilarity_quarterRad, self.dissimilarity_halfRad, self.dissimilarity_threeQuaRad] = dissimilarity[0,0], dissimilarity[0,1], dissimilarity[0,2], dissimilarity[0,3]
        
        correlation = texture_correlation
        print correlation  
        [self.correlation_zero, self.correlation_quarterRad, self.correlation_halfRad, self.correlation_threeQuaRad] = correlation[0,0], correlation[0,1], correlation[0,2], correlation[0,3]
        
        ASM = texture_ASM
        print ASM  
        [self.ASM_zero, self.ASM_quarterRad, self.ASM_halfRad, self.ASM_threeQuaRad] = ASM[0,0], ASM[0,1], ASM[0,2], ASM[0,3]
        
        energy = texture_energy
        print energy  
        [self.energy_zero, self.energy_quarterRad, self.energy_halfRad, self.energy_threeQuaRad] = energy[0,0], energy[0,1], energy[0,2], energy[0,3]

              
        ##################################################
        # orgamize into dataframe
        self.textureT2Features = DataFrame( data=array([[ self.contrast_zero, self.contrast_quarterRad, self.contrast_halfRad, self.contrast_halfRad, 
                                                   self.homogeneity_zero, self.homogeneity_quarterRad, self.homogeneity_halfRad, self.homogeneity_threeQuaRad,
                                                   self.dissimilarity_zero, self.dissimilarity_quarterRad, self.dissimilarity_halfRad, self.dissimilarity_threeQuaRad,
                                                   self.correlation_zero, self.correlation_quarterRad, self.correlation_halfRad, self.correlation_threeQuaRad,
                                                   self.ASM_zero, self.ASM_quarterRad, self.ASM_halfRad, self.ASM_threeQuaRad,
                                                   self.energy_zero, self.energy_quarterRad, self.energy_halfRad, self.energy_threeQuaRad ]]), 
        columns=['T2texture_contrast_zero', 'T2texture_contrast_quarterRad', 'T2texture_contrast_halfRad', 'T2texture_contrast_threeQuaRad', 'T2texture_homogeneity_zero', 'T2texture_homogeneity_quarterRad', 'T2texture_homogeneity_halfRad', 'T2texture_homogeneity_threeQuaRad', 'T2texture_dissimilarity_zero', 'T2texture_dissimilarity_quarterRad', 'T2texture_dissimilarity_halfRad', 'T2texture_dissimilarity_threeQuaRad', 'T2texture_correlation_zero', 'T2texture_correlation_quarterRad', 'T2texture_correlation_halfRad', 'T2texture_correlation_threeQuaRad', 'T2texture_ASM_zero', 'T2texture_ASM_quarterRad', 'T2texture_ASM_halfRad', 'T2texture_ASM_threeQuaRad', 'T2texture_energy_zero', 'T2texture_energy_quarterRad', 'T2texture_energy_halfRad', 'T2texture_energy_threeQuaRad'])


        
        return self.textureT2Features