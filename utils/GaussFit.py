#!/opt/local/bin/python

__author__ = "Andrew G. Clark"
__date__ = "7 May 2014"
__maintainer__ = "Andrew G. Clark"
__email__ = "andrew.clark@curie.fr"
__JupyterNotebookBinding__ = "Varun Kapoor", "21 December 2018"

""" This script analyzes linescans and extracts cortex thickness and density from actin/membrane linescan pairs.
The script can be run in a 'pair' mode (to analyze a single linescan pair)
or 'batch' mode (to analyze multiple directories full of linescan pairs).
The mode can be specified at the bottom ("main" function).
For batch mode:
Your parent directory should contain a file called 'dir_list.dat'
with the following information in row/column form, with only space as delimiters:
sub_dir  px_size  category  ch_actin  sigma_actin
stk_1    0.05     control   1         0.119
stk_2    0.04     siRNA     2         0.220
...
The first row must contain the column headers as shown
Definitions of input parameters:
sub_dir: The name of the sub-directory containing the linescan pairs (linescan pairs must end in '...average.dat')
px_size: The pixel size for the linescans in the given sub_dir
category: The category of the experiment in each sub_dir (can be used for plotting later)
ch_actin: The actin channel (either '1' or '2'; used for extracting cortex thickness/i_c)
sigma_actin: The sigma of the point spread function for the actin channel (used for extracting h/i_c)
Note: For the sub_dir entries in the dir_list, only those directories NOT appearing in 'completed_list_v4_1.dat' will be analyzed
Output:
In each sub-directory, a list called '.../ls_data/ls_fit_data.dat' will be created containing linescan and thickness data
    -The columns are labeled according to channel number (ch1/ch2)
    -delta is always the position of the peak intensity of channel 2 (ch2.x_peak) minus ch1.x_peak
In each sub-directory, plots of the linescans and the linescans with fits (if applicable) will be saved in '.../ls_plots/'
At the end, a master list of all of the data combined is be created in the parent_directory
For 'manual' mode:
When running the script, windows will pop up sequentially to request the following information:
-Channel 1 average linescan file
-Channel 2 average linescan file
-Pixel Size
-Actin Channel
-Sigma (Actin)
These parameters are defined above.
"""



import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize, stats
import os
import math
from copy import deepcopy
from scipy.optimize import minimize

import scipy
from pathlib import Path
import pylab

def MakePath(targetdir):
    
    p = Path(targetdir)
    if not p.is_dir():
        p.mkdir(parents=True,exist_ok=True)
    
    
        
    
def ReadFit(X, I, membraneX, membraneI, Fitaround
             , psf, inisigmaguess, showaftertime,Thickness, Time, i):
    
    
        membraneimageGaussFit = Linescan(membraneX,membraneI, Fitaround, inisigmaguess)
       
        
        GaussFit = Linescan(X,I, Fitaround, inisigmaguess)
       
        
        
        
        CortexThickness = Cortex(membraneimageGaussFit,GaussFit,psf, 2)  
        CortexThickness.get_h_i_c()
        PeakActin = GaussFit.gauss_params[2]
        PeakMembrane = membraneimageGaussFit.gauss_params[2]
        PeakDiff = PeakActin - PeakMembrane 
        if CortexThickness.h is not None :
            Time.append(i)    
            if i%showaftertime == 0: 
               print('time:', i) 
               print("Membrane Fit: (Amp, Sigma, PeakPos, C)", membraneimageGaussFit.gauss_params )
               print("Actin Fit:", GaussFit.gauss_params ) 
               CortexThickness.plot_lss()
               CortexThickness.plot_fits()
               print("Thickness (nm), center cortex , cortical actin intensity (from fit)",1000*abs(CortexThickness.h), (CortexThickness.X_c), (CortexThickness.i_c))
               
               
            Thickness.append(abs(CortexThickness.h)) 
                
                  
            
        else:
          Thickness.append(0) 
          Time.append(i)  
             
           
        return Thickness, Time 

    

    
    
def ShiftFit(Block_Actin, Block_Membrane,BlockAverageActin,BlockAverageMembrane, Time_unit, Xcalibration, Fitaround
             , psf, inisigmaguess, showaftertime,Thickness, Intensity,   Time,  t):
    

  
    
    
    Shift_Actin = []
    Shift_Membrane = []
    

    if len(Block_Actin) == 0 and len(Block_Membrane) == 0:    
     
            Thickness.append(0) 
            Intensity.append(0)
            Time.append(0)   
    
    for i in range(0, len(Block_Membrane)): 
        Membrane_param, Membrane_X, Membrane_I = Block_Membrane[i]
        Actin_param, Actin_X, Actin_I = Block_Actin[i]
        
        shift = Membrane_param[2]
        
        Membrane_X = Membrane_X - shift
        Actin_X = Actin_X - shift
        
        Shift_Membrane.append([Membrane_X, Membrane_I])  
        Shift_Actin.append([Actin_X, Actin_I])
        

    
    Shift_Membrane = np.asarray(Shift_Membrane)
    oneDMembrane = np.mean(Shift_Membrane, axis = 0)
    BlockAverageMembrane.append(oneDMembrane)
    Shift_Actin = np.asarray(Shift_Actin)
    oneDActin = np.mean(Shift_Actin, axis = 0)
    BlockAverageActin.append(oneDActin)
    
    if len(Block_Actin) > 0 and len(Block_Membrane) > 0: 
     plt.plot(oneDActin[0],oneDActin[1])
     plt.plot(oneDMembrane[0],oneDMembrane[1])
     plt.title('Mean Membrane-Actin Shifted')
    
    
    
     membraneimageGaussFit = Linescan(oneDMembrane[0],oneDMembrane[1], Fitaround, inisigmaguess)
       
        
     GaussFit = Linescan(oneDActin[0],oneDActin[1], Fitaround, inisigmaguess)
    
    
    
     CortexThickness = Cortex(membraneimageGaussFit,GaussFit,psf, 2)  
     CortexThickness.get_h_i_c()
   
        

        
        
     if CortexThickness.h is not None and abs(CortexThickness.h) < 1.0E100:
            
             
               print("Membrane Fit: (Amp, Sigma, PeakPos, C)", membraneimageGaussFit.gauss_params )
               print("Actin Fit:", GaussFit.gauss_params ) 
               CortexThickness.plot_lss()
               CortexThickness.plot_fits()
               print("Thickness (nm), center cortex , cortical actin intensity (from fit)",1000*abs(CortexThickness.h), abs(CortexThickness.X_c), (CortexThickness.i_c))
               
               
               Time.append(t)
               Thickness.append(abs(CortexThickness.h)) 
               Intensity.append((oneDActin[1].max()))   
               
                     
   
     else:
                  Thickness.append(0) 
                  Intensity.append(0)
                  
                  Time.append(0)
           
        
def takeSecond(elem):
    return elem[1]


def takeFirst(elem):
    return elem[0]

def MegaBlock(Block):
    
        allImages = Block[0]
        for i in range(1, len(Block)):
          allImages = np.hstack((allImages, Block[i]))
    
        return allImages

def SelectScan(membraneimage, image, Xcalibration, N):
    
    assert(image.shape == membraneimage.shape)
    
    Scanindex = []
    Measureindex = []
    for i in range(image.shape[1]):
        X = []
        I = []
        strip = image[:image.shape[0],i]
        for j in range(strip.shape[0]):
           X.append(j * Xcalibration)
           I.append(strip[j])
        
        
        X = np.asarray(X)
        I = np.asarray(I)
        
        Scanindex.append([i, np.amax(I)])
        
    sortedList = sorted(Scanindex, key = takeSecond, reverse = True)    
    
  
    
    for i in range(N):
        index,value = sortedList[i]
        Measureindex.append([index, value])
        
        
    SortedMeasureindex =  sorted(Measureindex, key = takeSecond, reverse = True)     
    
    return SortedMeasureindex
    
 
    
def MegaFit(membraneimage, image, N, Time_unit, Xcalibration,showaftertime, Fitaround, psf, inisigmaguess, 
            Thickness, Intensity,Peak_Actin, Block_Actin, Peak_Membrane, Block_Membrane,BlockAverageActin,BlockAverageMembrane,Time, ID ):
    
    SortedMeasureindex = SelectScan(membraneimage, image, Xcalibration, N)
    
    
    assert(image.shape == membraneimage.shape)
    for i, maxintensity in SortedMeasureindex:
        X = []
        I = []
        membraneimageX = []
        membraneimageI = []
        strip = image[:image.shape[0],i]
        membraneimagestrip = membraneimage[:membraneimage.shape[0],i]
        for j in range(strip.shape[0]):
           X.append(j * Xcalibration)
           I.append(strip[j])
        
        
        X = np.asarray(X)
        I = np.asarray(I)
        for j in range(membraneimagestrip.shape[0]):    
           membraneimageX.append(j * Xcalibration)
           membraneimageI.append(membraneimagestrip[j]) 
           
        membraneimageX = np.asarray(membraneimageX)
        membraneimageI = np.asarray(membraneimageI)
     
        ID.append(int(i / Time) + 1)
        print('ID:', int(i / Time) + 1, 'Maxint:', maxintensity)
        membraneimageGaussFit = Linescan(membraneimageX,membraneimageI, Fitaround, inisigmaguess)
       
        
        GaussFit = Linescan(X,I, Fitaround, inisigmaguess)
       
      
        Block_Actin.append([GaussFit.gauss_params, X, I] )
        Block_Membrane.append([membraneimageGaussFit.gauss_params, membraneimageX, membraneimageI] )
        
        
        Shift_Actin = []
        Shift_Membrane = []
    


    
        for i in range(0, len(Block_Membrane)): 
            Membrane_param, Membrane_X, Membrane_I = Block_Membrane[i]
            Actin_param, Actin_X, Actin_I = Block_Actin[i]
        
            shift = Membrane_param[2]
        
            Membrane_X = Membrane_X - shift
            Actin_X = Actin_X - shift
        
            Shift_Membrane.append([Membrane_X, Membrane_I])  
            Shift_Actin.append([Actin_X, Actin_I])
        

    
    Shift_Membrane = np.asarray(Shift_Membrane)
    
    
    
    Shift_Actin = np.asarray(Shift_Actin)
   
    

    
    for i in range(0, len(Shift_Actin)):
      membraneimageGaussFit = Linescan(Shift_Membrane[i][0],Shift_Membrane[i][1], Fitaround, inisigmaguess)
       
        
      GaussFit = Linescan(Shift_Actin[i][0],Shift_Actin[i][1], Fitaround, inisigmaguess)
    
    
    
      CortexThickness = Cortex(membraneimageGaussFit,GaussFit,psf, 2)  
      CortexThickness.get_h_i_c()
    
        

        
        
      if CortexThickness.h is not None and abs(CortexThickness.h) < 1.0E100:
            
             if i%showaftertime == 0:   
               print("Membrane Fit: (Amp, Sigma, PeakPos, C)", membraneimageGaussFit.gauss_params )
               print("Actin Fit:", GaussFit.gauss_params ) 
               CortexThickness.plot_lss()
               CortexThickness.plot_fits()
               print("Thickness (nm), center cortex , cortical actin intensity (from fit)",1000*abs(CortexThickness.h), abs(CortexThickness.X_c), (CortexThickness.i_c))
               
               
             Thickness.append(abs(CortexThickness.h)) 
             Intensity.append((Shift_Actin[i][1].max()))   
               
                     
   
      else:
                  Thickness.append(0) 
                  Intensity.append(0)
                  
        
 
    
    oneDMembrane = np.mean(Shift_Membrane, axis = 0)
  
    oneDActin = np.mean(Shift_Actin, axis = 0)
    
    plt.plot(oneDActin[0],oneDActin[1])
    plt.plot(oneDMembrane[0],oneDMembrane[1])
    plt.title('Final Mean Membrane-Actin Shifted')
    
    
    
    membraneimageGaussFit = Linescan(oneDMembrane[0],oneDMembrane[1], Fitaround, inisigmaguess)
       
        
    GaussFit = Linescan(oneDActin[0],oneDActin[1], Fitaround, inisigmaguess)
    
    
    
    CortexThickness = Cortex(membraneimageGaussFit,GaussFit,psf, 2)  
    CortexThickness.get_h_i_c()
    CortexThickness.plot_lss()
    CortexThickness.plot_fits()
    print("Final Answer : Thickness (nm), center cortex , cortical actin intensity (from fit)",1000*abs(CortexThickness.h), abs(CortexThickness.X_c), (CortexThickness.i_c))    
    
def fit_func(x, a, sigma, mu, c ):
    """Definition of gaussian function used to fit linescan peaks.
    p = [a, sigma, mu, c].
    """
    g = a / (sigma * math.sqrt(2 * math.pi)) * scipy.exp(-(x - mu)**2 / (2 * sigma**2)) + c
    return g
    
def gauss_func(p, x):
    """Definition of gaussian function used to fit linescan peaks.
    p = [a, sigma, mu, c].
    """
    a, sigma, mu, c = p #unpacks p (for readability)
    g = a / (sigma * math.sqrt(2 * math.pi)) * scipy.exp(-(x - mu)**2 / (2 * sigma**2)) + c
    return g


def convolved(p,x):
    """Defines convolved linescan. Args: x: float or list/iterable of floats,
    the position for which convolved intensity is calculated; p: list/iterable
    of floats, linecan parameters (p=[i_in, i_c, i_out, h, x_c, sigma]).
    Returns: i: float, intensity at x.
    """
    i_in, i_c, i_out, h, x_c, sigma = p #unpacks p (for readability)

    i = (i_in + (i_c - i_in) * stats.norm.cdf((x - x_c) + h / 2., 0., sigma) +
         (i_out - i_c) * stats.norm.cdf((x - x_c) - h / 2., 0., sigma))

    return i

def unconvolved(p,x):
    """Defines unconvolved linescan. Args: x: float or list/iterable of floats,
    the position for which intensity is calculated; p: list/iterable of floats,
    linecan parameters (p=[i_in, i_c, i_out, h, x_c]). Returns: i: float,
    intensity at x.
    """

    i_in, i_c, i_out, h, x_c = p #unpacks p (for readability)

    i = np.zeros(len(x))

    for j in range(len(x)):
        if x[j] < x_c - h / 2.:
            i[j] = i_in
        if x[j] >=  x_c - h / 2. and x[j] <  x_c + h / 2.:
            i[j] = i_c
        if x[j] >= x_c + h / 2.:
            i[j] = i_out

    return i
class Linescan():
    """Linescan object with methods to extract important parameters
    from linescans.
    """

    def __init__(self,x,i, Fitaround, inisigmaguess):
        """Initializes linescan.
        Args:
            x (list of numbers): the position values
            i (list of numbers): the intensity values
        """
        #populate linescan position/intensity
        self.x = x #position list as NumPy array of floats
        self.i = i #intensity list as NumPy array of floats
        self.inisigmaguess = inisigmaguess
        #detminere a few easy parameters from position/intensity
        self.H = self.x[-1] - self.x[0]
        self.i_tot = np.trapz(self.i,self.x)
        self.Fitaround = Fitaround
        
        #populate other attributes
        self.dist_to_x_in_out = 1. #specifies how far away x_in is from the peak (in um)
        self.gauss_params = None #parameter list from Gaussian fit to find peak
        self.x_peak = None #linescan peak position
        self.i_peak = None #linescan peak intensity
        self.i_in = None #intracellular intensity
        self.i_out = None #extracellular intensity
        self.max_idx = None #index of point near linescan center with highest intensity
        self.x_fit = None #position list used for peak fitting
        self.i_fit = None #intensity list used for peak fitting
        self.i_in_x_list = None #position list used to determine self.i_in
        self.i_in_i_list = None #intensity list used to determine self.i_in
        self.i_out_x_list = None #position list used to determine self.i_out
        self.i_out_i_list = None #intensity list used to determine self.i_out
        self.x_in_upper_index = None #the index at the upper end of the region where x_in is calculated
        self.x_out_lower_index = None #the index at the lower end of the region where x_out is calculated
        self.fwhm = None #full width at half-max

        #initializes linescans and determines linescan parameters
        self.extract_ls_parameters()

    def convert_px_to_um(self):
        """Multiplies list of coordinates by pixel_size."""

        self.x = np.array([a * self.px_size for a in self.x])

    def extract_ls_parameters(self):
        """Extracts intensity and position information from linescan"""

        
        self.get_peak()
        self.get_i_in_out()
        self.get_fwhm()
    def plot_gauss(self):
        x_gauss_fit =self.x_fit
        i_gauss_fit = gauss_func(self.gauss_params,self.x_fit)
        if self.gauss_params[1] > 0:
         plt.plot(x_gauss_fit,i_gauss_fit,'b')
         
         plt.show()
    def get_peak(self):
        """Finds the peak position and intensity of a linescan by fitting
        a Gaussian near the peak.
        """
        length = len(self.i)
        #restricts fitting to near the center of the linescan
        self.max_idx = np.argmax(self.i[int(length/2- self.Fitaround):int(length/2+ self.Fitaround)])+ int(length/2)- self.Fitaround 
        self.x_fit = self.x[self.max_idx-int(self.Fitaround/2):self.max_idx+int(self.Fitaround/2)+1]
        self.i_fit = self.i[self.max_idx-int(self.Fitaround/2):self.max_idx+int(self.Fitaround/2)+1]

        #picks reasonable starting values for fit
        self.i_in_guess = np.mean(self.i[:int(self.max_idx-self.Fitaround )])
        a = (self.i[self.max_idx] - self.i_in_guess)/ 2
        sigma = self.inisigmaguess
        mu = self.x[self.max_idx]
        b = self.i_in_guess

        #perform fit with starting values
        #p0 = [a, sigma, mu, b]
        p1, sucess  = optimize.curve_fit(fit_func,self.x_fit, self.i_fit,  p0 = [a, sigma, mu, b], maxfev = 1000000)
        #optimize.leastsq(self.residuals_gauss,p0,
                                      #  args=(self.x_fit, self.i_fit),
                                       # maxfev = 1000000)
        #optimize.curve_fit(fit_func,self.x_fit, self.i_fit,  p0 = [a, sigma, mu, b], maxfev = 1000000)
        
      
        self.gauss_params = p1
        
        self.x_peak = p1[2]
        self.i_peak = gauss_func(p1, self.x_peak)

    def get_i_in_out(self):
        """Gets values for intracellular intensity (self.i_in) and
        extracellular intensity (self.i_out). The left of the linescan
        (nearer zero) is always assumed to be the intracellular side.
        Note: the i_in and i_out values are calculated to be the average value
        of the ten points out from the distance between the peak and position x away
        from the peak, where x is given by self.dist_to_x_in_out (defined in __init__).
        """

        length = len(self.i)
        x_in_upper = self.x_peak - self.dist_to_x_in_out
        x_in_upper_index = np.argmin(abs(self.x - x_in_upper))
        self.x_in_upper_index = x_in_upper_index #for use in finding total intensity for density calculation
        self.i_in_x_list = self.x[int(x_in_upper_index-10):x_in_upper_index]
        self.i_in_i_list = self.i[int(x_in_upper_index-10):x_in_upper_index]
        self.i_in = np.mean(self.i_in_i_list)

        x_out_lower = self.x_peak + self.dist_to_x_in_out
        x_out_lower_index = np.argmin(abs(self.x - x_out_lower))
        self.x_out_lower_index = x_out_lower_index #for use in finding total intensity for density calculation
        self.i_out_x_list = self.x[x_out_lower_index:int(x_out_lower_index+10)]
        self.i_out_i_list = self.i[x_out_lower_index:int(x_out_lower_index+10)]
        self.i_out = np.mean(self.i_out_i_list)

    def residuals_gauss(self,p,x,x_data):
        """Returns residuals for Gaussian fit of the intensity peak.
        Possible values for fit parameters are constrained to avoid
        overestimation of peak intensity.
        Args:
            p (list): fit parameters, [a, sigma, mu, c]
            x (list): position values
            x_data (list): intensity values
        Returns:
            residuals (list): residuals for fit
             -or-
            fail_array (list): in place of residuals if the fit fails
        """

        a, sigma, mu, c = p #unpacks p (for readability)

        i_peak_guess = gauss_func(p, mu)

        fail_array = np.ones(len(x)) * 99999.

        if all([sigma >= 0.1,
               abs(i_peak_guess - self.i[self.max_idx]) < 0.5 * self.i[self.max_idx]]):
    
            residuals = gauss_func(p,x) - x_data
            
            return residuals

        else:
            return fail_array

    def get_fwhm(self):
        """Calculates the full-width at half maximum (FWHM) of the linescan peak"""

        #determines half-max
        hm = (self.i_in + self.i_peak) / 2.
        # print hm

        # finds points closest to hm to the left of the peak
        search = self.i[:self.max_idx]
        self.left_index = (np.abs(search - hm)).argmin()
        if hm > self.i[self.left_index]:
            self.left_index_left = deepcopy(self.left_index)
            self.left_index_right = self.left_index_left
        else:
            self.left_index_right = deepcopy(self.left_index)
            self.left_index_left = self.left_index_right - 1

        #gets interpolated intensity (linear interpolation between 2 surrounding points
        m_left = (self.i[self.left_index_right] - self.i[self.left_index_left]) /  (self.x[self.left_index_right] - self.x[self.left_index_left])
        b_left = self.i[self.left_index_right] - m_left * self.x[self.left_index_right]
        x_fwhm_left = (hm - b_left) / m_left
        self.fwhm_left = [x_fwhm_left,hm]

        #finds point closest to hm to the right of the peak
        search = self.i[self.max_idx:]
        self.right_index = (np.abs(search - hm)).argmin() + self.max_idx
        if hm < self.i[self.right_index]:
            self.right_index_left = deepcopy(self.right_index)
            self.right_index_right = self.right_index_left 
        else:
            self.right_index_right = deepcopy(self.right_index)
            self.right_index_left = self.right_index_right - 1

        #gets interpolated intensity (linear interpolation between 2 surrounding points
        m_right = (self.i[self.right_index_right] - self.i[self.right_index_left]) / (self.x[self.right_index_right] - self.x[self.right_index_left])
        b_right = self.i[self.right_index_right] - m_right * self.x[self.right_index_right]
        x_fwhm_right = (hm - b_right) / m_right
        self.fwhm_right = [x_fwhm_right,hm]

        self.fwhm = x_fwhm_right - x_fwhm_left
    
    
class Cortex():
    """A Class for a cortex, with actin and membrane linescans and
     methods to determine cortex thickness and density.
    """
    def __init__(self,ch1,ch2,sigma_actin,ch_actin=2):
        """Initializes linescan pairs and remaining attributes.
            Args:
                ch1 (Linescan class): the ch1 linescan
                ch2 (Linescan class): the ch2 linescan
                sigma_actin (float): the sigma of the PSF for the actin channel
            Kwargs:
                ch_actin (int): says which channel is actin
        """
        self.ch1 = ch1
        self.ch2 = ch2
        self.sigma_actin = sigma_actin
        self.ch_actin = ch_actin

        self.delta = self.ch2.gauss_params[2] - self.ch1.gauss_params[2] #separation between ch2 and ch1 peaks

    
        self.actin = self.ch2
        self.memb = self.ch1

        
        self.h_max = 5* self.delta #maximum cortex thickness (for constraining fit)
        self.i_c_max = 5* (self.actin.i.max()) #maximum cortex intensity (for constraining fit)
        self.h = None #cortex thickness (from fit)
        self.i_c = None #cortical actin intensity (from fit)
        self.density = None #cortical actin density
        self.X_c = None #background-independent center position of the cortical actin (from fit)
        self.solution = None #solution from actin cortex thickness fit

    def get_h_i_c(self):
       

          delta = abs(self.delta)
          
          #SET STARTING VALUES FOR ROOTS AND SOLUTIONS
          self.solution = 2e+20

          #only try fitting if the peak is higher than both i_in and i_out
          if ((self.actin.i_out - self.actin.i_peak) /
                (self.actin.i_in - self.actin.i_peak))>=0:

            #loops through several different starting values for i_c and h
            for i_c_factor in np.arange(0.2,3.5,0.5):
                for h_factor in np.arange(0.2,3.5,0.5):

                    i_c_start = self.actin.i_peak * i_c_factor
                    delta_start = ((self.sigma_actin**2 / delta*2) *
                                   np.log(((self.actin.i_out - i_c_start) /
                                          (self.actin.i_in - i_c_start  ))))
                    h_start = 2 * (delta - delta_start ) * h_factor
                    #print(delta, delta_start, h_start, i_c_start)
                    #performs fit
                    p0 = [h_start, i_c_start]

                    try:
                        result = optimize.leastsq(self.residuals, p0,
                                                  maxfev=10000, full_output=1)

                        solution_temp = np.sum([x**2 for x in result[2]['fvec']])

                        if solution_temp < self.solution:
                            self.solution = deepcopy(solution_temp)
                            p1 = result[0]

                    except TypeError:
                        pass


           
                self.h, self.i_c = p1
                actin_ls_mean = np.mean(self.actin.i[:self.actin.x_out_lower_index+10])
                self.density = (self.i_c - self.actin.i_in) / actin_ls_mean
                
                if self.memb.x_peak > self.actin.x_peak:
                  self.X_c = self.memb.x_peak - self.h / 2
                else:
                  self.X_c = self.memb.x_peak + self.h / 2 
                
                
    
    def residuals(self,p):
        """Calculates residuals for cortex linescan fit to extract cortex
        thickness and intensity values
        Args:
            p (list of floats): [thickness, cortex_intensity]
        Returns:
            residuals (list of floats): [residual1, residual2]
            -or-
            fail_array (list of floats): [1000000., 1000000.]
             (returned only if fitting fails)
        """

        fail_array = [1000000., 1000000.]

        #constrains fit and ensures log term is positive
        if all([self.h_max>p[0]>0,
               self.i_c_max>p[1]>self.actin.i_in,
               (self.actin.i_out - p[1]) / (self.actin.i_in - p[1]) > 0]):

            #X_c is the position of the center of the cortex
            #x_c is the position of the cortex peak
            if self.memb.x_peak > self.actin.x_peak:
             X_c_try = self.memb.x_peak - p[0] / 2
            else:
             X_c_try = self.memb.x_peak + p[0] / 2   
            delta_try = (self.sigma_actin**2 / p[0]) * np.log((self.actin.i_out - p[1]) / (self.actin.i_in - p[1]))
            x_c_try = X_c_try - delta_try
            i_peak_try = convolved([self.actin.i_in, p[1], self.actin.i_out, p[0], X_c_try, self.sigma_actin], x_c_try)

            #residuals are difference between calculated peak position/intensity and values from data
            residuals = [x_c_try - self.actin.x_peak, i_peak_try - self.actin.i_peak]
            return residuals

        else:
            return fail_array

    def plot_lss(self):
        """Plots linescans"""

        fig = pylab.figure()
        ax = fig.add_subplot(1,1,1)

        #plots raw data
        pylab.plot(self.ch1.x,self.ch1.i,'r',label="Ch. 1")
        pylab.plot(self.ch2.x,self.ch2.i,'g',label="Ch. 2")

        #plots points used for determining i_in and i_out
        pylab.plot(self.ch1.i_in_x_list,self.ch1.i_in_i_list,'yo',label=r"$i_{\rm{in}}$, $i_{\rm{out}}$")
        pylab.plot(self.ch2.i_in_x_list,self.ch2.i_in_i_list,'yo')
        pylab.plot(self.ch1.i_out_x_list,self.ch1.i_out_i_list,'yo')
        pylab.plot(self.ch2.i_out_x_list,self.ch2.i_out_i_list,'yo')

        #plots points used to calculate fwhm and shows the fwhm
        # pylab.plot(self.ch1.x[self.ch1.left_index_left],self.ch1.i[self.ch1.left_index_left],'ko',label="fwhm points")
        # pylab.plot(self.ch1.x[self.ch1.left_index_left],self.ch1.i[self.ch1.left_index_left],'ko')
        # pylab.plot(self.ch1.x[self.ch1.left_index_right],self.ch1.i[self.ch1.left_index_right],'ko')
        # pylab.plot(self.ch1.x[self.ch1.right_index_left],self.ch1.i[self.ch1.right_index_left],'ko')
        # pylab.plot(self.ch1.x[self.ch1.right_index_right],self.ch1.i[self.ch1.right_index_right],'ko')
        #
        # pylab.plot(self.ch2.x[self.ch2.left_index_left],self.ch2.i[self.ch2.left_index_left],'ko')
        # pylab.plot(self.ch2.x[self.ch2.left_index_right],self.ch2.i[self.ch2.left_index_right],'ko')
        # pylab.plot(self.ch2.x[self.ch2.right_index_left],self.ch2.i[self.ch2.right_index_left],'ko')
        # pylab.plot(self.ch2.x[self.ch2.right_index_right],self.ch2.i[self.ch2.right_index_right],'ko')

        x_fwhm1, i_fwhm1 = zip(self.ch1.fwhm_left,self.ch1.fwhm_right)
        x_fwhm2, i_fwhm2 = zip(self.ch2.fwhm_left,self.ch2.fwhm_right)

        pylab.plot(x_fwhm1, i_fwhm1,'r',ls='-',marker='x',label="fwhm")
        pylab.plot(x_fwhm2, i_fwhm2,'g',ls='-',marker='x',label='fwhm')

        # x_fwhm1 = [self.ch1.x[self.ch1.left_index],self.ch1.x[self.ch1.right_index]]
        # y_fwhm1 = (self.ch1.i[self.ch1.left_index] + self.ch1.i[self.ch1.right_index]) / 2.
        # i_fwhm1 = [y_fwhm1,y_fwhm1]
        # pylab.plot(x_fwhm1,i_fwhm1,'g-',label="fwhm")
        #
        # x_fwhm2 = [self.ch2.x[self.ch2.left_index],self.ch2.x[self.ch2.right_index]]
        # y_fwhm2 = (self.ch2.i[self.ch2.left_index] + self.ch2.i[self.ch2.right_index]) / 2.
        # i_fwhm2 = [y_fwhm2,y_fwhm2]
        # pylab.plot(x_fwhm2,i_fwhm2,'r-',label="fwhm")

        #plots gaussian fit curve
        x_gauss_fit_ch1 = np.linspace(self.ch1.x_fit[0],self.ch1.x_fit[-1],100)
        i_gauss_fit_ch1 = gauss_func(self.ch1.gauss_params,x_gauss_fit_ch1)
        pylab.plot(x_gauss_fit_ch1,i_gauss_fit_ch1,'b',label="Peak fit")

        x_gauss_fit_ch2 = np.linspace(self.ch2.x_fit[0],self.ch2.x_fit[-1],100)
        i_gauss_fit_ch2 = gauss_func(self.ch2.gauss_params,x_gauss_fit_ch2)
        pylab.plot(x_gauss_fit_ch2,i_gauss_fit_ch2,'b')

        #finish plot
        y_min, y_max = ax.get_ylim()
        
        pylab.ylim = (y_max/2,y_max)
        x_min, x_max = ax.get_xlim()
        #pylab.axis([x_max/2-2,x_max-1, y_max/2, y_max])
        pylab.xlabel("Position ($\mu$m)")
        pylab.ylabel("Intensity (AU)")
        pylab.legend(loc='upper right')
        pylab.gcf().subplots_adjust(bottom=0.15)
        
        pylab.show()
    def plot_fits(self):
        """Plots linescan pair with fitted cortex thickness"""

        fig = pylab.figure()
        ax = fig.add_subplot(1,1,1)

        if self.ch_actin==1 or self.ch_actin=="1":
            color_actin = 'g'
            color_memb = 'r'
        elif self.ch_actin==2 or self.ch_actin=="2":
            color_actin = 'g'
            color_memb = 'r'
        else:
            raise ValueError("Please specify ch_actin as <<1>>, <<2>> for plotting fit!")

        #plots raw data
        pylab.plot(self.memb.x,self.memb.i,'o',color=color_memb,label="Memb. (raw)")
        pylab.plot(self.actin.x,self.actin.i,'o',color=color_actin,label="Actin (raw)")

        #plots unconvolved and extracted actin linescans from fits
        x_actin_hd = np.linspace(self.actin.x[0],self.actin.x[-1],1000)
        i_actin_unconv = unconvolved([self.actin.i_in, self.i_c,
                                       self.actin.i_out, self.h, self.X_c],
                                      x_actin_hd)
        i_actin_conv = convolved([self.actin.i_in, self.i_c,
                                   self.actin.i_out, self.h, self.X_c, self.sigma_actin],
                                  x_actin_hd)

        pylab.plot(x_actin_hd,i_actin_unconv,ls='-',color=color_actin)
        pylab.plot(x_actin_hd,i_actin_conv,ls='--',color=color_actin)

        pylab.axvline(x=self.memb.x_peak, color=color_memb, ls='--')

        #finishes plot
        y_min, y_max = ax.get_ylim()
        pylab.ylim = (y_max/2,y_max)
        x_min, x_max = ax.get_xlim()
        #pylab.axis([x_max/2-2,x_max-1, y_max/2, y_max])
        pylab.xlabel("Position ($\mu$m)")
        pylab.ylabel("Intensity (AU)")
        pylab.legend(loc='upper right')
        pylab.gcf().subplots_adjust(bottom=0.15)
        pylab.show()
def write_master_list(parent_dir,version):
    """Writes a master data lis in the parent directory for batch mode.
    Args:
        parent_dir (string): path of the parent directory
        version (string): the version of the software (for naming output file)
    """

    dir_list_path = parent_dir + '/dir_list.dat'
    subdir_list = [_[0] for _ in uf.read_file(dir_list_path)][1:]

    master_data = []
    for i in range(len(subdir_list)):
        data_dir = parent_dir + '/' + subdir_list[i]
        data = uf.read_file(data_dir + '/ls_data/ls_data.dat')
        if i==0:
            for line in data:
                master_data.append(line)
        else:
            for line in data[1:]:
                master_data.append(line)

    # print master_data
    uf.save_data_array(master_data, parent_dir + '/master_list_v%s.dat'%version)

def load_ls(ls_path,px_size=1.):
    """Loads a linescan file
    Args:
        ls_path (str): path of the average linescan file to be loaded
        px_size (float): pixel size in microns
    Returns:
        x (numpy array): the positions (in microns)
        i (numpy array): the intensities
    """

    ls_data = uf.read_file(ls_path)
    x = np.array([float(_[0]) for _ in ls_data]) * px_size
    i = np.array([float(_[1]) for _ in ls_data])
    return x,i

def analyze_cortex(file_ch1,file_ch2,px_size,ch_actin,sigma_actin):

    """Extracts linescan parameters and coretx thickness/density
    for a pair of linescans
    Args:
        file_ch1 (str): the filepath for the first linescan
        file_ch2 (str): the filepath for the second linescan
        px_size (float): the pixel size for the linescans (for the whole directory)
        ch_actin (int): the channel of the actin linescan (1 or 2)
        sigma_actin (float): the sigma of the PSF for the actin channel
    Kwargs:
        category (str): used to keep track of different conditions in the output data file
    Returns:
        cortex (Cortex class): the cortex with associated attributes
    """

    x_ch1, i_ch1 = load_ls(file_ch1,px_size=px_size)
    x_ch2, i_ch2 = load_ls(file_ch2,px_size=px_size)
    x = deepcopy(x_ch1) #the x values should be the same for both linescans!

    basename = file_ch1.split('/')[-1][:-4]
    print('Analyzing file pair for:', basename)

    # extracts data
    actin = Linescan(x,i_ch1)
    memb = Linescan(x,i_ch2)
    cortex = Cortex(actin, memb, sigma_actin, ch_actin=ch_actin)

    if ch_actin==1 or ch_actin==2:
        cortex.get_h_i_c()
    elif ch_actin == "None":
        pass
    else:
        raise ValueError("Please specify ch_actin as <<1>> or <<2>> for %s!"%file_ch1)

    print('h =', cortex.h)
    return cortex

def analyze_ls_pair(file_ch1,file_ch2,px_size,ch_actin,sigma_actin,version):
    """Analyzes linescans to extract cortex thickness/density
    for a single linescan pair. Data and plots are generated and saved
    to a new folder with same name as file_ch1
    Args:
        file_ch1 (str): the filepath for the first linescan
        file_ch2 (str): the filepath for the second linescan
        px_size (float): the pixel size for the linescans (for the whole directory)
        ch_actin (int): the channel of the actin linescan (1 or 2)
        sigma_actin (float): the sigma of the PSF for the actin channel
    """

    # makes directory in data_dir for saving
    save_dir = file_ch1[:-4] + '_ls_data'
    uf.make_dir(save_dir)

    # makes a list of parameters to extract from cortex data
    data_to_write = [['basename', 'category',
                      'delta', 'h', 'i_c', 'density', 'X_c', 'solution',
                      'ch1.i_tot', 'ch1.H', 'ch1.x_peak', 'ch1.i_peak', 'ch1.i_in', 'ch1.i_out', 'ch1.fwhm',
                      'ch2.i_tot', 'ch2.H', 'ch2.x_peak', 'ch2.i_peak', 'ch2.i_in', 'ch2.i_out', 'ch2.fwhm'
                      ]]

    basename = file_ch1.split('/')[-1][:-4]
    category = 'pair'

    #gets cortex and linescan data
    cortex = analyze_cortex(file_ch1, file_ch2, px_size, ch_actin, sigma_actin)

    # plots raw linescans
    cortex.plot_lss()
    pylab.savefig(save_dir + "/" + basename + ".png")
    pylab.close()

    # plots linescans with h fits
    if cortex.h != None:
        cortex.plot_fits()
        pylab.savefig(save_dir + "/" + basename + "_fit.png")
        pylab.close()

    # gets extracted linescan data
    data_temp = [basename, category]
    for param in data_to_write[0][2:]:
        data_temp.append(eval("cortex.%s" % param))
    data_to_write.append(data_temp)

    # print data_to_write
    uf.save_data_array(data_to_write, save_dir + "/ls_data.dat")

def analyze_dir(data_dir,px_size,category,ch_actin,sigma_actin,version):
    """ Analyzes all linescan pairs in a directory full of linescans
    Args:
        data_dir (str): the directory containing the linescans
        px_size (float): the pixel size for the linescans (for the whole directory)
        category (str): the category for the experiment
        ch_actin (int): the channel of the actin linescan (1 or 2)
        version (str): version number (for output filenames)
    """

    #makes necessary directories in data_dir for saving
    save_dir = data_dir + '/ls_data'
    uf.make_dir(save_dir)

    #makes a list of parameters to extract from cortex data
    data_to_write = [['basename','category',
                      'delta', 'h', 'i_c', 'density', 'X_c', 'solution',
                      'ch1.i_tot','ch1.H','ch1.x_peak','ch1.i_peak','ch1.i_in','ch1.i_out','ch1.fwhm',
                      'ch2.i_tot','ch2.H','ch2.x_peak','ch2.i_peak','ch2.i_in','ch2.i_out','ch2.fwhm'
                      ]]

    #gets and sorts list of average linescans
    linescan_list = [x for x in os.listdir(data_dir) if 'average.dat' in x]

    for _ in linescan_list:
        
        print(re.search('frame' + '_([0-9]+)_', _).group(1))
    linescan_list = sort_ls_list(linescan_list)


    #extracts linescan parameters and thickness/density
    for i in range(len(linescan_list)/2):

        file_ch1 = data_dir + '/' + linescan_list[2*i]
        file_ch2 = data_dir + '/' + linescan_list[2*i + 1]
        basename = file_ch1.split('/')[-1][:-4]

        cortex = analyze_cortex(file_ch1,file_ch2,px_size,ch_actin,sigma_actin)

        # plots raw linescans
        cortex.plot_lss()
        pylab.savefig(save_dir + "/" + basename + ".png")
        pylab.close()

        # plots linescans with h fits
        if cortex.h != None:
            cortex.plot_fits()
            pylab.savefig(save_dir + "/" + basename + "_fit.png")
            pylab.close()

        # gets extracted linescan data
        data_temp = [basename,category]
        for param in data_to_write[0][2:]:
            data_temp.append(eval("cortex.%s"%param))
        data_to_write.append(data_temp)

    # print data_to_write
    uf.save_data_array(data_to_write,save_dir + "/ls_data.dat")

    
    