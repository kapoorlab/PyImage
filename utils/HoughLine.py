import numpy as np

import math
from skimage.transform import (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
from matplotlib import cm
import matplotlib.pyplot as plt
from skimage.measure import LineModelND, ransac
from numpy.linalg import norm
from scipy import ndimage as ndi
from skimage.filters import roberts, sobel, scharr, prewitt, gaussian
from skimage.morphology import watershed
from skimage.feature import peak_local_max
from Normalize import save_tiff_imagej_compatible
from skimage.segmentation import find_boundaries,find_boundaries, relabel_sequential
from skimage.morphology import remove_small_objects, binary_erosion
from skimage.filters import threshold_otsu, threshold_mean
from skimage.exposure import rescale_intensity
import os
from scipy.ndimage import gaussian_filter
from sklearn.preprocessing import PolynomialFeatures 
from bokeh.io import export_png, output_notebook
from bokeh.plotting import figure, output_file, show
from sklearn.cluster import KMeans
from skimage import feature, filters
from scipy import ndimage



def show_poly_regression(X, Y, degree = 2):    
 
 z = np.polyfit(X, Y, degree)
 ffit = np.poly1d(z)
 xp = np.linspace(-2, 6, 100)
 print('Coefficients (High -> Low)', z) 
 x_new = np.linspace(X[0], X[-1], 50)
 y_new = ffit(x_new)
 plt.plot(X,Y, x_new, y_new)
 plt.title('Polynomial Fit') 
 plt.xlabel('Time') 
 plt.ylabel('Divisions') 
  
 plt.show()


 p = figure(title='Division Counter', x_axis_label='Time', y_axis_label='Division Number')
 p.line(x_new, y_new, legend = "Divisions-Time", line_width = 2)
   
 output_notebook()
 show(p)


def rgb2gray(rgb):
    return np.dot(rgb[..., :3], [0.299, 0.587, 0.114]).astype(np.uint8)

def show_plot(points,  ymin, ymax):

    
    fig, ax = plt.subplots() 
    ax.plot(points[:, 0], points[: , 1], '.b', alpha=0.6,
        label='Inlier data')
    x_min, x_max = ax.get_xlim()
    ax.axis([x_min,x_max, ymin, ymax])
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Thickness (um)')
    plt.show()
 
def show_intensity_plot(points,save_dir,name,   ymin, ymax, num_clusters, title = None  ):

    
    fig, ax = plt.subplots() 
    ax.plot(points[:, 1], points[: , 0], '.b', alpha=0.6,
        label='Inlier data')
    x_min, x_max = ax.get_xlim()
    ax.axis([x_min,x_max, ymin, ymax])
    if title is None:
     ax.set_xlabel('Intensity')
    else:
     ax.set_xlabel(title)   
    ax.set_ylabel('Thickness (um)')
    X = np.column_stack([points[:, 1], points[: , 0]])  
    kmeans = KMeans(n_clusters=num_clusters) # You want cluster the passenger records into 2
    kmeans.fit(X)
    y_kmeans = kmeans.predict(X)
    plt.scatter(X[:,0], X[:,1], c=y_kmeans, s=50, cmap='viridis')
    centers = kmeans.cluster_centers_
    plt.scatter(centers[:, 0], centers[:, 1], c='black', s=200, alpha=0.5);
    for i in range(0, len(centers)):
     print('X:', centers[i, 0], 'Y: ', centers[i, 1])
     distances = compute_distance(X, centers, num_clusters)
     print('Standard deviation:', np.mean(distances))
    
    #if title is not None:
     #plt.savefig(save_dir + "/" + name + " Thickness-" + title + '.png')
    #else:
     #plt.savefig(save_dir + "/" + name + " Thickness-" + "Intensity" + '.png')  
    plt.show()
    
    
def compute_distance(X, centroids, n_clusters):
        distance = np.zeros((X.shape[0], n_clusters))
        for k in range(n_clusters):
            row_norm = norm(X - centroids[k, :], axis=1)
            distance[:, k] = np.square(row_norm)
        return distance
def Correlation_plot(pointsA, pointsB,x_min,x_max,y_min,y_max, id):
    
    fig, ax = plt.subplots()
    
    ax.axis([x_min,x_max, y_min, y_max])
    ax.plot(pointsB, pointsA, '.b', alpha = 0.6, label = 'Correlation plot')
    for i in range(0, len(id)):
     currentid = id[i]
     try:
           x = pointsB[i]
           y = pointsA[i]
          
     except IndexError:
            continue
     if(x < x_max and x > x_min and y< y_max and y > y_min ):       
      ax.text(pointsB[i], pointsA[i], str(currentid)) 
     ax.set_xlabel('Intensity')
     ax.set_ylabel('Thickness (um)')
    
    plt.show()

def Peak_Data(pointsA, pointsB, title):
    
    

    
    
    fig, ax = plt.subplots()
    
    ax.plot(pointsB, pointsA, '.b', alpha = 0.6, label = 'Peak-Data' + title)
    ax.set_xlabel('Data')
    ax.set_ylabel('Peak'+ title)
    
    plt.show()    
    
def show_ransac_points_line(points,  min_samples=2, residual_threshold=0.1, max_trials=1000, Xrange = 100, displayoutlier= False):
   
    # fit line using all data
 model = LineModelND()
 if(len(points) > 2):
  model.estimate(points)
 
  fig, ax = plt.subplots()   

  # robustly fit line only using inlier data with RANSAC algorithm
  model_robust, inliers = ransac(points, LineModelND, min_samples=min_samples,
                               residual_threshold=residual_threshold, max_trials=max_trials)
  slope , intercept = model_robust.params
 
  outliers = inliers == False
  # generate coordinates of estimated models
  line_x = np.arange(0, Xrange)
  line_y = model.predict_y(line_x)
  line_y_robust = model_robust.predict_y(line_x)
 
  #print('Model Fit' , 'yVal = ' , line_y_robust)
  #print('Model Fit', 'xVal = ' , line_x)
  ax.plot(points[inliers, 0], points[inliers, 1], '.b', alpha=0.6,
        label='Inlier data')
  if displayoutlier:
   ax.plot(points[outliers, 0], points[outliers, 1], '.r', alpha=0.6,
        label='Outlier data')
  #ax.plot(line_x, line_y, '-r', label='Normal line model')
  ax.plot(line_x, line_y_robust, '-g', label='Robust line model')
  ax.legend(loc='upper left')
   
  ax.set_xlabel('Time (s)')
  ax.set_ylabel('Thickness (um)')
  print('Ransac Slope = ', str('%.3e'%((line_y_robust[Xrange - 1] - line_y_robust[0])/ (Xrange)) )) 
  print('Regression Slope = ', str('%.3e'%((line_y[Xrange - 1] - line_y_robust[0])/ (Xrange)) )) 
  print('Mean Thickness (After outlier removal) = ', str('%.3f'%(sum(points[inliers, 1])/len(points[inliers, 1]))), 'um')   
  plt.show()
 
    
    

def show_ransac_line(img, Xcalibration, Time_unit, maxlines, min_samples=2, residual_threshold=0.1, max_trials=1000):
    points = np.array(np.nonzero(img)).T

    f, ax = plt.subplots(figsize = (10, 10))

    points = points[:, ::-1]

    for i in range(maxlines):
  
     
    # robustly fit line only using inlier data with RANSAC algorithm
     model_robust, inliers = ransac(points, LineModelND,  min_samples=min_samples,
                               residual_threshold=residual_threshold, max_trials=max_trials)
     slope , intercept = model_robust.params
 
     points = points[~inliers]   

     print("Estimated Wave Velocity by Ransac : " , np.abs(slope[0])* (Xcalibration / Time_unit)) 
     x0 = np.arange(img.shape[1])   
 
     y0 =  model_robust.predict_y(x0)
     y_min, y_max = ax.get_ylim()
     x_min, x_max = ax.get_xlim()
     
     ax.axis = ([x_min,x_max, y_max/2, y_max])
     plt.plot(x0, model_robust.predict_y(x0), '-r')
 
    ax.imshow(img)
    
def watershed_binary(image, size, gaussradius, kernel, peakpercent):
 
 
 distance = ndi.distance_transform_edt(image)

 gauss = gaussian_filter(distance, gaussradius)

 local_maxi = peak_local_max(gauss, indices=False, footprint=np.ones((kernel, kernel)),
                            labels=image)
 markers = ndi.label(peakpercent * local_maxi)[0]
 labels = watershed(-distance, markers, mask=image)


 nonormimg = remove_small_objects(labels, min_size=size, connectivity=4, in_place=False)
 nonormimg, forward_map, inverse_map = relabel_sequential(nonormimg)    
 labels = nonormimg

    
    
 
 #plt.imshow(labels)
 #plt.title('Watershed labels')   
 #plt.show()
 return labels   
def watershed_image(image, size, targetdir, Label, Filename, Xcalibration,Time_unit):
 distance = ndi.distance_transform_edt(image)
 
 plt.imshow(distance)
 plt.title('Distance transform')   
 plt.show()  
 
 local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((1, 1)),
                            labels=image)
 markers = ndi.label(local_maxi)[0]
 labels = watershed(-distance, markers, mask=image)

 nonormimg = remove_small_objects(labels, min_size=size, connectivity=4, in_place=False)
 nonormimg, forward_map, inverse_map = relabel_sequential(nonormimg)    
 labels = nonormimg

    
    
 
 plt.imshow(labels)
 plt.title('Watershed labels')   
 plt.show()
 print('Doing Hough in +' , np.unique(labels) , 'labels')
 Velocity = []
 Images = []
 Besty0 = []
 Besty1 = []
 # loop over the unique labels returned by the Watershed
 # algorithm
 for label in np.unique(labels):
      
      if label== 0:
            continue
     
      mask = np.zeros(image.shape, dtype="uint8")
      mask[labels == label] = 1
     
          
      h, theta, d = hough_line(mask)  
      img, besty0, besty1, velocity = show_hough_linetransform(mask, h, theta, d, Xcalibration, 
                               Time_unit,targetdir, Filename[0])

      if np.abs(velocity) > 1.0E-5:  
       Velocity.append(velocity)
       Images.append(img)
       Besty0.append(besty0)
       Besty1.append(besty1)
 return Velocity, Images, Besty0, Besty1    

    
def show_hough_linetransform(img, accumulator, thetas, rhos, Xcalibration, Tcalibration,  save_path=None, File = None):
    import matplotlib.pyplot as plt

    #fig, ax = plt.subplots(1, 2, figsize=(10, 10))

   

    #ax[0].imshow(
        #accumulator, cmap=cm.gray,
        #extent=[np.rad2deg(thetas[-1]), np.rad2deg(thetas[0]), rhos[-1], rhos[0]])
    
    #ax[0].set_title('Hough transform')
    #ax[0].set_xlabel('Angles (degrees)')
    #ax[0].set_ylabel('Distance (pixels)')
    #ax[0].axis('image')
    #ax[1].imshow(img, cmap=cm.gray)
    
    bestpeak = 0
    bestslope = 0
    besty0 = 0
    besty1 = 0
    Est_vel = []
    for _, angle, dist in zip(*hough_line_peaks(accumulator, thetas, rhos)):
     y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
     y1 = (dist - img.shape[1] * np.cos(angle)) / np.sin(angle)
    
     pixelslope =   -( np.cos(angle) / np.sin(angle) )
     pixelintercept =  dist / np.sin(angle)  
     slope =  -( np.cos(angle) / np.sin(angle) )* (Xcalibration / Tcalibration)
     
    #Draw high slopes
     peak = 0;
     for index, pixel in np.ndenumerate(img):
            x, y = index
            vals = img[x,y]
            if  vals > 0:
                peak+=vals
                if peak >= bestpeak:
                    bestpeak = peak
                    bestslope = slope
                    besty0 = y0
                    besty1 = y1
   
    
    
    #ax[1].plot((0, img.shape[1]), (besty0, besty1), '-r')
    
    #ax[1].set_xlim((0, img.shape[1]))
    #ax[1].set_ylim((img.shape[0], 0))
    #ax[1].set_axis_off()
    #ax[1].set_title('Detected lines')

    # plt.axis('off')
    if save_path is not None and File is not None:
       plt.savefig(save_path + 'HoughPlot' + File + '.png')
    if save_path is not None and File is None:
        plt.savefig(save_path + 'HoughPlot' + '.png')
  
    
       #plt.show()

    return (img,besty0, besty1, bestslope)    