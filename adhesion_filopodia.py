"""iCLOTS is a free software created for the analysis of common hematology workflow image data

Author: Meredith Fay, Lam Lab, Georgia Institute of Technology and Emory University
Last updated: 2022-07-29
This script corresponds to tools available in version 1.0b1, more recent implementations of tools
may be available within the iCLOTS software and in source code at github.com/iCLOTS

Script that analyzes static, fluorescence microscopy images of cells (platelets or WBCs)
adhered to some surface specifically for a single-cell resolution filopodia count and characterization

Script relies heavily on:
--skimage library region property analysis, tutorial available at:
https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_regionprops.html
--skimage library Harris corner detection, tutorial available at:
https://scikit-image.org/docs/stable/auto_examples/features_detection/plot_corner.html

Input files:
--Script is designed to work with a folder of image files (.jpg, .png, and/or .tif)
----The same input parameters are applied to each image

Input parameters:
--umpix: The ratio of microns (1e-6 m) to pixels for the image
----Use umpix = 1 for no conversion
--main_channel: 'r' for red/layer 1, 'g' for green/layer 2, 'b' for blue/layer 3
----Value must match 'r', 'g', or 'b' string format exactly
----Typically represents area of cell, like a membrane stain
----In the current version of this script, no quantification of a secondary "functional stain" is included
-----but this may be added in future versions of iCLOTS
--main_threshold: integer between 0 (black), 255 (white) to be used for main channel threshold
----Any value below this threshold becomes 0, or background
----Any value greater than or equal to this threshold becomes 1, or signal to further quantify
--min_area: Minimum area (pixels) of a region (ideally, cell) to be quantified
----Can be used to filter out obvious noise
--max_area: Maximum area (pixels) of a region to be quantified
----Can be used to filter out debris or cell clusters

Harris corner detection parameters
--k (a.u.): Harris detector free parameter
----Sharpness parameter ranging from 0 to 2, with 0 indicating you'd like the most defined filopodia only
--tr (a.u.): Corner peak finding relative threshold
----Minimum intensity of peaks, calculated as max(image) * thresh_relative
--min_distance: Minimum distance between detected filopodia
----Used with a "peak finding" algorithm


Output files
--Each original image with each detected cell labeled with an index
--Each original image with threshold applied, each detected cell is labeled with the same index
--In each image type each filopdia counted is indicated with a small circle

--A corresponding .xlsx sheet containing:
----Area, circularity, texture, filopdia count, min/mean/max/stdev length of individual filopodia (if any)
------Length of filopodia is calculated as distance of detected end point to center point of cell
--------Future versions of this code will give length of each as a vector
------For each cell - one sheet/image
------For all cells - one sheet (this combined sheet is best for analyzing replicates)
----Descriptive statistics (min, mean, max, standard deviation for all metrics)
------For each image and for all images combined
------Descriptive statistics sheet also includes a density measurement (n cells/mm^2)
----Parameters used and time/date analysis was performed, for reference

--Pairplot
----Three types:
------For each image
------For all images, combined, one color represents all
------For all images, combined, each color represents a different image

Some tips from the iCLOTS team:

Computational and experimental methods:
--We suggest high microscope magnification for this application, iCLOTS was tested on 100x magnification images
--For all fluorescent applications, each stain to quantify must be solely in one RGB channel, no other colors
----See export options on your microscopy acquisition software
--After application of the thresholds, skimage analyzes each interconnected region of 1/signal as a cell
--Analysis methods cannot distinguish between overlapping cells
----If cells are significantly overlapping, repeat experiment with a lower cell concentration
--Searching for number filopdia can be computationally expensive

Choosing parameters:
--Be sure to use microns-to-pixel ratio, not pixel-to-micron ratio
--Sometimes cells (especially platelets) have a high-intensity "body" and low-intensity spreading or protrusions
----Choose a low threshold since by counting filopodia you're primarily quantifying the morphology of cells
--Depending on the quality of your image, you may want to use an adaptive threshold rather than a binary threshold
----Line to edit indicated within script
----See OpenCV tutorial: https://docs.opencv.org/4.x/d7/d4d/tutorial_py_thresholding.html
--Err on the high side of max_area and low side of min_area parameters
---unless data is particularly noisy or there's a large amount of debris
----By sorting in excel, you can verify that small/large objects are/are not cells
--It can be tricky to adjust all 3 parameters to get a roughly accurate filopodia count
----We suggest doing robust sensitivity analysis, or you can edit them with real-time feedback in iCLOTS

Output files:
--Analysis files are named after the folder containing all images (.xlsx) or image names (.png)
----Avoid spaces, punctuation, etc. within file names
--.xlsx and pairplot data includes a sheet/graph with all images combined
----Only use this when analyzing replicates of the same sample
----No fluorescence intensity metrics are reported from the main color, as this color should indicate morphology only

"""

# Import statements
#    File management
from tkinter import filedialog
import os
import glob
import datetime
#   Computer vision, image analysis
import cv2 # Computer vision/image processing
from skimage import measure, img_as_float  # Region analysis
from skimage.feature import corner_harris, corner_peaks
import math
#   Number, file, and image management
import numpy as np
import pandas as pd  # For dataframe management
#   Plotting results
import matplotlib.pyplot as plt
import seaborn as sns  # Pairplots

cvfont = cv2.FONT_HERSHEY_SIMPLEX  # Not an import, but used for labeling images
top, bottom, left, right = [10] * 4  # Used for creating border around individual cell images


# IMPORTANT: PARAMETERS TO EDIT
umpix = 0.166  # 1 = no conversion
main_channel = 'r'  # 'r' for red/layer 1, 'g' for green/layer 2, 'b' for blue/layer 3
main_threshold = 50  # Integer between 0 (black), 255 (white)
min_area = 10
max_area = 2000
min_distance = 10  # Minimum distance between detected filopodia
k = 0.05  # Sharpness parameter ranging from 0 to 2, with 0 indicating you'd like the most defined filopodia only
tr = 0.5  # Relative threshold, minimum intensity of peaks, calculated as max(image) * threshold_rel
kernel = np.ones((5,5))  # Clean up image using morphological closing

# Select directory of files
dirpath = filedialog.askdirectory()

# Create a directory for saved results including time at which operation was performed
now = datetime.datetime.now()
# Create strings to indicate operations performed
output_folder = os.path.join(dirpath, 'Analysis, ' + now.strftime("%m_%d_%Y, %H_%M_%S"))
os.mkdir(output_folder)
os.chdir(output_folder)

# Create a list of all image files
imglist_png = sorted(glob.glob(dirpath + "/*.png"))
imglist_jpg = sorted(glob.glob(dirpath + "/*.jpg"))
imglist_tif = sorted(glob.glob(dirpath + "/*.tif"))
imglist = imglist_png + imglist_jpg + imglist_tif

# Set up method to write final excel files
dir_name = os.path.basename(dirpath)  # Also used for final graphs
excel_name = dir_name + '_analysis.xlsx'
writer = pd.ExcelWriter(excel_name, engine='openpyxl')

# Analyze images
df_all = pd.DataFrame()  # To combine all data points into one sheet, useful for replicates
df_summary = pd.DataFrame()  # For descriptive statistics

def descriptive_statistics(df_input, img_size):
    """Function to calculate descriptive statistics for each population, represented as a dataframe"""

    dict = {'n cells': len(df_input),
                  u'Min. area (\u03bcm\u00b2)': df_input[u'Area (\u03bcm\u00b2)'].min(),
                  u'Mean area (\u03bcm\u00b2)': df_input[u'Area (\u03bcm\u00b2)'].mean(),
                  u'Max. area (\u03bcm\u00b2)': df_input[u'Area (\u03bcm\u00b2)'].max(),
                  u'Stdev, area (\u03bcm\u00b2)': df_input[u'Area (\u03bcm\u00b2)'].std(),
                  'Min. area (pix)': df_input['Area (pix)'].min(),
                  'Mean area (pix)': df_input['Area (pix)'].mean(),
                  'Max. area (pix)': df_input['Area (pix)'].max(),
                  'Stdev, area (pix)': df_input['Area (pix)'].std(),
                  'Min. circularity (a.u.)': df_input['Circularity (a.u.)'].min(),
                  'Mean circularity (a.u.)': df_input['Circularity (a.u.)'].mean(),
                  'Max. circularity (a.u.)': df_input['Circularity (a.u.)'].max(),
                  'Stdev, circularity (a.u.)': df_input['Circularity (a.u.)'].std(),
                  'Min. texture (a.u.)': df_input['Texture (a.u.)'].min(),
                  'Mean texture (a.u.)': df_input['Texture (a.u.)'].mean(),
                  'Max. texture (a.u.)': df_input['Texture (a.u.)'].max(),
                  'Stdev, texture (a.u.)': df_input['Texture (a.u.)'].std(),
                  'Min. filopodia (n)': df_input['Filopodia (n)'].min(),
                  'Mean filopodia (n)': df_input['Filopodia (n)'].mean(),
                  'Max. filopodia (n)': df_input['Filopodia (n)'].max(),
                  'Stdev, filopodia (n)': df_input['Filopodia (n)'].std(),
                  'Field of view (mm\u00b2)': img_size,
                  'Cell density (n/mm\u00b2)': len(df_input)/img_size
                  }
    dict_df = pd.DataFrame(dict, index=[0])

    return dict_df

# Set up combined and summary dataframes
df_all = pd.DataFrame()
df_summary = pd.DataFrame()

# For each image
total_area = 0  # For calculating final density measurement
for imgname in imglist:

    # Read image
    filename = os.path.basename(imgname).split(".")[0]  # Name of individual file
    img = cv2.imread(imgname)  # Image for labeling

    img_size = img.size / 3 * umpix * umpix / 1E6  # Convert area of image (one layer) to mm2
    total_area += img_size  # Record total area of all images for final density calculation

    # Choose correct colors for analysis based on input parameters
    # Main color selection, typically membrane stain
    # OpenCV uses a "BGR" color scheme
    if main_channel is 'r':
        main_img = img[:, :, 2]
        main_color = [0, 0, 255]
    elif main_channel is 'g':
        main_img = img[:, :, 1]
        main_color = [0, 255, 0]
    elif main_channel is 'b':
        main_img = img[:, :, 0]
        main_color = [255, 0, 0]

    # Apply thresholds
    # Here, a binary threshold is applied - consider editing to use adaptive thresholds
    thm, main_img_t = cv2.threshold(main_img, main_threshold, 255, cv2.THRESH_BINARY)
    main_img_t = cv2.morphologyEx(main_img_t, cv2.MORPH_CLOSE, kernel)

    # Calculate size, circularity of each region/cell event
    # Label cell events - detect primary stain "blobs"
    label_img = measure.label(main_img_t)  # Create a labeled image as an input
    main_props = measure.regionprops_table(label_img, properties=('centroid', 'area', 'filled_area', 'image',
                                                                  'convex_area', 'convex_image',
                                                                  'bbox', 'eccentricity', 'coords'))  # Image for count

    # Convert to dataframe, filter
    df_props = pd.DataFrame(main_props)
    if df_props is not None:  # If any cells found
        # Filter by min, max size (two step)
        df_filt = df_props[df_props['area'] > min_area]
        df_filt = df_filt[df_filt['area'] < max_area]

    # Calculate additional properties of cells
    filopodia_count_vector = []  # n filopodia
    min_length = []  # minimum filopodia length
    mean_length = []  # mean length of filopodia
    max_length = []  # maximum filopodia length
    stdev_length = []  # standard deviation of length of filopodia
    texture_vector = []  # Texture, a membrane property

    # Set up images to label
    # Flip layer orientation of original image from BGR to RGB
    img_tolabel = np.dstack((img[:, :, 2], img[:, :, 1], img[:, :, 0]))

    # Create thresholded image
    t_tolabel = np.zeros((img.shape[0], img.shape[1], 3))  # Base - rows, columns, 3 layers
    t_tolabel[np.where(main_img_t == 255)] = main_color  # Primary color


    for i in range(len(df_filt)):
        indices = np.array(df_filt['coords'].iloc[i]).astype(int)  # Pixels comprising single cell

        # Texture
        texture = np.std(main_img[indices[:, 0], indices[:, 1]])  # Standard deviation of original signal
        texture_vector.append(texture)

        # Filopodia count
        # Convex area is used so that inner corners don't also get counted - just outermost points
        # # This could result in some points within the convex shape being missed
        convex = df_filt['convex_image'].iloc[i]
        convex_image = np.asarray(convex * 255).astype(np.uint8)  # Convert to uint8 image for openCV
        # Border allows outermost points to be counted - corner detection doesn't work on points at edge
        image_with_border = cv2.copyMakeBorder(convex_image, top, bottom, left, right, cv2.BORDER_CONSTANT, value=0)

        # Find coordinates of corners
        coords = corner_peaks(corner_harris(image_with_border, k=k), threshold_rel=tr, min_distance=min_distance)

        if coords.any():
            filopodia_count_vector.append(len(coords))
            cell_center = [image_with_border.shape[1]/2, image_with_border.shape[0]/2]  # Find center of cell
            distances = []  # Distance of coordinates from center
            for pt in coords:  # This probably doesn't need to be a loop, would appreciate github pull requests
                # Would also appreciate github pull requests for saving distances as a list within pandas dataframe
                distances.append(math.sqrt((cell_center[0] - pt[0]) ** 2
                                           + (cell_center[1] - pt[1]) ** 2) * umpix)

                # Label filopodia on original and threshold image
                pt1 = int(df_filt['centroid-1'].iloc[i] - convex_image.shape[1] / 2 + pt[1] - 10)
                pt0 = int(df_filt['centroid-0'].iloc[i] - convex_image.shape[0] / 2 + pt[0] - 10)

                cv2.circle(img_tolabel, tuple((pt1, pt0)), 1, (0, 255, 0), 2)
                cv2.circle(t_tolabel, tuple((pt1, pt0)), 1, (0, 255, 0), 2)


            min_length.append(np.min(distances))  # minimum filopodia length
            mean_length.append(np.mean(distances))  # mean length of filopodia
            max_length.append(np.max(distances))  # maximum filopodia length
            stdev_length.append(np.std(distances))  # standard deviation of length of filopodia

        # If not, append 0 to indicate no signal or N/A
        else:
            filopodia_count_vector.append(0)
            min_length.append(0)
            mean_length.append(0)
            max_length.append(0)
            stdev_length.append(0)

    # Append vectors to dataframe as column
    df_filt['Texture (a.u.)'] = texture_vector
    df_filt['Filopodia (n)'] = filopodia_count_vector
    df_filt['Min. filopodia length (\u03bcm)'] = min_length
    df_filt['Mean filopodia length (\u03bcm)'] = mean_length
    df_filt['Max. filopodia length (\u03bcm)'] = max_length
    df_filt['Stdev. filopodia length (\u03bcm)'] = stdev_length

    # Calculated values: area (um)
    df_filt[u'Area (\u03bcm\u00b2)'] = df_filt['area'] * umpix * umpix

    # Add index to resultant dataframe
    index = range(len(df_filt))
    df_filt.insert(0, 'Index', index)

    # Rename additional columns to be saved
    df_filt = df_filt.rename(columns={'centroid-0': 'y', 'centroid-1': 'x',
                                      'area': 'Area (pix)', 'eccentricity': 'Circularity (a.u.)'})

    # Select and reorder columns
    df = df_filt[['Index', 'x', 'y', 'Area (pix)', u'Area (\u03bcm\u00b2)',
                           'Circularity (a.u.)', 'Texture (a.u.)', 'Filopodia (n)', 'Min. filopodia length (\u03bcm)',
                           'Mean filopodia length (\u03bcm)', 'Max. filopodia length (\u03bcm)',
                           'Stdev. filopodia length (\u03bcm)']]

    # Label images - original and threshold (indices on each correspond to the same cell)

    # Write ID text on saved image (cyan)
    for k in range(len(df)):
        # Original image
        cv2.putText(
            img_tolabel,
            str(df['Index'].iloc[k]),
            (int(df['x'].iloc[k]), int(df['y'].iloc[k])),
            cvfont,
            fontScale=0.3,
            color=(255, 0, 255),
            thickness=1)

        # Threshold
        cv2.putText(
            t_tolabel,
            str(df['Index'].iloc[k]),
            (int(df['x'].iloc[k]), int(df['y'].iloc[k])),
            cvfont,
            fontScale=0.3,
            color=(255, 0, 255),
            thickness=1)

    # Save images
    img_tolabel_rgb = cv2.cvtColor(img_tolabel, cv2.COLOR_BGR2RGB)  # Convert color order
    cv2.imwrite(filename + '_labeled-original.png', img_tolabel_rgb)
    cv2.imwrite(filename + '_labeled-threshold.png', t_tolabel)
    cv2.destroyAllWindows()

    # Create pairplots
    df_subset = df[[u'Area (\u03bcm\u00b2)',
                           'Circularity (a.u.)', 'Texture (a.u.)', 'Filopodia (n)', 'Min. filopodia length (\u03bcm)',
                           'Mean filopodia length (\u03bcm)', 'Max. filopodia length (\u03bcm)',
                           'Stdev. filopodia length (\u03bcm)']]
    sns.pairplot(df_subset)
    plt.savefig(filename + '_pairplot.png', dpi=300)
    plt.close()

    # Save individual image dataframe to .xlsx file
    df.to_excel(writer, sheet_name=filename[:30], index=False)  # Filename cropped to prevent excel error

    # Add descriptive statistics for image to summary sheet
    df_image = descriptive_statistics(df, img_size)
    df_image.insert(0, 'Image', filename)
    df_summary = df_summary.append(df_image, ignore_index=True)

    # Add individual image dataframe to df_all
    df.insert(0, 'Image', filename)
    df_all = df_all.append(df, ignore_index=True)

# Create pairplots (with and without functional stain intensity data)
# One color
df_all_subset = df_all[['Image', u'Area (\u03bcm\u00b2)',
                           'Circularity (a.u.)', 'Texture (a.u.)', 'Filopodia (n)', 'Min. filopodia length (\u03bcm)',
                           'Mean filopodia length (\u03bcm)', 'Max. filopodia length (\u03bcm)',
                           'Stdev. filopodia length (\u03bcm)']]
sns.pairplot(df_all_subset)
plt.savefig(dir_name + '_pairplot.png', dpi=300)
plt.close()


# One color per image
sns.pairplot(df_all_subset, hue='Image')
plt.savefig(dir_name + '_multicolor_pairplot.png', dpi=300)
plt.close()

# After all images have been analyzed, write additional data sheets
df_all.to_excel(writer, sheet_name='All data points', index=False)# All data points

# Update summary sheet with summary of all images
dict_df_final = descriptive_statistics(df_all, total_area)
dict_df_final.insert(0, 'Image', 'All images')
df_summary = df_summary.append(dict_df_final, ignore_index=True)
df_summary.to_excel(writer, sheet_name='Descriptive statistics', index=False)

# Save parameters
param_df = pd.DataFrame({u'Ratio, \u03bcm-to-pixels': umpix,
                         'Main channel': main_channel,
                         'Main threshold': main_threshold,
                         'Minimum area': min_area,
                         'Maximum area': max_area,
                         'Minimum distance between filopodia': min_distance,
                         'Sharpness (k)': k,
                         'Analysis date': now.strftime("%D"),
                         'Analysis time': now.strftime("%H:%M:%S")}, index=[1])
param_df.to_excel(writer, sheet_name='Parameters used', index=False)

# Finish code
writer.save()
writer.close()