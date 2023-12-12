"""iCLOTS is a free software created for the analysis of common hematology workflow image data

Author: Meredith Fay, Lam Lab, Georgia Institute of Technology and Emory University
Last updated: 2022-07-12
This script corresponds to tools available in version 1.0b1, more recent implementations of tools
may be available within the iCLOTS software and in source code at github.com/iCLOTS

Script that analyzes static, brightfield images of cells (platelets, RBCs, or WBCs)
adhered to some surface
May also be suitable for preliminary digital pathology approaches

Script relies heavily on Trackpy python library, documentation available at:
http://soft-matter.github.io/trackpy/v0.5.0/

Input files:
--Script is designed to work with a folder of image files (.jpg, .png, and/or .tif)
----The same input parameters are applied to each image
Input parameters:
--umpix: The ratio of microns (1e-6 m) to pixels for the image
----Use umpix = 1 for no conversion
--max_diameter: Maximum diameter of a cell to be considered, MUST be an odd integer
--min_mass: Minimum summed intensity of the cell
--invert: Binary value
----True: searching for dark features on a light background (most appropriate for brightfield)
----False: searching for light features on a dark background (may be appropriate for WBCs)

Output files
--Each original image with each detected cell labeled with an index

--A corresponding .xlsx sheet containing:
----Mass, area, and circularity for each cell - one sheet/image
----Mass, area, and circularity for all cells - one sheet (this combined sheet is best for analyzing replicates)
----Descriptive statistics (min, mean, max, standard deviation for area and circularity)
------For each image and for all images combined
------Descriptive statistics sheet also includes a density measurement (n cells/mm^2)
----Parameters used and time/date analysis was performed, for reference

--Pairplot including area and circularity metrics
----For each image
----For all images, combined, one color represents all
----For all images, combined, each color represents a different image

Some tips from the iCLOTS team:

Computational and experimental methods:
--Trackpy searches for particles represented by image regions with Gaussian-like distributions
---of pixel brightness
----Cells with an especially heterogenous appearance, such as bi-concave RBCs or actived WBCs, may
-----be hard to detect. While count will be accurate, take care in interpreting area or circularity
--Analysis methods cannot distinguish between overlapping cells
----If cells are significantly overlapping, repeat experiment with a lower cell concentration
--Owing to the heterogenous appearance of blood cell types, brightfield analysis is always challenging
----Consider using a fluorescent membrane stain coupled with our fluorescent adhesion applications
-----if this does not conflict with your experimental goals, especially for WBCs/neutrophils

Choosing parameters:
--Be sure to use microns-to-pixel ratio, not pixel-to-micron ratio
--Err on the high side of max_diameter and low side of min_mass parameters
---unless data is particularly noisy or there's a large amount of debris
----By sorting in excel, you can verify that small/large objects are/are not cells

Output files:
--Analysis files are named after the folder containing all images (.xlsx) or image names (.png)
----Avoid spaces, punctuation, etc. within file names
--.xlsx and pairplot data includes a sheet/graph with all images combined
----Only use this when analyzing replicates of the same sample

"""

# Import statements
#    File management
from tkinter import filedialog
import os
import glob
import datetime
#   Computer vision, image analysis
import cv2  # Computer vision/image processing
import trackpy as tp  # Particle tracking
#   Number, file, and image management
import pandas as pd  # For dataframe management
#   Plotting results
import matplotlib.pyplot as plt
import seaborn as sns  # Pairplots
from PIL import Image, ImageDraw  # For labeling images
#   Misc.
from math import pi

# IMPORTANT: PARAMETERS TO EDIT
umpix = 1  # 1 = no conversion
max_diameter = 23  # Err on high side, MUST be an odd integer
min_mass = 2200  # Err on low side
invert=False  # True is dark on light, False is light on dark

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
                  'Field of view (mm\u00b2)': img_size,
                  'Cell density (n/mm\u00b2)': len(df_input)/img_size
                  }
    dict_df = pd.DataFrame(dict, index=[0])

    return dict_df

# For each image
total_area = 0  # For calculating final density measurement
for imgname in imglist:

    # Read image
    filename = os.path.basename(imgname).split(".")[0]  # Name of individual file
    img_color = cv2.imread(imgname)  # Image for labeling
    img = cv2.imread(imgname, 0)  # Image for analysis, 0 flag denotes grayscale
    img_size = img.size * umpix * umpix / 1E6  # Convert area of image to mm2
    total_area += img_size  # Record total area of all images for final density calculation

    # Locate particles (ideally, cells) using Trackpy
    # See walkthrough: http://soft-matter.github.io/trackpy/dev/tutorial/walkthrough.html
    f = tp.locate(img, max_diameter, minmass=min_mass, invert=invert)

    # Add index to resultant dataframe
    index = range(len(f))
    f.insert(0, 'Index', index)

    # Take most useful subset
    f = f[['Index', 'x', 'y', 'mass', 'size', 'ecc']]

    # Calculate additional metrics
    f['Area (pix)'] = f['size'] * f['size'] * pi  # pi*r^2
    f[u'Area (\u03bcm\u00b2)'] = f['Area (pix)'] * umpix * umpix

    # Rename columns
    f = f.rename(columns={'size': 'Radius (pix)',
                          'ecc': 'Circularity (a.u.)',
                          'mass': 'Mass (a.u.)'
                          })

    # Label image
    PILimg = Image.fromarray(img_color)  # Pillow library prepares image for drawing
    drawimg = ImageDraw.Draw(PILimg)
    for i in range(len(f)):
        drawimg.text((f['x'].iloc[i], f['y'].iloc[i]), str(i), fill="#ff0000")
    PILimg.save(filename + '_labeled.png')

    # Create pairplot
    f_subset = f[[u'Area (\u03bcm\u00b2)', 'Circularity (a.u.)']]
    sns.pairplot(f_subset)
    plt.savefig(filename + '_pairplot.png', dpi=300)
    plt.close()

    # Save individual image dataframe to .xlsx file
    f.to_excel(writer, sheet_name=filename[:30], index=False)  # Filename cropped to prevent excel error

    # Add descriptive statistics for image to summary sheet
    df_image = descriptive_statistics(f, img_size)
    df_image.insert(0, 'Image', filename)
    df_summary = df_summary.append(df_image, ignore_index=True)

    # Add individual image dataframe to df_all
    f.insert(0, 'Image', filename)
    df_all = df_all.append(f, ignore_index=True)

# Create pairplots
# One color
df_all_subset = df_all[['Image', u'Area (\u03bcm\u00b2)', 'Circularity (a.u.)']]
sns.pairplot(df_all_subset)
plt.savefig(dir_name + '_singlecolor_pairplot.png', dpi=300)
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
                         'Maximum cell diameter': max_diameter,
                         'Minimum mass': min_mass,
                         'Invert': str(invert),
                         'Analysis date': now.strftime("%D"),
                         'Analysis time': now.strftime("%H:%M:%S")}, index=[1])
param_df.to_excel(writer, sheet_name='Parameters used', index=False)

# Finish code
writer.save()
writer.close()