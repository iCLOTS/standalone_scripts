"""iCLOTS is a free software created for the analysis of common hematology workflow image data

Author: Meredith Fay, Lam Lab, Georgia Institute of Technology and Emory University
Last updated: 2022-07-13
This script corresponds to tools available in version 1.0b1, more recent implementations of tools
may be available within the iCLOTS software and in source code at github.com/iCLOTS

Script that analyzes static, fluorescence microscopy images of cells (platelets, RBCs, or WBCs)
adhered to some surface

Script relies heavily on skimage library region property analysis, tutorial available at:
https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_regionprops.html

Input files:
--Script is designed to work with a folder of image files (.jpg, .png, and/or .tif)
----The same input parameters are applied to each image
Input parameters:
--umpix: The ratio of microns (1e-6 m) to pixels for the image
----Use umpix = 1 for no conversion
--main_channel: 'r' for red/layer 1, 'g' for green/layer 2, 'b' for blue/layer 3
----Value must match 'r', 'g', or 'b' string format exactly
----Typically represents area of cell, like a membrane stain
--fn_channel: like main_channel, but this is a "functional" stain
----Typically represents some activity or characteristic
----Optional additional value 'n' for no functional stain, all functional stain measurements will be 0
--main_threshold: integer between 0 (black), 255 (white) to be used for main channel threshold
----Any value below this threshold becomes 0, or background
----Any value greater than or equal to this threshold becomes 1, or signal to further quantify
--fn_threshold: like main_threshold, but for the function/characteristic stain
--min_area: Minimum area (pixels) of a region (ideally, cell) to be quantified
----Can be used to filter out obvious noise
--max_area: Maximum area (pixels) of a region to be quantified
----Can be used to filter out debris or cell clusters
--min_distance: Minimum distance between maxima of regions of functional stain
----Used with a "peak finding" algorithm


Output files
--Each original image with each detected cell labeled with an index
--Each original image with threshold applied, each detected cell is labeled with the same index

--A corresponding .xlsx sheet containing:
----Area, circularity, texture
----Presence of secondary functional stain, total fluorescence intensity, regions of fluoresence intensity
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
----Each type exported with and without functional stain metrics

Some tips from the iCLOTS team:

Computational and experimental methods:
--For all fluorescent applications, each stain to quantify must be solely in one RGB channel, no other colors
----See export options on your microscopy acquisition software
--After application of the thresholds, skimage analyzes each interconnected region of 1/signal as a cell
--Analysis methods cannot distinguish between overlapping cells
----If cells are significantly overlapping, repeat experiment with a lower cell concentration
--Red blood cells can be difficult to stain
----Antibody staining signal is weak and we've found membrane stains can affect mechanical properties
----Consider using our brightfield adhesion application if this does not conflict with your experimental goals
--Searching for number of regions of functional stain is computationally expensive
----If you don't need this functionality and your dataset is large, consider commenting out related code

Choosing parameters:
--Be sure to use microns-to-pixel ratio, not pixel-to-micron ratio
--Sometimes cells (especially platelets) have a high-intensity "body" and low-intensity spreading or protrusions
----Choose a high threshold if you're primarily quantifying the number of cells
----Choose a low threshold if you're primarily quantifying the morphology of cells
--Depending on the quality of your image, you may want to use an adaptive threshold rather than a binary threshold
----Line to edit indicated within script
----See OpenCV tutorial: https://docs.opencv.org/4.x/d7/d4d/tutorial_py_thresholding.html
--Err on the high side of max_area and low side of min_area parameters
---unless data is particularly noisy or there's a large amount of debris
----By sorting in excel, you can verify that small/large objects are/are not cells

Output files:
--Analysis files are named after the folder containing all images (.xlsx) or image names (.png)
----Avoid spaces, punctuation, etc. within file names
--.xlsx and pairplot data includes a sheet/graph with all images combined
----Only use this when analyzing replicates of the same sample
--Secondary stain metrics are reported in three ways:
----Signal (binary): 0 = no, 1 = yes binary value for presence of secondary stain
------Useful for calculating a percent expression
----Fn. stain intensity (a.u.): summed value of all functional stain pixels within the membrane stain shape
------Take care interpreting, as range of intensity can vary image-to-image due to changes in laser power, etc.
----Fn. stain regions (n): "Peaks" of signal, could represent nuclei lobes, RNA spot, etc.
----No intensity metrics are reported from the main color, as this color should indicate morphology only

"""

# Import statements
#    File management
from tkinter import filedialog
import os
import glob
import datetime
#   Computer vision, image analysis
import cv2  # Computer vision/image processing
from skimage import measure, img_as_float  # Region analysis
from skimage.feature import peak_local_max
#   Number, file, and image management
import numpy as np
import pandas as pd  # For dataframe management
#   Plotting results
import matplotlib.pyplot as plt
import seaborn as sns  # Pairplots
cvfont = cv2.FONT_HERSHEY_SIMPLEX  # Not an import, but used for labeling images


# IMPORTANT: PARAMETERS TO EDIT
umpix = 0.166  # 1 = no conversion
main_channel = 'r'  # 'r' for red/layer 1, 'g' for green/layer 2, 'b' for blue/layer 3
fn_channel = 'g'  # " "
main_threshold = 50  # Integer between 0 (black), 255 (white)
fn_threshold = 30  # " "
min_area = 10
max_area = 2000
min_distance = 5  # Minimum distance between "peaks" in signal, for counting features within cell

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
                  'Fn. stain positive cells (%)': df_input['Signal (binary)'].sum()/len(df_input) * 100,
                  'Min. fn. stain intensity (a.u.)': df_input['Fn. stain intensity (a.u.)'].min(),
                  'Mean fn. stain intensity (a.u.)': df_input['Fn. stain intensity (a.u.)'].mean(),
                  'Max. fn. stain intensity (a.u.)': df_input['Fn. stain intensity (a.u.)'].max(),
                  'Stdev, fn. stain intensity (a.u.)': df_input['Fn. stain intensity (a.u.)'].std(),
                  'Min. fn. stain regions (n)': df_input['Fn. stain regions (n)'].min(),
                  'Mean fn. stain regions (n)': df_input['Fn. stain regions (n)'].mean(),
                  'Max. fn. stain regions (n)': df_input['Fn. stain regions (n)'].max(),
                  'Stdev, fn. stain regions (n)': df_input['Fn. stain regions (n)'].std(),
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

    # Secondary color selection, typically some stain representing a cell function
    if fn_channel is 'r':
        fn_img = img[:, :, 2]
        fn_color = [255, 0, 0]
    elif fn_channel is 'g':
        fn_img = img[:, :, 1]
        fn_color = [0, 255, 0]
    elif fn_channel is 'b':
        fn_img = img[:, :, 0]
        fn_color = [0, 0, 255]
    elif fn_channel is 'n':
        fn_img = np.zeros(img[:, :, 0].shape)
        fn_color = [0, 0, 0]

    # Apply thresholds
    # Here, a binary threshold is applied - consider editing to use adaptive thresholds
    thm, main_img_t = cv2.threshold(main_img, main_threshold, 255, cv2.THRESH_BINARY)
    thf, fn_img_t = cv2.threshold(fn_img, fn_threshold, 255, cv2.THRESH_BINARY)

    # Calculate size, circularity of each region/cell event
    # Label cell events - detect primary stain "blobs"
    label_img = measure.label(main_img_t)  # Create a labeled image as an input
    main_props = measure.regionprops_table(label_img, properties=('centroid', 'area',
                                                                  'bbox', 'eccentricity', 'coords'))

    # Convert to dataframe, filter
    df_props = pd.DataFrame(main_props)
    if df_props is not None:  # If any cells found
        # Filter by min, max size (two step)
        df_filt = df_props[df_props['area'] > min_area]
        df_filt = df_filt[df_filt['area'] < max_area]

    # Calculate additional properties of cells, including functional stain measurements
    binary_int_vector = []  # 0 (no stain) or 1 (stain present)
    int_vector = []  # Summed raw intensity of stain within cell
    n_int_vector = []  # Regions of stain within cell
    texture_vector = []  # Texture, a membrane property

    for i in range(len(df_filt)):
        indices = np.array(df_filt['coords'].iloc[i]).astype(int)  # Pixels comprising single cell

        texture = np.std(main_img[indices[:, 0], indices[:, 1]])  # Standard deviation of original signal
        texture_vector.append(texture)

        # If any functional stain present
        if np.any(fn_img_t[indices[:, 0], indices[:, 1]] > 0):
            binary_int_vector.append(1)  # 1: colocalization yes/no
            int_vector.append(np.sum(fn_img[indices[:, 0], indices[:, 1]]))  # Intensity of signal within cell

            # Count regions of stain
            fn_img_subset = fn_img[df_filt['bbox-0'].iloc[i]:df_filt['bbox-2'].iloc[i],
                            df_filt['bbox-1'].iloc[i]:df_filt['bbox-3'].iloc[i]]
            subset_float = img_as_float(fn_img_subset)
            subset_float_blur = cv2.blur(subset_float, (3,3))

            max_coordinates = peak_local_max(subset_float_blur, min_distance=min_distance)
            n_int_vector.append(len(max_coordinates))  # Save count

        # If not, append 0 to indicate no signal or N/A
        else:
            binary_int_vector.append(0)
            int_vector.append(0)
            n_int_vector.append(0)

    # Append vectors to dataframe as column
    df_filt['Signal (binary)'] = binary_int_vector
    df_filt['Fn. stain intensity (a.u.)'] = int_vector
    df_filt['Fn. stain regions (n)'] = n_int_vector
    df_filt['Texture (a.u.)'] = texture_vector

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
                           'Circularity (a.u.)', 'Texture (a.u.)', 'Signal (binary)',
                           'Fn. stain intensity (a.u.)', 'Fn. stain regions (n)']]


    # Label images - original and threshold (indices on each correspond to the same cell)

    # Flip layer orientation of original image from BGR to RGB
    img_tolabel = np.dstack((img[:, :, 2], img[:, :, 1], img[:, :, 0]))

    # Create thresholded image
    t_tolabel = np.zeros((img.shape[0], img.shape[1], 3))  # Base - rows, columns, 3 layers
    t_tolabel[np.where(main_img_t == 255)] = main_color  # Primary color
    t_tolabel[np.where(fn_img_t == 255)] = fn_color  # Secondary stain

    # Write ID text on saved image (cyan)
    for j in range(len(df)):
        # Original image
        cv2.putText(
            img_tolabel,
            str(df['Index'].iloc[j]),
            (int(df['x'].iloc[j]), int(df['y'].iloc[j])),
            cvfont,
            fontScale=0.3,
            color=(255, 0, 255),
            thickness=1)

        # Threshold
        cv2.putText(
            t_tolabel,
            str(df['Index'].iloc[j]),
            (int(df['x'].iloc[j]), int(df['y'].iloc[j])),
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
    df_subset = df[[u'Area (\u03bcm\u00b2)', 'Circularity (a.u.)', 'Texture (a.u.)',
                            'Signal (binary)', 'Fn. stain intensity (a.u.)', 'Fn. stain regions (n)']]
    sns.pairplot(df_subset)
    plt.savefig(filename + '_pairplot.png', dpi=300)
    plt.close()

    df_subset_noint = df[[u'Area (\u03bcm\u00b2)', 'Circularity (a.u.)', 'Texture (a.u.)']]
    sns.pairplot(df_subset_noint)
    plt.savefig(filename + '_pairplot_no-int.png', dpi=300)
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
df_all_subset = df_all[['Image', u'Area (\u03bcm\u00b2)', 'Circularity (a.u.)', 'Texture (a.u.)',
                        'Signal (binary)', 'Fn. stain intensity (a.u.)', 'Fn. stain regions (n)']]
sns.pairplot(df_all_subset)
plt.savefig(dir_name + '_pairplot.png', dpi=300)
plt.close()
df_all_subset_noint = df_all[['Image', u'Area (\u03bcm\u00b2)', 'Circularity (a.u.)', 'Texture (a.u.)']]
sns.pairplot(df_all_subset_noint)
plt.savefig(dir_name + '_pairplot_no-int.png', dpi=300)
plt.close()


# One color per image
sns.pairplot(df_all_subset, hue='Image')
plt.savefig(dir_name + '_multicolor_pairplot.png', dpi=300)
plt.close()
sns.pairplot(df_all_subset_noint, hue='Image')
plt.savefig(dir_name + '_multicolor_pairplot_no-int.png', dpi=300)
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
                         'Functional channel': fn_channel,
                         'Functional threshold': fn_threshold,
                         'Minimum area': min_area,
                         'Maximum area': max_area,
                         'Minimum peak distance': min_distance,
                         'Analysis date': now.strftime("%D"),
                         'Analysis time': now.strftime("%H:%M:%S")}, index=[1])
param_df.to_excel(writer, sheet_name='Parameters used', index=False)

# Finish code
writer.save()
writer.close()