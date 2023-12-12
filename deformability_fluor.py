"""iCLOTS is a free software created for the analysis of common hematology workflow image data
Author: Meredith Fay, Lam Lab, Georgia Institute of Technology and Emory University
Last updated: 2022-07-28

This script corresponds to tools available in version 1.0b1, more recent implementations of tools
may be available within the iCLOTS software and in source code at github.com/iCLOTS

Script that analyzes fluorescent videomicroscpy of cells transiting
the biophysical flow cytometer device, initial description of assay available at:
https://pubmed.ncbi.nlm.nih.gov/18584080/

Functionality as-is is not available within iCLOTS.
A more generalized version of this script exists as sct_fluor (sct indicates "single cell tracking")
This script tracks single cells through any microfluidic device
Single cell tracking applications are available in iCLOTS version 1.0b1


Script relies heavily on Trackpy python library, documentation available at:
http://soft-matter.github.io/trackpy/v0.5.0/

Input files:
--Script is designed to work a directory of videomicroscopy files (.avi)
----If your data is saved as a series of frames, please see the suite of video editing tools to convert to .avi

Input parameters:
--umpix: The ratio of microns (1e-6 m) to pixels for the videos
----Use umpix = 1 for no conversion
--fps: Frames per second, the rate of imaging
----Note that FPS values pulled directly from videos can be inaccurate, especially if the video
-----has been resized or edited in any way
--max_diameter: The maximum diameter of the tracked cells, MUST be an odd integer
--labelimg: Boolean variable (True or False) indicating if you'd like labeled image data

--min_mass for trackpy algorithms is set artificially low such that all cells are counted, even if dim

Output files
--If labelimg is true, each original frame of the provided video with each detected cell labeled with an index
--A corresponding .xlsx sheet containing:
----sDI, area, summed cell fluorescence intensity
------For each cell - one sheet/video
------For all cells - one sheet (this combined sheet is best for analyzing replicates)
----Descriptive statistics (min, mean, max, standard deviation for all metrics)
------For each video and for all videos combined
----Parameters used and time/date analysis was performed, for reference

--Pairplot
----Three types:
------For each video
------For all videos, combined, one color represents all
------For all videos, combined, each color represents a different video

Some tips from the iCLOTS team:

Computational and experimental methods:
--Choose the height of your biophysical flow cytometer mask carefully
----Cells must necessarily deform, but they cannot obstruct
------Try ~2 um tall for RBCs, ~12-14 um tall for WBCs (channel width ~6 um)
----WBCs can be challenging as they are "stickier"
------Be sure to coat the channel with a bovine serum albumin solution to prevent non-specific binding
----For best results, do not reuse devices for multiple experimental runs
--Trackpy searches for particles represented by image regions with Gaussian-like distributions
---of pixel brightness
--Analysis methods cannot distinguish between cells transiting the channel together, in contact
----If cells are significantly contacting, repeat experiment with a lower cell concentration
--It is very hard to get a bright signal from a cell moving through the BFC device
----Should cells be very dim, use the edit contrast tool
------Set alpha as 255 / the maximum intensity of the video
------After analysis, divide by alpha to get the original raw intensity value
--------Alpha will bias values, making bright values appear brighter
--The "map" to choose an ROI presented is based on all signal present in the video
----Building this map is computationally expensive
----Do your best choosing an ROI
----If the signal was too weak to form a good map, try again with a higher cell concentration
--If the analysis is taking an unacceptably long time, you can resize videos to be smaller
----This may cause you to miss the signal and/pr smallest cells

Choosing parameters:
--Be sure to use microns-to-pixel ratio, not pixel-to-micron ratio
--Err on the high side of max_diameter and low side of min_mass parameters
---unless data is particularly noisy or there's a large amount of debris

Output files:
--Analysis files are named after the folder containing all videos (.xlsx) or video names (.png)
----Avoid spaces, punctuation, etc. within file names
--Library used to write excel sheets crops filenames for returned data sheets to avoid corrupted files
----Make sure the first 15 characters of each video name are distinct
--.xlsx and pairplot data includes a sheet/graph with all videos combined
----Only use this when analyzing replicates of the same sample

Some quality metrics are set with the biophysical flow cytometer in mind
--Minimum number of frames (3)
--Minimum distance traveled (1/3 the width of the rOI)
--Search range (1/3 the width of the ROI)
--Distance is calculated as x-direction only
You may want to edit these parameters for your specific usage

"""

# Import libraries
# File management
from tkinter import filedialog
import os
import glob
import shutil
import datetime
# Number, file, and image management
import cv2  # Computer vision/image processing
import numpy as np  # For array management
import pandas as pd  # For database management
import math
# Particle tracking
import trackpy as tp
import warnings
warnings.filterwarnings("ignore", module="trackpy")
# Labeling and plotting results
from PIL import Image, ImageDraw
import matplotlib.pyplot as plt
import seaborn as sns


# IMPORTANT: PARAMETERS TO EDIT
umpix = 1  # 1 = no conversion
fps = 25  # Frames per second rate of imaging
max_diameter = 41  # Maximum diameter of tracked cells
min_mass = 100  # Minimum intensity of a tracked cell
# If you'd like graphical data and the images labeled with the tracked cells, set as "True"
labelimg = True  # Recommended

# Select directory of files
dirpath = filedialog.askdirectory()

# Create a directory for saved results including time at which operation was performed
now = datetime.datetime.now()
# Create strings to indicate operations performed
output_folder = os.path.join(dirpath, 'Analysis, ' + now.strftime("%m_%d_%Y, %H_%M_%S"))
os.mkdir(output_folder)
os.chdir(output_folder)

# Create a list of all video files
video_list = glob.glob(dirpath + "/*.avi")

# Set up method to write final excel files
dir_name = os.path.basename(dirpath)  # Also used for final graphs
excel_name = dir_name + '_analysis.xlsx'
writer = pd.ExcelWriter(excel_name, engine='openpyxl')

# Create background subtractor for analysis
fgbg = cv2.createBackgroundSubtractorMOG2(detectShadows=False)

def descriptive_statistics(df_input):
    """Function to calculate descriptive statistics for each population, represented as a dataframe"""

    dict = {'n cells': len(df_input),
            u'Min. sDI (\u03bcm/s)': df_input['sDI (\u03bcm/s)'].min(),
            u'Mean sDI (\u03bcm/s)': df_input['sDI (\u03bcm/s)'].mean(),
            u'Max. sDI (\u03bcm/s)': df_input['sDI (\u03bcm/s)'].max(),
            u'Stdev, sDI (\u03bcm/s)': df_input['sDI (\u03bcm/s)'].std(),
            u'Min. area (pix)': df_input['Area (pix)'].min(),
            u'Mean area (pix)': df_input['Area (pix)'].mean(),
            u'Max. area (pix)': df_input['Area (pix)'].max(),
            u'Stdev, area (pix)': df_input['Area (pix)'].std(),
            u'Min. fl. int. (a.u.)': df_input['Fl. int. (a.u.)'].min(),
            u'Mean fl. int. (a.u.)': df_input['Fl. int. (a.u.)'].mean(),
            u'Max. fl. int. (a.u.)': df_input['Fl. int. (a.u.)'].max(),
            u'Stdev, fl. int. (a.u.)': df_input['Fl. int. (a.u.)'].std()
            }
    dict_df = pd.DataFrame(dict, index=[0])

    return dict_df

# Set up combined and summary dataframes
df_all = pd.DataFrame()
df_summary = pd.DataFrame()

# For each video
for video in video_list:

    filename = os.path.basename(video).split(".")[0]  # Name of individual file
    os.chdir(output_folder)  # Return to original analysis folder

    # Read video into an array using openCV
    cap = cv2.VideoCapture(video)  # Capture video
    n_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))  # Record number of frames within video
    w = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))  # Record width of video frame (pixels)
    h = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))  # Record height of video frame (pixels)
    # Create an empty array to save frames with background removed
    # Helps with removing channel walls, ignoring debris
    frarray = np.empty((h, w, n_frames), np.dtype('uint8'))
    gray_frarray = np.empty((h, w, n_frames), np.dtype('uint8'))

    fgbg = cv2.createBackgroundSubtractorMOG2(detectShadows=False)  # Create background subtractor
    kernel = np.ones((15, 15), np.uint8)
    fc = 0  # Initialize frame count for video processing loop
    ret = True  # Initialize "true" value for video processing loop
    while fc < n_frames and ret:
        ret, frame = cap.read()  # Read the video capture object
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)  # Read as gray/one layer
        mask = fgbg.apply(gray)  # Apply background subtraction
        closing = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
        frarray[:, :, fc] = closing  # Save one layer image to array
        gray_frarray[:, :, fc] = gray  # Save original image
        fc += 1  # Update count after frame is processed

    # Generate map of all signal to assist with fluorescent ROI selection
    map = np.sum(gray_frarray, axis=2)  # Sum along layer axis
    map_scaled = map * 255 // map.max()  # Scale max. value to 255
    map_scaled_img = map_scaled.astype(np.uint8)  # Convert to uint8

    # Apply a threshold
    threshold = np.percentile(map_scaled_img, 70)  # Create a threshold by selecting signal above the 70th percentile
    d_c_ret, img_for_roi = cv2.threshold(map_scaled_img, threshold, 255,
                                    cv2.THRESH_BINARY)  # Apply binary threshold to create final map
    # Something like an Otsu threshold won't work because most of the signal is very light

    # Choose ROI from last frame (often initial frames have changes in illumination)
    fromCenter = False  # Set up to choose as a drag-able rectangle
    r = cv2.selectROI("Image", img_for_roi, fromCenter)  # Choose ROI
    ROI_x = int(r[0])  # Take result of selectROI and place into a variable
    ROI_y = int(r[1])  # " "
    ROI_w = int(r[2])  # " "
    ROI_h = int(r[3])  # " "

    ROI_w_um = ROI_w * umpix  # Width converted to microns

    # Create a small kernel for morphological operations
    kernel = np.ones((5, 5), np.uint8)

    crop_img = []  # To save cropped frames
    for i in range(n_frames - 1):
        scaled_img = frarray[:, :, i] / 255  # Convert to 0/1
        frame = gray_frarray[:, :, i].astype(float)  # Original intensity values
        img = frame * scaled_img  # Intensity only in regions where cell is detected

        crop = img[ROI_y: (ROI_y + ROI_h), ROI_x: (ROI_x + ROI_w)]

        crop_img.append(crop.copy())

    # Begin trackpy tracking analysis
    tp.quiet()
    f = tp.batch(crop_img[:n_frames], max_diameter, minmass=1000, invert=False); # Detect particles/cells
    # If f exists, link particles, cells into dataframe format
    if len(f) > 0:
        # Search range criteria: must travel no further than 1/3 the channel length in one frame
        # Memory here signifies a particle/cell cannot "disappear" for more than one frame
        tr = tp.link_df(f, search_range=ROI_w/3, memory=1, adaptive_stop=1, adaptive_step=0.95)
        # Filter stubs criteria requires a particle/cell to be present for at least three frames
        t_final = tp.filter_stubs(tr, 3)

        # If linked cells are found
        if len(t_final) > 0:

            # Series of vectors for final results dataframe
            p_i = []  # Particle index
            f_start = []  # Start frame, frame where cell first detected
            f_end = []  # End frame, frame where cell last detected
            dist = []  # Distance traveled
            time = []  # Time for travel
            sizes = []  # Cell size
            fl_int = []  # Cell intensity
            t_sdi = pd.DataFrame()  # Create dataframe
            # For each particle, calculate RDI and save data for results dataframe:
            for p in range(t_final['particle'].iloc[-1]):
                df_p = tr[tr['particle'] == p]  # Region of trackpy dataframe corresponding to individual particle index
                x_0 = df_p['x'].iloc[0]  # First x-position
                x_n = df_p['x'].iloc[-1]  # Last x-position
                f_0 = df_p['frame'].iloc[0]  # First frame number
                f_n = df_p['frame'].iloc[-1]  # Last frame number
                s = df_p['size'].mean() * df_p['size'].mean() * math.pi  # Area of cell (pi*r^2)
                m = df_p['mass'].mean()  # Intensity
                d = (x_n - x_0) * umpix  # Distance (microns) - x-direction only
                t = (f_n - f_0) / fps  # Time (seconds)
                # Criteria to save cells as a valid data point:
                # Must travel no less than 1/3 the length of channel
                # Must travel no further than length of channel
                if d < ROI_w_um and d > ROI_w_um / 3:
                    t_sdi = t_sdi.append(df_p, ignore_index=True)  # Save trackpy metrics
                    # Append data for particle/cell
                    p_i.append(p)
                    f_start.append(f_0)
                    f_end.append(f_n)
                    dist.append(d)
                    time.append(t)
                    sizes.append(s)  # Background subtractor changes size of cell, size is a relative measurement
                    fl_int.append(m)

            # Calculate sDI by dividing distance by time (um/sec)
            sdi = []
            sdi = np.asarray([u / v for u, v in zip(dist, time)])

            # Organize time, location, and RDI data in a list format
            df_video = pd.DataFrame(
                {'Particle': p_i,
                 'Start frame': f_start,
                 'End frame': f_end,
                 'Transit time (s)': time,
                 'Distance traveled (\u03bcm)': dist,
                 'sDI (\u03bcm/s)': sdi,
                 'Area (pix)': sizes,
                 'Fl. int. (a.u.)': fl_int
                })

            # Renumber particles 0 to n
            df_video['Particle'] = np.arange(len(df_video))
            uniqvals = t_sdi['particle'].unique()
            newvals = np.arange(len(uniqvals))
            for val in newvals:
                uniqval = uniqvals[val]
                t_sdi['particle'] = t_sdi['particle'].replace(uniqval, val)

            # Final data to excel
            # Filename cropped to prevent excel errors caused by a sheet name over 30 characters
            if len(filename) > 14:
                filename = filename[:14]
            df_video.to_excel(writer, sheet_name = filename + ', all', index=False)  # RDI
            t_sdi.to_excel(writer, sheet_name= filename +', trackpy', index=False)  # Trackpy outputs

            # Add descriptive statistics for video to summary sheet
            df_file = descriptive_statistics(df_video)
            df_file.insert(0, 'Video', filename)
            df_summary = df_summary.append(df_file, ignore_index=True)

            # Add individual video dataframe to df_all
            df_video.insert(0, 'Video', filename)
            df_all = df_all.append(df_video, ignore_index=True)

            # Pairplot
            df_subset = df_video[['sDI (\u03bcm/s)', 'Area (pix)', 'Fl. int. (a.u.)']]
            sns.pairplot(df_subset)
            plt.savefig(filename + '_pairplot.png', dpi=300)
            plt.close()

            # If user would like graphical data and frames labeled with particle numbers
            # Computationally expensive but useful, so recommended
            if labelimg is True:
                # Create a directory for that image
                dir_folder = output_folder + '/' + filename
                if os.path.exists(dir_folder):
                    shutil.rmtree(dir_folder)  # Delete if already exists, prevents error
                os.makedirs(dir_folder)
                os.chdir(dir_folder)

                # Read video as an OpenCV video object
                cap = cv2.VideoCapture(video)
                # Label
                success, image = cap.read()
                count = 1
                while success:
                    # New filename: original, frame, frame number
                    image_name = filename + '_frame_' + str(count).zfill(5)
                    success, image = cap.read()
                    if image is not None:
                        f = t_sdi[t_sdi['frame'] == count]
                        # Set up image to label, including cropping
                        PILimg = Image.fromarray(image[ROI_y: (ROI_y + ROI_h), ROI_x: (ROI_x + ROI_w)])
                        drawimg = ImageDraw.Draw(PILimg)  # " "
                        for i in range(len(f)):
                            drawimg.text((f['x'].iloc[i], f['y'].iloc[i]), str(f['particle'].iloc[i]),
                                         fill="#ff0000")  # Label
                        PILimg.save(image_name + "_labeled.png")  # Save image
                    count += 1

os.chdir(output_folder)  # Return to original analysis folder

# Save and create summary excel sheets and pairplots
# Create pairplots (with and without functional stain intensity data)
# One color
df_all_subset = df_all[['Video', 'sDI (\u03bcm/s)', 'Area (pix)', 'Fl. int. (a.u.)']]
sns.pairplot(df_all_subset)
plt.savefig(dir_name + '_pairplot.png', dpi=300)
plt.close()

# One color per video
sns.pairplot(df_all_subset, hue='Video')
plt.savefig(dir_name + '_multicolor_pairplot.png', dpi=300)
plt.close()

# After all videos have been analyzed, write additional data sheets
df_all.to_excel(writer, sheet_name='All data points', index=False)# All data points

# Update summary sheet with summary of all videos
dict_df_final = descriptive_statistics(df_all)
dict_df_final.insert(0, 'Video', 'All videos')
df_summary = df_summary.append(dict_df_final, ignore_index=True)
df_summary.to_excel(writer, sheet_name='Descriptive statistics', index=False)

# Save parameters
param_df = pd.DataFrame({u'Ratio, \u03bcm-to-pixels': umpix,
                         'FPS': fps,
                         'Max. diameter': max_diameter,
                         'Min. intensity': 1000,
                         'Label image': labelimg,
                         'Analysis date': now.strftime("%D"),
                         'Analysis time': now.strftime("%H:%M:%S")}, index=[1])
param_df.to_excel(writer, sheet_name='Parameters used', index=False)

# Save and close excel file writer
writer.save()
writer.close()