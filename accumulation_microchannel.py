"""iCLOTS is a free software created for the analysis of common microfluidic and hematology workflow image data

Author: Meredith Fay, Lam Lab, Georgia Institute of Technology and Emory University
Last updated: 2022-10-29
This script corresponds to tools available in version 1.0b1, more recent implementations of tools
may be available within the iCLOTS software and in source code at github.com/iCLOTS

Microchannel accumulation/occlusion analysis GUI analyzes single or timeseries images
of cells adhered to one or many straight-channel portions of a microfluidic device
Accomodates red, blue, and/or green channels

Application relies heavily on skimage library region property analysis, tutorial available at:
https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_regionprops.html

Input files:
--Application is designed to work with a folder of images (.jpg, .png, and/or .tif)
----The same input parameters are applied to each image
----This code will work on a single image
----It's important that frames are labeled sequentially in timeseries order (i.e. 01, 02.. 10 vs. 1, 2... 10)
----Each image should consist of one or many straight portions of a microfluidic device
----After uploading one or several images, the user is prompted to choose an ROI from the first image
------This ROI should contain the straight channel portions
------The same ROI is applied to all images, take care that all images represent the same field of view
--The algorithm relies on left-to-right indexing to form the channel regions to analyze
----As such, channels should be perfectly horizontal
------iCLOTS provides a video-editing rotation tool that does not affect aspect ratio
--In order to create a complete channel area to analyze, some fluorescence signal must be present
---at every y pixel of the channel
----Staining the channels, or some feature of the channel, like a cell layer, helps with this

Input parameters:
--umpix: micron to pixel ratio: The ratio of microns (1e-6 m) to pixels for the image
----Use value = 1 for no conversion
--Presence of fluorescent signal (rchannel, gchannel, bchannel) - boolean stating whether or not to analyze a channel
--Red, green, and/or blue threshold(s) (rthresh, gthresh, bthresh) - [0, 255] threshold value for each channel

Output files
--Single pixel resolution
--ROI from each image
--The "map" used - the channel region(s) as detected
--Images with image processing (threshold) steps applied
----For each frame, each selected color as detected by the set threshold overlaid on the channel map

--A corresponding .xlsx sheet containing:
----For each selected channel:
------Raw data: A percent y-occlusion for very frame, channel, x-position within the channel
--------Percent y-occlusion indicates what percentage of the height of the microchannel contains signal
------Channel data: Area of signal and signal accumulation (pixels, um2)
-------and %y occlusion for each channel in each frame
------Frame data: mean signal area, signal accumulation, and %y occlusion per frame (all channels)
------Conversion notes:
--------To convert accumulation per frame into per timepoint, divide frame # by FPS imaging rate
--------To convert x-pixel coordinate to a measurement, multiple by \u03bcm-to-pixel ratio)

--Future versions of iCLOTS will include obstruction, the percent of the y-direction for each channel
---containing some signal

--Occlusion/accumulation graph:
----For the time series, a line graph showing:
------Occlusion (titled, left) for each channel (light lines) and mean (dark lines) for each color
------Accumulation (titled, right) - " "

Some tips from the iCLOTS team:

Computational and experimental methods:
--See input requirements: a time series, in the same field of view, with "complete" y-height horizontal channels
----The left-to-right indexing to form the channels requires some signal at every height point in that channel
------Consider staining the microfluidic channels
------We are planning a brightfield/fluoresence microscopy combo. app for future iCLOTS releases
--Time series images must be in the proper alphabetical/numerical order
----If image names contain numbers, use preceding zeros to order properly
------i.e. 01, 02... 10 instead of 1, 2... 10
--The Lam lab has developed these methods on an "endothelialized" branching microfluidic device
----See "Endothelialized Microfluidics for Studying Microvascular Interactions in Hematologic Diseases"
-----manuscript by Myers and Sakurai et al., 2012, JOVE
----We are happy to share a detailed endothelialization protocol upon request
----We are happy to share the mask design files and instructions for fabrication upon request

Choosing parameters:
--Be sure to use micron-to-pixel ratio, not pixel-to-micron ratio
--Depending on the threshold you set, while the "trend" of accumulation/occlusion will stay constant,
---the degree of accumulation/occlusion will decrease as threshold increases
----If you are comparing conditions, make sure they were taken with the same imaging settings
-----and use the same threshold values
------Ideally these experiments are direct control-to-experimental comparisons taken on the same day

Output files:
--Analysis files are named after the folder containing all images (.xlsx) or image names (.png)
----Avoid spaces, punctuation, etc. within file names

"""

# Import statements
#    File management
from tkinter import filedialog
import os
import glob
#   Image processing and numerical data
import cv2
import numpy as np
import pandas as pd
from skimage import measure
from scipy.ndimage import label
import matplotlib.pyplot as plt
import datetime
import shutil

# IMPORTANT: PARAMETERS TO EDIT
umpix = 0.5  # micron-to-pixel ratio
rchannel = True  # Present channels
gchannel = False
bchannel = True
rthresh = 50  # Channel thresholds
gthresh = 50
bthresh = 10

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
filelist = imglist_png + imglist_jpg + imglist_tif

# Set up method to write final excel files
foldername = os.path.basename(dirpath)  # Also used for final graphs
excel_name = foldername + '_analysis.xlsx'
writer = pd.ExcelWriter(excel_name, engine='openpyxl')

# Set up dataframes
df_all = pd.DataFrame()
df_image = pd.DataFrame()
df_raw = pd.DataFrame()
df_data_byframe = pd.DataFrame()
df_data_raw = pd.DataFrame()
df = pd.DataFrame()  # For raw data
df_summary_all = pd.DataFrame()  # For summary data
df_summary_frame = pd.DataFrame()  # For summary data, summarized into frames

df_colors = pd.DataFrame()
df_img = pd.DataFrame()

def analysis_math(df_img, mapbin_ext, filelist, umpix, layer, threshold, x, y, w, h):
    """Calculates occlusion and accumulation for desired RGB channel"""

    # Create a numbered list of channels
    lbl, nlbls = label(mapbin_ext)
    rp = measure.regionprops_table(lbl, properties=('label', 'bbox', 'centroid'))
    rp_df = pd.DataFrame(rp)
    lab = rp_df['label']
    x1 = rp_df['bbox-0']
    x2 = rp_df['bbox-2']

    data_raw = []
    data = []
    a0 = np.zeros(len(lab))  # List of channels zero long
    # Column index titles
    # rng = range(self.w.get())
    rng = range(w)
    ypix_names = [str(a) for a in rng]
    col_names = ['Frame', 'Channel'] + ypix_names

    df_img = df_img.append({'name': 'map', 'color': 'map', 'img': mapbin_ext}, ignore_index=True)

    for i in range(len(filelist)):
        img = cv2.imread(filelist[i])
        imgbasename = os.path.basename(filelist[i].split(".")[0])
        crop = img[y:(y + h), x:(x + w), :]  # Create cropped image

        df_img = df_img.append({'name': imgbasename, 'color': 'full', 'img': crop}, ignore_index=True)

        img_channel = crop[:, :, layer]

        # Add image to array
        ret, img_thresh = cv2.threshold(img_channel, threshold, 255, cv2.THRESH_BINARY)  # Red

        # Layer
        map_save = np.dstack((mapbin_ext, mapbin_ext, mapbin_ext))
        color = [0, 0, 0]
        color[layer] = 255
        map_save[np.where(img_thresh == 255)] = color

        df_img = df_img.append({'name': imgbasename, 'color': str(layer), 'img': map_save}, ignore_index=True)

        for j in range(len(lab)):  # For each channel
            # Index channel out of threshold image
            channel = np.array(img_thresh[x1[j]:x2[j]][:])/255
            # Sum along one axis and divide by height for a percent
            occ_vector = np.sum(channel, axis=0)/channel.shape[0] * 100

            # Mean percent occlusion across channel
            occ_mean = np.mean(occ_vector)
            # Max occlusion across channel
            occ_max = np.max(occ_vector)
            # Total area of signal
            area_channel = np.sum(channel)
            # Accumulation from previous frame
            acc_channel = area_channel - a0[j]  # Subtract previous area

            # Convert numbers to microns
            area_um = area_channel * umpix * umpix
            acc_um = acc_channel * umpix * umpix
            data_raw.append([i, j] + occ_vector.tolist())
            data.append([i, j, occ_mean, occ_max, area_channel, acc_channel, area_um, acc_um])

            # Reset
            a0[j] = area_channel

    # Save and return as dataframes
    df_data_raw = pd.DataFrame(data_raw, columns=col_names)
    df_data = pd.DataFrame(data, columns=['Frame', 'Channel', 'Mean occlusion (%)',
                                          'Max. occlusion (%)', 'Area (pix)', 'Accumulation (pix)',
                                          u'Area (\u03bcm\u00b2)', u'Accumulation (\u03bcm\u00b2)'])
    df_data_byframe = pd.DataFrame(df_data.groupby(['Frame']).mean(), columns=['Mean occlusion (%)',
                                                                               'Max. occlusion (%)',
                                                                               'Area (pix)',
                                                                               'Accumulation (pix)',
                                                                               u'Area (\u03bcm\u00b2)',
                                                                               u'Accumulation (\u03bcm\u00b2)'])
    framelist = np.linspace(0, len(df_data_byframe)-1, len(df_data_byframe))
    df_data_byframe.insert(0, column='Frame', value=framelist)

    return df_data_raw, df_data, df_data_byframe, df_img

# Choose ROI from last image
# Read image
imgarray = cv2.imread(filelist[-1])
fromCenter = False  # Set up to choose as a drag-able rectangle rather than a rectangle chosen from center
r = cv2.selectROI("Image", imgarray, fromCenter)  # Choose ROI function from cv2 - opens a window to choose
x = int(r[0])  # Take result of selectROI and put into a variable
y = int(r[1])  # " "
w = int(r[2])  # " "
h = int(r[3])  # " "
cv2.destroyAllWindows()  # Destroy window when ready - requires any keyboard input


# Read last image
img = cv2.imread(filelist[-1])

crop = img[y:(y + h), x:(x + w), :]  # Create cropped image

# Generate map
# Create and display map
# Convert images to binary with thresholds, automatically uses all three colors
ret, img_th_red = cv2.threshold(crop[:, :, 2], rthresh, 255, cv2.THRESH_BINARY)  # Red
ret, img_th_green = cv2.threshold(crop[:, :, 1], gthresh, 255, cv2.THRESH_BINARY)  # Green
ret, img_th_blue = cv2.threshold(crop[:, :, 0], bthresh, 255, cv2.THRESH_BINARY)  # Blue

# Set up holder
# mapbin = np.zeros((self.h.get(), self.w.get()))
mapbin = np.zeros((h, w))
mapbin_ext = mapbin.copy()
# Threshold map
layered_arr = np.array([img_th_red, img_th_green, img_th_blue]).sum(axis=0)
mapbin[layered_arr >= 1] = 1
# Compress into one line
mapbin_1d = np.sum(mapbin, axis=1)
# Extend to width of channel
mapbin_ext[np.hstack(mapbin_1d * w) > 1] = 255

graphs = plt.figure()

# For each present color, run analysis_math for by-color spatial dataframe and mean, max dataframes
# Add to graph
if rchannel is True:
    df_data_raw, df_data, df_data_byframe, df_img = analysis_math(df_img, mapbin_ext, filelist, umpix, 2, rthresh, x, y,
                                                                  w, h)

    # Add column with color, append to overall dataframe
    df_data_raw.insert(0, 'Color', 'red')
    df_data.insert(0, 'Color', 'red')
    df_data_byframe.insert(0, 'Color', 'red')

    df = df.append(df_data_raw, ignore_index=True)
    df_summary_all = df_summary_all.append(df_data, ignore_index=True)
    df_summary_frame = df_summary_frame.append(df_data_byframe, ignore_index=True)
    df_colors = df_colors.append(df_img, ignore_index=True)

    # Add to graph
    # Plot each channel as light color
    for k in range(df_data['Channel'].max()):
        df_graph = df_data[df_data['Channel'] == k]
        # Occlusion
        plt.subplot(1, 2, 1)
        plt.plot(df_graph['Frame'], df_graph['Mean occlusion (%)'], color='salmon')
        # Accumulation
        plt.subplot(1, 2, 2)
        plt.plot(df_graph['Frame'], df_graph[u'Accumulation (\u03bcm\u00b2)'], color='salmon')
    # Plot mean as darker color
    # Occlusion
    plt.subplot(1, 2, 1)
    plt.plot(df_data_byframe['Frame'], df_data_byframe['Mean occlusion (%)'], color='red', linewidth='3',
             marker='3')
    # Accumulation
    plt.subplot(1, 2, 2)
    plt.plot(df_data_byframe['Frame'], df_data_byframe[u'Accumulation (\u03bcm\u00b2)'], color='red',
             linewidth='3', marker='3')
    # Titles, xlabels, ylabels
    # Occlusion
    plt.subplot(1, 2, 1)
    plt.title('Mean occlusion')
    plt.xlabel('Time point (n)')
    plt.ylabel('Percent channel occluded')
    # Accumulation
    plt.subplot(1, 2, 2)
    plt.title('Accumulation')
    plt.xlabel('Time point (n)')
    plt.ylabel(u'Accumulation (\u03bcm\u00b2)')

if gchannel is True:
    df_data_raw, df_data, df_data_byframe, df_img = analysis_math(df_img, mapbin_ext, filelist, umpix, 1, gthresh, x, y,
                                                                  w, h)

    # Add column with color, append to overall dataframe
    df_data_raw.insert(0, 'Color', 'green')
    df_data.insert(0, 'Color', 'green')
    df_data_byframe.insert(0, 'Color', 'green')
    df = df.append(df_data_raw, ignore_index=True)
    df_summary_all = df_summary_all.append(df_data, ignore_index=True)
    df_summary_frame = df_summary_frame.append(df_data_byframe, ignore_index=True)
    df_colors = df_colors.append(df_img, ignore_index=True)

    # Add to graph
    # Plot each channel as light color
    for k in range(df_data['Channel'].max()):
        df_graph = df_data[df_data['Channel'] == k]
        # Occlusion
        plt.subplot(1, 2, 1)
        plt.plot(df_graph['Frame'], df_graph['Mean occlusion (%)'], color='palegreen')
        # Accumulation
        plt.subplot(1, 2, 2)
        plt.plot(df_graph['Frame'], df_graph[u'Accumulation (\u03bcm\u00b2)'], color='palegreen')
    # Plot mean as darker color
    # Occlusion
    plt.subplot(1, 2, 1)
    plt.plot(df_data_byframe['Frame'], df_data_byframe['Mean occlusion (%)'], color='green', linewidth='3')
    # Accumulation
    plt.subplot(1, 2, 2)
    plt.plot(df_data_byframe['Frame'], df_data_byframe[u'Accumulation (\u03bcm\u00b2)'], color='green', linewidth='3')
    # Titles, xlabels, ylabels
    # Occlusion
    plt.subplot(1, 2, 1)
    plt.title('Mean occlusion')
    plt.xlabel('Time point (n)')
    plt.ylabel('Percent channel occluded')
    # Accumulation
    plt.subplot(1, 2, 2)
    plt.title('Accumulation')
    plt.xlabel('Time point (n)')
    plt.ylabel(u'Accumulation (\u03bcm\u00b2)')
    plt.tight_layout()

if bchannel is True:
    df_data_raw, df_data, df_data_byframe, df_img = analysis_math(df_img, mapbin_ext, filelist, umpix,
                                                                  0, bthresh, x, y, w, h)

    # Add column with color, append to overall dataframe
    df_data_raw.insert(0, 'Color', 'blue')
    df_data.insert(0, 'Color', 'blue')
    df_data_byframe.insert(0, 'Color', 'blue')
    df = df.append(df_data_raw, ignore_index=True)
    df_summary_all = df_summary_all.append(df_data, ignore_index=True)
    df_summary_frame = df_summary_frame.append(df_data_byframe, ignore_index=True)
    df_colors = df_colors.append(df_img, ignore_index=True)

    # Add to graph
    # Plot each channel as light color
    for k in range(df_data['Channel'].max()):
        df_graph = df_data[df_data['Channel'] == k]
        # Occlusion
        plt.subplot(1, 2, 1)
        plt.plot(df_graph['Frame'], df_graph['Mean occlusion (%)'], color='lightskyblue')
        # Accumulation
        plt.subplot(1, 2, 2)
        plt.plot(df_graph['Frame'], df_graph[u'Accumulation (\u03bcm\u00b2)'], color='lightskyblue')
    # Plot mean as darker color
    # Occlusion
    plt.subplot(1, 2, 1)
    plt.plot(df_data_byframe['Frame'], df_data_byframe['Mean occlusion (%)'], color='blue', linewidth='3')
    # Accumulation
    plt.subplot(1, 2, 2)
    plt.plot(df_data_byframe['Frame'], df_data_byframe[u'Accumulation (\u03bcm\u00b2)'], color='blue', linewidth='3')
    # Titles, xlabels, ylabels
    # Occlusion
    plt.subplot(1, 2, 1)
    plt.title('Mean occlusion')
    plt.xlabel('Time point (n)')
    plt.ylabel('Percent channel occluded')
    # Accumulation
    plt.subplot(1, 2, 2)
    plt.title('Accumulation')
    plt.xlabel('Time point (n)')
    plt.ylabel(u'Accumulation (\u03bcm\u00b2)')

# Create naming convention for excel sheet
if len(filelist) == 1:  # Single file
    nameconvention = os.path.basename(filelist[0]).split(".")[0]
elif len(filelist) > 1:  # Directory of files
    nameconvention = os.path.dirname(filelist[0]).split("/")[-1]

# Create writer to save results to
writer = pd.ExcelWriter(nameconvention[0:14] + '_analysis.xlsx', engine='openpyxl')  # Crop to avoid excel error

# Calculate summary metrics, write individual sheets to excel file writer
# Individual image names
unique_colors = []
if rchannel is True:
    unique_colors.append('red')
if gchannel is True:
    unique_colors.append('green')
if bchannel is True:
    unique_colors.append('blue')
# unique_colors = df.color.unique()

for uc in unique_colors:
    # Find all rows corresponding to unique name, three dataframes
    df_a = df[df['Color'] == uc]
    df_a.to_excel(writer, sheet_name=uc + ' raw data', index=False)  # Crop name to prevent errors
    df_summary_all_a = df_summary_all[df_summary_all['Color'] == uc]
    df_summary_all_a.to_excel(writer, sheet_name=uc + ' channel data', index=False)
    df_summary_frame_a = df_summary_frame[df_summary_frame['Color'] == uc]
    df_summary_frame_a.to_excel(writer, sheet_name=uc + ' frame data', index=False)

now = datetime.datetime.now()
# Print parameters to a sheet
param_df = pd.DataFrame({'Ratio, \u03bcm-to-pixels': umpix,
              'Red channel': rchannel,
              'Threshold, red channel': rthresh,
              'Green channel': gchannel,
              'Threshold, green channel': gthresh,
              'Blue channel': bchannel,
              'Threshold, blue channel': bthresh,
              'X coordinate (top)': x,
              'Y coordinate (top)': y,
              'ROI width': w,
              'ROI height': h,
              'Analysis date': now.strftime("%D"),
              'Analysis time': now.strftime("%H:%M:%S")}, index=[1])
param_df.to_excel(writer, sheet_name='Parameters used', index=False)

writer.save()
writer.close()

current_dir = os.getcwd()  # Select filepath

if current_dir.split('/')[-1] == 'Results, labeled image data':
    os.chdir(os.path.dirname(current_dir))

# Convert image to BGR
plt.savefig('Analysis_graph.png', dpi=300)

current_dir = os.getcwd()  # Select filepath

# if 'Results' in current_dir.split('/')[-1]:
#     current_dir = os.path.dirname(current_dir)
img_folder = os.path.join(current_dir, 'Results, labeled image data')

if os.path.exists(img_folder):
    shutil.rmtree(img_folder)

os.mkdir(img_folder)
os.chdir(img_folder)

for j in range(len(df_colors)):
    array = (df_colors['img'].iloc[j]).astype('uint8')
    cv2.imwrite(df_colors['name'].iloc[j] + '_' + df_colors['color'].iloc[j] + '_image.png', array)

# Finish code
writer.save()
writer.close()