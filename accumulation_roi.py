"""iCLOTS is a free software created for the analysis of common microfluidic and hematology workflow image data

Author: Meredith Fay, Lam Lab, Georgia Institute of Technology and Emory University
Last updated: 2022-10-29
This script corresponds to tools available in version 1.0b1, more recent implementations of tools
may be available within the iCLOTS software and in source code at github.com/iCLOTS

Region of interest (ROI) accumulation/occlusion analysis GUI analyzes single or timeseries images
of cells adhered to a rectangular surface of portion of a microfluidic device
Accomodates red, blue, and/or green channels

Application relies heavily on skimage library region property analysis, tutorial available at:
https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_regionprops.html

Input files:
--Application is designed to work with a folder of images (.jpg, .png, and/or .tif)
----The same input parameters are applied to each image
----This code will work on a single image
----It's important that frames are labeled sequentially in timeseries order (i.e. 01, 02.. 10 vs. 1, 2... 10)
----After uploading one or several images, the user is prompted to choose an ROI from the first image
------The same ROI is applied to all images, take care that all images represent the same field of view

Input parameters:
--umpix: micron to pixel ratio: The ratio of microns (1e-6 m) to pixels for the image
----Use value = 1 for no conversion
--Presence of fluorescent signal (rchannel, gchannel, bchannel) - boolean stating whether or not to analyze a channel
--Red, green, and/or blue threshold(s) (rthresh, gthresh, bthresh) - [0, 255] threshold value for each channel

Output files
--ROI from last image
--Images with image processing (threshold) steps applied

--A corresponding .xlsx sheet containing:
----For each selected channel:
------Frame data: mean signal area, signal accumulation (all channels)
------Conversion notes:
--------To convert accumulation per frame into per timepoint, divide frame # by FPS imaging rate

--Future versions of iCLOTS will include obstruction, the percent of the y-direction for each channel
---containing some signal

--Occlusion/accumulation graph:
----For the time series, a line graph showing:
------Occlusion (titled, left) for each channel (light lines) and mean (dark lines) for each color
------Accumulation (titled, right) - " "

Some tips from the iCLOTS team:

Computational and experimental methods:
--See input requirements: a time series, in the same field of view
----We are planning a brightfield/fluoresence microscopy combo. app for future iCLOTS releases
--Time series images must be in the proper alphabetical/numerical order
----If image names contain numbers, use preceding zeros to order properly
------i.e. 01, 02... 10 instead of 1, 2... 10
--The Lam lab has developed these methods using an ibidi slide
----We are happy to share ordering information upon request

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
import matplotlib.pyplot as plt
import datetime

# IMPORTANT: PARAMETERS TO EDIT
umpix = 0.5  # micron-to-pixel ratio
rchannel = True  # Present channels
gchannel = True
bchannel = False
rthresh = 5  # Channel thresholds
gthresh = 5
bthresh = 5

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

# Choose ROI
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

crop = imgarray[y:(y + h), x:(x + w), :]  # Create cropped image

cv2.imwrite(foldername + '_ROI.png', crop)

def occ_acc(name, thresh, layer):
    occlusion = [0]  # init
    accumulation = []
    time = [0]

    for img in filelist:

        imgname = os.path.basename(img).split(".")[0]
        channelimg = cv2.imread(img)

        channelimg = channelimg[:, :, layer]  # Pull out correct layer

        crop = channelimg[y:(y + h), x:(x + w)]  # Create cropped image

        ret, array_bin = cv2.threshold(crop, thresh, 255, cv2.THRESH_BINARY)
        cv2.imwrite(imgname + '_' + name + '_threshold-applied.png', array_bin)
        occ = np.sum(array_bin/255)

        time.append(imgname)
        occlusion.append(occ)

    # Calculate accumulation as change in occlusion between frames
    for i in range(len(occlusion)-1):
        accumulation.append(occlusion[i+1] - occlusion[i])

    occdata = pd.DataFrame(
        {'Image': time,
         'Occlusion (pix)': occlusion,
         'Occlusion (\u03bcm\u00b2)': np.asarray(occlusion) * umpix * umpix
         })
    occdata.to_excel(writer, sheet_name=name + ' occlusion', index=False)

    accdata = pd.DataFrame(
        {'Image': time[1:],
         'Accumulation (pix/timepoint)': accumulation,
         'Accumulation (\u03bcm\u00b2/timepoint)': np.asarray(accumulation) * umpix * umpix
         })
    accdata.to_excel(writer, sheet_name=name + ' accumulation', index=False)

    # Make some graphs
    # Occlusion
    fig = plt.figure()
    timevec = range(len(time))
    plt.plot(timevec, np.asarray(occlusion) * umpix * umpix, color=name)
    plt.title(name + ' occlusion over time')
    plt.xlabel('Timepoint (n)')
    plt.xticks(rotation=90)
    plt.ylabel('Occlusion (\u03bcm\u00b2)')
    plt.tight_layout()

    graphname = foldername + '_' + name + '_occlusion_graph.png'
    plt.savefig(graphname, dpi=300)

    plt.close()

    fig = plt.figure()
    timevec = range(len(time[1:]))
    plt.plot(timevec, np.asarray(accumulation) * umpix * umpix, color=name)
    plt.title(name + ' accumulation over time')
    plt.xlabel('Timepoint (n)')
    plt.xticks(rotation=90)
    plt.ylabel('Accumulation (\u03bcm\u00b2)/timepoint')
    plt.tight_layout()

    graphname = foldername + '_' + name + '_accumulation_graph.png'
    plt.savefig(graphname, dpi=300)
    plt.close()


if rchannel is True:
    occ_acc('red', rthresh, 2)  # Name, threshold, layer
if gchannel is True:
    occ_acc('green', gthresh, 1)
if bchannel is True:
    occ_acc('blue', bthresh, 0)

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