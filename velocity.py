"""iCLOTS is a free software created for the analysis of common hematology workflow image data

Author: Meredith Fay, Lam Lab, Georgia Institute of Technology and Emory University
Last updated: 2022-07-27
This script corresponds to tools available in version 1.0b1, more recent implementations of tools
may be available within the iCLOTS software and in source code at github.com/iCLOTS

Script that analyzes videomicroscopy of cell suspensions under flow
A separate application (single cell tracking) analyzes individual cells under flow

Script relies heavily on OpenCV optical flow methods, documentation available at:
https://docs.opencv.org/3.4/d4/dee/tutorial_optical_flow.html

Input files:
--Script is designed to work a directory of videomicroscopy files (.avi)
----If your data is saved as a series of frames, please see the suite of video editing tools to convert to .avi

Input parameters, general:
--umpix: The ratio of microns (1e-6 m) to pixels for the video
----Use umpix = 1 for no conversion
--fps: Frames per second, the rate of imaging
----Note that FPS values pulled directly from videos can be inaccurate, especially if the video
-----has been resized or edited in any way
--n_bins: Number of bins to divide channel into for profile
----Typically aiming for roughly the size of a cell is best - e.g. enough bins so that each is 5 um
--max_corners: The maximum number of corners (features, typically a pattern of cells) to track
--labelimg_first: Boolean variable (True or False) indicating if you'd like labeled image data: first 100 frames
--labelimg_linspace: Boolean variable (True or False) indicating if you'd like labeled image data: every 100th frame

Input parameters, Shi-Tomasi corner finding:
----Choose highest n possible without adding unnecessary computational expense - I like ~500
--qual_level: The quality level - how "good" are the corners being tracked
----Choose a high "top percentage", like 0.01 (best 1%)
--min_dist_st: Minimum distance between Shi-Tomasi corners found
--block_size: Block size to compute eigenvalues over, I find roughly the dimensions of a cell is useful

Input parameters, Kanade-Lucas-Tomasi feature tracking:
--win_size: The most important parameter! Area (in pixels) to search for a feature in the subsequent frame (displacement)
----Err on the high side, as values to small can miss fastest moving features
----(set as y, x) Can reduce y to reduce computational expense if flow is uni-directional in x direction
----Keep in mind, even if flow is uni-directional, the algorithm will search in all directions
----This ties closely with the maximum pyramid level, below - e.g. for a pyramid level 2, divide window in mind by 2
--max_level: Typically set as 2. Reduces dimensions of image to reduce computational expense. See OpenCV docs.
--iter: Iterations to search for a feature. Typically 1 - want to see in each sequential frame.
--min_dist_klt: Distance a feature must travel in KLT algorithm. Typically 1 pixel is sufficient.

Output files
----If either/both labelimg parameter are true
------Specified original frames of the provided video with each detected event (cell cluster pattern)
-------labeled with an arrow indicating displacement

--A corresponding .csv sheet for each video containing location and velocity of all features tracked
----Typically a very large file

--A corresponding .xlsx sheet containing:
----Minimum, mean, and maximum velocity per frame - one sheet/video
----Descriptive statistics (mean, standard deviation for velocity) per video
------Minimum and maximum often aren't useful as some erroneously low and high values are always present
----Profile information: mean velocity and standard deviation of velocity (separate sheets)
------Calculated based on n_bins
------If profile seems erroneously blunted, increase window size - fastest values are being missed
----Parameters used and time/date analysis was performed, for reference

--Graphical data for each video
----Timecourse data, including minimum, maximum, and mean values per frame
----Profile data overlaid on each displacement data point

Some tips from the iCLOTS team:

Computational and experimental methods:
--All test data was taken at a frame rate of ~160 FPS, may need a fast camera to perform experiments
--Channel size was set artificially shallow (~15 um) - deeper channels are hard to get distinct features from

--Features are typically patterns of cells rather than a single individual cell
----As such, no label index is provided
--To reduce computational expense:
----Use relatively short video clips - 10 seconds is oftentimes sufficient to establish a profile
----Carefully choose window size - too high increases expense (but, too low misses critical displacement measurements)

Output files:
--Analysis files are named after the folder containing all videos (.xlsx) or video names (.png)
----Avoid spaces, punctuation, etc. within file names

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
# Labeling and plotting results
import matplotlib.pyplot as plt

# IMPORTANT: PARAMETERS TO EDIT
umpix = 70/80  # 1 = no conversion
fps = 160  # Frames per second rate of imaging
n_bins = 10  # Number of bins for profile
labelimg_first = False  # Save the first 100 frames with displacement labeled
labelimg_linspace = True  # Save every 100th frame with displacement labeled

# params for ShiTomasi corner detection
# to be used in line of code: cv2.goodFeaturesToTrack
# helpful: https://docs.opencv.org/master/d4/d8c/tutorial_py_shi_tomasi.html
max_corners = 500
qual_level = 0.01
min_dist_st = 5
block_size = 5
feature_params = dict( maxCorners = max_corners,  # n best corners to track, more = more computationally expensive
                       qualityLevel = qual_level,  # parameter characterizing the minimal accepted quality of image corners
                       minDistance = min_dist_st,  # minimum
                       # possible Euclidean distance between the returned corners (pix)
                       blockSize = block_size)  # size of an average block for computing a derivative covariation matrix over each pixel neighborhood

# Parameters for lucas kanade optical flow
# to be used in the line of code: cv2.calcOpticalFlowPyrLK
# helpful: https://docs.opencv.org/3.4/d4/dee/tutorial_optical_flow.html
win_size = (3, 10)
max_level = 2
iter = 1  # Numer of iterations (for search criteria)
min_dist_klt = 1  # Minimum distance that a corner must move (for search criteria)
lk_params = dict( winSize  = win_size,  # size of the search window at each pyramid level (w, h)
                  maxLevel = max_level,  # 0-based maximal pyramid level number, 0 = original, 1 = 2 levels used, 2 = 3 levels used
                  criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, iter, min_dist_klt))  # termination criteria of the iterative search algorithm

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

def descriptive_statistics(df_input):
    """Function to calculate descriptive statistics for each population, represented as a dataframe"""

    dict = {'n events tracked': len(df_input),
                  u'Mean velocity (\u03bcm/s)': df_input['Velocity (\u03bcm/s)'].mean(),
                  u'Stdev, velocity (\u03bcm/s)': df_input['Velocity (\u03bcm/s)'].std()
                  }
    dict_df = pd.DataFrame(dict, index=[0])

    return dict_df

def calculate_profile(df_input, h, n_bins):
    """"Function to calculate mean and stdev velocity values for a specified number of evenly distributed bins"""

    # Create a profile
    bins = np.linspace(0, h, n_bins)
    bins_um = bins * umpix  # For graphing
    digitized = np.digitize(df_input['Channel pos. (pix)'], bins)  # Bins
    profile = [df_input['Velocity (\u03bcm/s)'][digitized == i].mean() for i in range(1, len(bins))]

    # Create a profile (standard deviation)
    bins = np.linspace(0, h, n_bins)
    digitized = np.digitize(data_all['Channel pos. (pix)'], bins)  # Bins
    profile_stdev = [data_all['Velocity (\u03bcm/s)'][digitized == i].std() for i in range(1, len(bins))]

    return bins_um, profile, profile_stdev


# Set up combined, summary, and profile dataframes
df_summary = pd.DataFrame()
df_profiles = pd.DataFrame()
df_profiles_stdev = pd.DataFrame()

# For each video
for video in video_list:

    filename = os.path.basename(video).split(".")[0]  # Name of individual file
    os.chdir(output_folder)  # Return to original analysis folder
    if labelimg_first or labelimg_linspace is True:  # Make new output folder if images are being saved
        # Create a directory for that image
        dir_folder = output_folder + '/' + filename
        if os.path.exists(dir_folder):
            shutil.rmtree(dir_folder)  # Delete if already exists, prevents error
        os.makedirs(dir_folder)
        os.chdir(dir_folder)

    # Read video
    cap = cv2.VideoCapture(video)
    n_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    w = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))  # Record width of video frame (pixels)
    h = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))  # Record height of video frame (pixels)

    ret, first_frame = cap.read()

    # # Choose ROI from last frame (often initial frames have changes in illumination)
    # fromCenter = False  # Set up to choose as a drag-able rectangle
    # # r = cv2.selectROI("Image", first_frame, fromCenter)  # Choose ROI
    # ROI_x = int(r[0])  # Take result of selectROI and place into a variable
    # ROI_y = int(r[1])  # " "
    # ROI_w = int(r[2])  # " "
    # ROI_h = int(r[3])  # " "

    ROI_x = int(0)  # Take result of selectROI and place into a variable
    ROI_y = int(0)  # " "
    ROI_w = int(w)  # " "
    ROI_h = int(h)  # " "

    init_frame = first_frame[ROI_y:(ROI_y + ROI_h), ROI_x:(ROI_x + ROI_w), :]  # Create cropped image
    init_gray = cv2.cvtColor(init_frame, cv2.COLOR_BGR2GRAY)

    p0 = cv2.goodFeaturesToTrack(init_gray, mask=None, **feature_params)  # Initial points
    count = 0  # Count initial value

    # Initial lists to save frame, position, velocity values (dataframe too computationally expensive)
    frames = []
    positions = []
    velocities = []

    cap_ret = True
    while count < n_frames - 1 and cap_ret:
        cap_ret, frame = cap.read()  # Read
        frame_crop = frame[ROI_y:(ROI_y + ROI_h), ROI_x:(ROI_x + ROI_w), :]  # Create cropped image
        frame_gray = cv2.cvtColor(frame_crop, cv2.COLOR_BGR2GRAY)  # One layer

        # Select initial points
        p0 = cv2.goodFeaturesToTrack(frame_gray, mask=None, **feature_params)
        if p0 is not None:
            # Track points
            p1, st, err = cv2.calcOpticalFlowPyrLK(init_gray, frame_gray, p0, None, **lk_params)
            if p1.any():  # If points found
                good_new = p1[st==1]
                good_old = p0[st==1]

                for i,(new,old) in enumerate(zip(good_new, good_old)):
                    a,b = new.ravel()
                    c,d = old.ravel()

                    # Calculate displacement/velocity
                    vel = np.sqrt((new[0] - old[0])**2 + (new[1] - old[1])**2) * fps * umpix  # With y displacement
                    # vel = (new[0] - old[0]) * fps * umpix  # Ignore y displacement

                    # Save frame, position, displacement
                    frames.append(count)
                    positions.append(new[1])
                    velocities.append(vel)

                    # If user has indicated they'd like to save the first 100 frames
                    if labelimg_first is True:
                        if count < 100:
                            image = cv2.arrowedLine(frame_crop, (int(c),int(b)),(int(a),int(d)), (255, 255, 0), 1)  # Cyan arrow
                            cv2.imwrite(filename + '_frame_' + str(count).zfill(5) + '.png', frame_crop)

                    # If user has indicated they'd like to save every 100th frame
                    if labelimg_linspace is True:
                        if count % 100 == 0:
                            image = cv2.arrowedLine(frame_crop, (int(c),int(b)),(int(a),int(d)), (255, 255, 0), 1)  # Cyan arrow
                            cv2.imwrite(filename + '_frame_' + str(count).zfill(5) + '.png', frame_crop)

        count +=1

        # Now update the previous frame
        init_gray = frame_gray.copy()
        p0 = good_new.reshape(-1, 1, 2)

    # Create a dataframe with frame, position, and displacement data
    dict = {'Frame': frames, 'Channel pos. (pix)': positions, 'Velocity (\u03bcm/s)': velocities}
    data_all = pd.DataFrame(dict)  # Convert to dictionary

    # Save all data to a .csv file (all values typically too large for excel)
    data_all.to_csv(filename + '_all_data.csv')

    # Save minimum, mean, and maximum values per frame to an excel sheet
    # Group all tracked points by frame, take mean
    by_frame_min = data_all.groupby('Frame').min()
    by_frame_mean = data_all.groupby('Frame').mean()
    by_frame_max = data_all.groupby('Frame').max()

    by_frame_timepoint = np.linspace(0, len(by_frame_min), len(by_frame_min)) / fps

    by_frame_dataframe = pd.DataFrame({'Time (s)': by_frame_timepoint,
                                       'Min. velocity (\u03bcm/s)': by_frame_min['Velocity (\u03bcm/s)'],
                                       'Mean velocity (\u03bcm/s)': by_frame_mean['Velocity (\u03bcm/s)'],
                                       'Max. velocity (\u03bcm/s)': by_frame_max['Velocity (\u03bcm/s)']})

    # Write to excel
    by_frame_dataframe.to_excel(writer, sheet_name=filename[:30], index=False)  # Crop videoname to prevent excel error

    # Calculate profiles based on all events, height of channel, n bins
    bins_um, profile, profile_stdev = calculate_profile(data_all, ROI_h, n_bins)
    # Save profile to dataframe
    df_profiles[filename] = profile
    # Save standard deviation profile to dataframe
    df_profiles_stdev[filename] = profile_stdev

    # Create graphs
    # Timecourse graph
    fig = plt.figure()
    plt.scatter(by_frame_timepoint, by_frame_min['Velocity (\u03bcm/s)'], color='springgreen', label='Min.')
    plt.scatter(by_frame_timepoint, by_frame_mean['Velocity (\u03bcm/s)'], color='dodgerblue', label='Mean')
    plt.scatter(by_frame_timepoint, by_frame_max['Velocity (\u03bcm/s)'], color='tomato', label='Max.')
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (\u03bcm/s)')
    plt.savefig(filename + '_timecourse.png', dpi=300)
    plt.close()

    # Profile graph
    fig = plt.figure()
    plt.scatter(data_all['Channel pos. (pix)'] * umpix, data_all['Velocity (\u03bcm/s)'], color='lightskyblue')
    plt.plot(bins_um[1:], profile, color='dodgerblue')
    plt.xlabel('Channel position (\u03bcm)')
    plt.ylabel('Velocity (\u03bcm/s)')
    plt.savefig(filename + '_profile.png', dpi=300)
    plt.close()

    # Add descriptive statistics for video to summary sheet
    df_video = descriptive_statistics(data_all)
    df_video.insert(0, 'Video', filename)
    df_summary = df_summary.append(df_video, ignore_index=True)

os.chdir(output_folder)

# Save profiles to individual sheets
df_profiles.to_excel(writer, sheet_name='Velocity profiles', index=False)
df_profiles_stdev.to_excel(writer, sheet_name='Velocity profiles, stdev.', index=False)

# Save parameters
param_df = pd.DataFrame({u'Ratio, \u03bcm-to-pixels': umpix,
                         'FPS': fps,
                         'Histogram bins (n)': n_bins,
                         'Max. corners': max_corners,
                         'Quality level': qual_level,
                         'Min. dist., Shi-Tomasi': min_dist_st,
                         'Block size': block_size,
                         'Window size': str(win_size[0]) + ", " + str(win_size[1]),
                         'KLT iterations': iter,
                         'Min. dist. KLT': min_dist_klt,
                         'Pyramid level': max_level,
                         'Label, first 100': labelimg_first,
                         'Label, each 100th': labelimg_linspace,
                         'Analysis date': now.strftime("%D"),
                         'Analysis time': now.strftime("%H:%M:%S")}, index=[1])
param_df.to_excel(writer, sheet_name='Parameters used', index=False)

# Save and close excel file writer
writer.save()
writer.close()