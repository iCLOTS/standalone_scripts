# iCLOTS_analysis_methods
Script files containing all iCLOTS analysis methods.
These scripts have been designed specifically for use with iCLOTS software, a Lam lab project available at iCLOTS.org, but you might find it useful to adapt them to suit the needs of your own projects.

Methods rely on: Pandas, Numpy, OpenCV, skimage, Trackpy, seaborn, and more.
Each script has information about inputs, parameters, and outputs and a "tips" section within the module docstrings. Briefly, users select a folder of .png, .jpg., .tif, and/or .avi files for analysis. Users edit indicated parameters, some image processing step is applied, metric calculations are made and all output files are returned in a new directory within the original.

## Scripts included in repository
- accumulation_map.py: calculate occlusion and accumulation metrics for multiple blood components adhered to some device or surface over time, code adapts to an any-dimension microfluidic device
- accumulation_microchannel.py: calculate occlusion and accumulation metrics for multiple blood components adhered to some device or surface over time, code adapts to a series of horizontally-oriented channels, spatial metrics are returned
- accumulation_roi.py: calculate occlusion and accumulation metrics for multiple blood components adhered to some device or surface over time
- adhesion_brightfield.py: calculate morphological metrics describing individual cells within brightfield microscopy images
- adhesion_fluor.py: calculate morphological and functional metrics describing individual cells within fluoresence microscopy images
- adhesion_filopodia.py: a specialized version of the fluoresence microscopy adhesion application designed to count filopodia of cells
- adhesion_video.py: calculate morphological metrics and transit time values describing individual cells as they transiently adhere to some surface while under flow
- sct_brightfield.py: single cell tracking application to quantify morphological characteristics and velocity of individual cell(s) traveling in any direction(s) within brightfield videomicroscopy
- sct_fluor.py: single cell tracking application to quantify morphological characteristics, summed fluoresence intensity, and velocity of individual cell(s) traveling in any direction(s) within fluoresence videomicroscopy
- deformability_brightfield.py: a specialized version of the single cell tracking application to be used to calculate measures of mechanical properties of individual cells within brightfield microscopy data taken from experiments performed with the "biophysical flow cytometer" device, may also be used for more generalized channel flow x-direction tracking
- deformability_fluor.py: a specialized version of the single cell tracking application to be used to calculate measures of mechanical properties of individual cells within fluoresence microscopy data taken from experiments performed with the "biophysical flow cytometer" device, may also be used for more generalized channel flow x-direction tracking
- velocity.py: application designed to track velocity of a cell suspension as a function of time, creates a velocity flow profile of a video clip from all calculated displacements

## Inputs, outputs, methods
Users are guided to choose a directory of .png, .jpg, .tif, and/or .avi files using a file dialog window. File type depends on application. A separate suite of video processing tools in the video_processing repository can help format your data into the correct file type.

Users should edit input parameters based on their own individual needs. All parameter values requiring user editing are directly under import statements. Sample (from adhesion_brightfield.py):

```
# IMPORTANT: PARAMETERS TO EDIT
umpix = 1  # 1 = no conversion
max_diameter = 21  # Err on high side, MUST be an odd integer
min_mass = 1000  # Err on low side
invert=True  # True is dark on light, False is light on dark

```

Outputs typically include:
- Original image or video data labeled with all detected events (e.g. cell). Each label is an index that corresponds to numerical data within an excel sheet.
- Excel file with numerical data. This includes: one sheet per file of all events within that file, one sheet with all events from all files, a sheet with descriptive statistics (e.g. min., mean, max., st.dev. for all calculated metrics), and a sheet containing the parameters used and time/date of analysis for reference.
- Pairplot graphical data for each file, all files combined, and all files combined with each file a different color.
- All outputs are stored in an "Analysis" folder created within the original directory.

## Machine learning Jupyter notebook included in repository
- iCLOTS_ML_clustering.ipynb: machine learning interpretation for iCLOTS-generated data. Users select a folder containing one or more Excel files, select desired features, select desired number of clusters, and receive labeled outputs from a K-means clustering algorithm.

## Help and contributing
Contributions are always welcome! Submit a pull request or contact me directly at meredith.fay@iclots.org. Sample data and/or all data from the iCLOTS manuscript available upon request. Microfluidic mask designs available upon request. You can also find me at www.iCLOTS.org.
