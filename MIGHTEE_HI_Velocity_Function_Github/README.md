Notebooks that layout how to Calculate the velocity function of MIGHTEE-HI Data.

Step 1: Calculate channel maps, moment 0 and moment 1 maps from the detections list. Calculate inclination from the HI disks.

Step 2: Extract the HI profiles from the data cubes. Write the spectral values to text files.

Step 3: Import spectra from step 2 and fit Busy Functions using PyMultiNest. Need to extract W50.

Step 4: Calculate the velocity function (volume density as a function of W50).

Step 5: Fit the velocity functions using a modified Schechter function.