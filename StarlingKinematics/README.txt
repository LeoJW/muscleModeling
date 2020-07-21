####——README——####

This is all the junk you need to use DLTdv7 (and easyWand5) to digitize flights from my Montana starling wind tunnel flights. 


Each folder in this head directory corresponds to a given bird (I’m dumb so I catalogued them by the date of the flight recordings ¯\_(ツ)_/¯). 

Also in this head directory is a text file called “camProfiles_profile.txt”. This contains the information easyWand5 needs to know about all the cameras, including their focal lengths and image sensor size and such.


Inside each folder there will be a ton of different files. 

Some of these are just for calibrations: If they have “calib” in the name, these are the calibration files to use when making a new file. Refer to the "calib_XX_dltCoefs.csv" files when asked for the calibration in DLTdv7.

The rest are data from individual flights, referred to as "runs". The .mat file stores all the data and metadata DLTdv7 needs to do its thing, and the .csv files are the exported output data from DLTdv7. The "rXX_xyzpts.csv" file is the one with the 3d kinematics data in the final form.

