##### INPUTS #####
detections_list ='all_detections_1330_1380.dat'  #example '/users/aycha/MIGHTEE/notebooks/1380_final_list.dat' 
                     #Please input the path to your detection file list in an ASCII table format (e.g xxx.dat) with columns RA DEC FREQ
                     #Where ra, dec = 02h17m40.0s, 04d33m40.0s or in decimal coordinates 33.44444, 04.55555
                     #frequency in GHz. 
path_to_cube='/idia/users/aycha/pbcorrected/XMMLSS12_1_2.1330.UVCON1.CUBE.R0p5.pbcorr.fits'      #e.g '/carta_share/groups/mightee/Sambatra/XMMLSS12_highres.1380.UVSUB.N5UVCONTSUB.NOGRID.ms.CUBE.PBCORR.FITS'      
                     #Please input the path to the FITS cube

#Extraction parameters

subcube_width=2 #in arcmin
subcube_width_in_text='2arcmin'
# However, we may encounter detections with bigger size, so you can select their index numbers:
index = []         # [1] for example in my detection, the detection with index 1 is bigger and should have a width of 4 arcminutes.
width_new=4
width_new_in_text='4arcmin'
velocity_convention='radio' #you can also change to 'optical'
half_freq_width=1e-3        #it is better to select a small frequency range for making moment maps, with this a 2MHz subcube will be extracted

#Beam convolution parameters

circular_beam_axis = 20 # 20'' x 20''

#Moment maps parameters
sigma=3                #Masking
mean_profile_width=300 #km/s, needed for moment 1 map
dirName = 'Moment_maps_1330_N1_R0p5_PBCORR_sorted' #Path of the Directory for moment maps to be created

print('Variable recorded successfully, ready to run')