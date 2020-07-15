import spectral_cube
from spectral_cube import SpectralCube
import numpy as np
import os
import pylab as pl
from astropy import units as u
from astropy import wcs
import aplpy
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

#########################################################################################################################################################
####                                                             INPUTS                                                                              ####
#########################################################################################################################################################
#1. Path to the continuum cubes
#e.g. '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547413577/V0_T29R00C02_3-MFS-image.fits'
#Continuum fields
c1 = '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547499081/V0_T29R03C01_3-MFS-image.fits'
c2 = '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547499081/V0_T29R03C03_3-MFS-image.fits'
c3 = '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547499081/V0_T29R03C05_3-MFS-image.fits'
c4 = '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547499081/V0_T29R04C02_3-MFS-image.fits'
c5 = '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547499081/V0_T29R04C04_3-MFS-image.fits'
c6 = '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547499081/V0_T29R04C06_3-MFS-image.fits'
c7 = '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547499081/V0_T30R00C02_3-MFS-image.fits'
c8 = '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547499081/V0_T30R00C04_3-MFS-image.fits'
c9 = '/idia/projects/vela/V0_GPS_results/0V0_continuum/0V0_continuum_1547499081/V0_T30R00C06_3-MFS-image.fits'
cube_cont1 = SpectralCube.read(c1)
cube_cont2 = SpectralCube.read(c2)
cube_cont3 = SpectralCube.read(c3)
cube_cont4 = SpectralCube.read(c4)
cube_cont5 = SpectralCube.read(c5)
cube_cont6 = SpectralCube.read(c6)
cube_cont7 = SpectralCube.read(c7)
cube_cont8 = SpectralCube.read(c8)
cube_cont9 = SpectralCube.read(c9)

print('Continuum cubes loaded')

#2. Path to the HI cubes
# e.g. '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547413577/0123456789_V0_T29R00C02_HI.image.fits'
#Continuum fields
h1 = '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547499081/0123456789_V0_T29R03C01_HI.image.fits'
h2 = '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547499081/0123456789_V0_T29R03C03_HI.image.fits'
h3 = '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547499081/0123456789_V0_T29R03C05_HI.image.fits'
h4 = '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547499081/0123456789_V0_T29R04C02_HI.image.fits'
h5 = '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547499081/0123456789_V0_T29R04C04_HI.image.fits'
h6 = '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547499081/0123456789_V0_T29R04C06_HI.image.fits'
h7 = '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547499081/0123456789_V0_T30R00C02_HI.image.fits'
h8 = '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547499081/0123456789_V0_T30R00C04_HI.image.fits'
h9 = '/idia/projects/vela/V0_GPS_results/1V0_HI_cube/1V0_HI_cube_1547499081/0123456789_V0_T30R00C06_HI.image.fits'
cube_h1 = SpectralCube.read(h1)
cube_h2 = SpectralCube.read(h2)
cube_h3 = SpectralCube.read(h3)
cube_h4 = SpectralCube.read(h4)
cube_h5 = SpectralCube.read(h5)
cube_h6 = SpectralCube.read(h6)
cube_h7 = SpectralCube.read(h7)
cube_h8 = SpectralCube.read(h8)
cube_h9 = SpectralCube.read(h9)

print('HI cubes loaded')

#3. Field names:
obs_name ='1547499081'
field = ['T29R03C01','T29R03C03','T29R03C05','T29R04C02','T29R04C04','T29R04C06','T30R00C02','T30R00C04','T30R00C06']


#######################################################################################################################################################


# READY TO RUN
#Output = PV images, continuum map images in png

#The cubes are all the same
#500 pixels in the center
sub_h1 = cube_h1[:,1200-250:1200+250,1200-250:1200+250]
sub_h2 = cube_h2[:,1200-250:1200+250,1200-250:1200+250]
sub_h3 = cube_h3[:,1200-250:1200+250,1200-250:1200+250]
sub_h4 = cube_h4[:,1200-250:1200+250,1200-250:1200+250]
sub_h5 = cube_h5[:,1200-250:1200+250,1200-250:1200+250]
sub_h6 = cube_h6[:,1200-250:1200+250,1200-250:1200+250]
sub_h7 = cube_h7[:,1200-250:1200+250,1200-250:1200+250]
sub_h8 = cube_h8[:,1200-250:1200+250,1200-250:1200+250]
sub_h9 = cube_h9[:,1200-250:1200+250,1200-250:1200+250]

print('500 pixels in the center of HI cubes selected')

## PV maps HI cubes
cube = [sub_h1,sub_h2,sub_h3,sub_h4,sub_h5,sub_h6,sub_h7,sub_h8,sub_h9]

fig, axs = plt.subplots(3, 3, figsize=(20, 20), sharex=True, sharey=True,
                        gridspec_kw={'hspace': -0.0, 'wspace': -0.1})

for i in range(len(axs.flat)):
    x = np.fliplr(cube[i][:,:, 250].array) #Inverting the declination with fliplr
    im = axs.flat[i].imshow(x.T, vmin=5e-6, vmax=5e-4) #RA fixed at the cube center, transpose x and y axis
    #500 pixels in the center for declination
    #vmin and vmax range of colors
    if i == 0 or i == 3:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
    elif i == 7 or i ==8:
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    elif i ==6:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    axs.flat[i].minorticks_on()
    axs.flat[i].tick_params(which='major', length=10, width=2, direction='in',color="white")
    axs.flat[i].tick_params(which='minor', length=4, width=2, direction='in',color="white")
    axs.flat[i].yaxis.set_ticks_position('both')
    axs.flat[i].xaxis.set_ticks_position('both')
    plt.setp(axs.flat[i].get_xticklabels(), visible=False)
    plt.setp(axs.flat[i].get_yticklabels(), visible=False)
    anchored_text = AnchoredText(field[i],loc=2, borderpad=0.5,frameon=False,
                                prop={'family': 'serif', 'size': 16, 'fontweight': 'normal', 'color':'white'})
    axs.flat[i].add_artist(anchored_text)
    #axs.flat[i].text(450, 25, field[i], size=30, rotation=0., family='serif', color = 'white',
         #ha="right", va="top",
         #bbox={'boxstyle': 'square','facecolor': 'black', 'alpha': 0.5}
         #)
pl.subplots_adjust(wspace=-0.5, hspace=0.00)
fig.suptitle('500 pix, vmin=5e-6, vmax = 5e-4', fontsize=25)
#pl.tight_layout()

cbar = fig.colorbar(im, ax=axs.flat,pad=0.05,orientation='horizontal')
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_xlabel(r'Flux (Jy/beam)',fontsize=20)

#Saving:
if not os.path.exists(obs_name):
    os.makedirs(obs_name)
    print("Directory " , obs_name,  " Created ")
    pl.savefig(obs_name+'/PV_map_HI_cube_500_pix_'+obs_name+'.png',overwrite=True)
else:
    pl.savefig(obs_name+'/PV_map_HI_cube_500_pix_'+obs_name+'.png',overwrite=True)
    
## PV map different levels
fig, axs = plt.subplots(3, 3, figsize=(20, 20), sharex=True, sharey=True,
                        gridspec_kw={'hspace': -0.0, 'wspace': -0.1})

for i in range(len(axs.flat)):
    x = np.fliplr(cube[i][:,:, 250].array) #Inverting the declination with fliplr
    im = axs.flat[i].imshow(x.T, vmin=-1e-3, vmax=1e-3) #RA fixed at the cube center, transpose x and y axis
    #500 pixels in the center for declination
    #vmin and vmax range of colors
    if i == 0 or i == 3:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
    elif i == 7 or i ==8:
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    elif i ==6:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    axs.flat[i].minorticks_on()
    axs.flat[i].tick_params(which='major', length=10, width=2, direction='in',color="white")
    axs.flat[i].tick_params(which='minor', length=4, width=2, direction='in',color="white")
    axs.flat[i].yaxis.set_ticks_position('both')
    axs.flat[i].xaxis.set_ticks_position('both')
    plt.setp(axs.flat[i].get_xticklabels(), visible=False)
    plt.setp(axs.flat[i].get_yticklabels(), visible=False)
    anchored_text = AnchoredText(field[i],loc=2, borderpad=0.5,frameon=False,
                                prop={'family': 'serif', 'size': 16, 'fontweight': 'normal', 'color':'white'})
    axs.flat[i].add_artist(anchored_text)
    #axs.flat[i].text(450, 25, field[i], size=30, rotation=0., family='serif', color = 'white',
         #ha="right", va="top",
         #bbox={'boxstyle': 'square','facecolor': 'black', 'alpha': 0.5}
         #)
pl.subplots_adjust(wspace=-0.5, hspace=0.00)
fig.suptitle('500 pix, vmin=-1e-3, vmax = 1e-3', fontsize=25)
#pl.tight_layout()

cbar = fig.colorbar(im, ax=axs.flat,pad=0.05,orientation='horizontal')
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_xlabel(r'Flux (Jy/beam)',fontsize=20)

pl.savefig(obs_name+'/PV_map_HI_cube_500_pix_-1e3_1e3'+obs_name+'.png',overwrite=True)
    
## PV map different levels
fig, axs = plt.subplots(3, 3, figsize=(20, 20), sharex=True, sharey=True,
                        gridspec_kw={'hspace': -0.0, 'wspace': -0.1})

for i in range(len(axs.flat)):
    x = np.fliplr(cube[i][:,:, 250].array) #Inverting the declination with fliplr
    im = axs.flat[i].imshow(x.T, vmin=-6e-4, vmax=6e-4) #RA fixed at the cube center, transpose x and y axis
    #500 pixels in the center for declination
    #vmin and vmax range of colors
    if i == 0 or i == 3:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
    elif i == 7 or i ==8:
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    elif i ==6:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    axs.flat[i].minorticks_on()
    axs.flat[i].tick_params(which='major', length=10, width=2, direction='in',color="white")
    axs.flat[i].tick_params(which='minor', length=4, width=2, direction='in',color="white")
    axs.flat[i].yaxis.set_ticks_position('both')
    axs.flat[i].xaxis.set_ticks_position('both')
    plt.setp(axs.flat[i].get_xticklabels(), visible=False)
    plt.setp(axs.flat[i].get_yticklabels(), visible=False)
    anchored_text = AnchoredText(field[i],loc=2, borderpad=0.5,frameon=False,
                                prop={'family': 'serif', 'size': 16, 'fontweight': 'normal', 'color':'white'})
    axs.flat[i].add_artist(anchored_text)
    #axs.flat[i].text(450, 25, field[i], size=30, rotation=0., family='serif', color = 'white',
         #ha="right", va="top",
         #bbox={'boxstyle': 'square','facecolor': 'black', 'alpha': 0.5}
         #)
pl.subplots_adjust(wspace=-0.5, hspace=0.00)
fig.suptitle('500 pix, vmin=-6e-4, vmax = 6e-4', fontsize=25)
#pl.tight_layout()

cbar = fig.colorbar(im, ax=axs.flat,pad=0.05,orientation='horizontal')
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_xlabel(r'Flux (Jy/beam)',fontsize=20)

pl.savefig(obs_name+'/PV_map_HI_cube_500_pix_-6e4_6e4'+obs_name+'.png',overwrite=True)

print('-6e4 to 6e4 levels plotted!')

## PV map different levels
fig, axs = plt.subplots(3, 3, figsize=(20, 20), sharex=True, sharey=True,
                        gridspec_kw={'hspace': -0.0, 'wspace': -0.1})

for i in range(len(axs.flat)):
    x = np.fliplr(cube[i][:,:, 250].array) #Inverting the declination with fliplr
    im = axs.flat[i].imshow(x.T, vmin=-9e-4, vmax=9e-4) #RA fixed at the cube center, transpose x and y axis
    #500 pixels in the center for declination
    #vmin and vmax range of colors
    if i == 0 or i == 3:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
    elif i == 7 or i ==8:
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    elif i ==6:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    axs.flat[i].minorticks_on()
    axs.flat[i].tick_params(which='major', length=10, width=2, direction='in',color="white")
    axs.flat[i].tick_params(which='minor', length=4, width=2, direction='in',color="white")
    axs.flat[i].yaxis.set_ticks_position('both')
    axs.flat[i].xaxis.set_ticks_position('both')
    plt.setp(axs.flat[i].get_xticklabels(), visible=False)
    plt.setp(axs.flat[i].get_yticklabels(), visible=False)
    anchored_text = AnchoredText(field[i],loc=2, borderpad=0.5,frameon=False,
                                prop={'family': 'serif', 'size': 16, 'fontweight': 'normal', 'color':'white'})
    axs.flat[i].add_artist(anchored_text)
    #axs.flat[i].text(450, 25, field[i], size=30, rotation=0., family='serif', color = 'white',
         #ha="right", va="top",
         #bbox={'boxstyle': 'square','facecolor': 'black', 'alpha': 0.5}
         #)
pl.subplots_adjust(wspace=-0.5, hspace=0.00)
fig.suptitle('500 pix, vmin=-9e-4, vmax = 9e-4', fontsize=25)
#pl.tight_layout()

cbar = fig.colorbar(im, ax=axs.flat,pad=0.05,orientation='horizontal')
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_xlabel(r'Flux (Jy/beam)',fontsize=20)

pl.savefig(obs_name+'/PV_map_HI_cube_500_pix_9e4_9e4'+obs_name+'.png',overwrite=True)
    
## PV map for the whole cube:
cube = [cube_h1,cube_h2,cube_h3,cube_h4,cube_h5,cube_h6,cube_h7,cube_h8,cube_h9]

fig, axs = plt.subplots(3, 3, figsize=(15, 40), sharex=True, sharey=True,
                        gridspec_kw={'hspace': -0.0, 'wspace': -0.5})

for i in range(len(axs.flat)):
    x = np.fliplr(cube[i][:,:, 1200].array) #Inverting the declination with fliplr
    im = axs.flat[i].imshow(x.T, vmin=5e-6, vmax=5e-4) #RA fixed, transpose x and y axis
    #500 pixels in the center for declination
    #vmin and vmax range of colors
    if i == 0 or i == 3:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
    elif i == 7 or i ==8:
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    elif i ==6:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
        axs.flat[i].set_xlabel('Spectral axis',fontsize=20,family='serif')
    axs.flat[i].minorticks_on()
    axs.flat[i].tick_params(which='major', length=10, width=2, direction='in',color="white")
    axs.flat[i].tick_params(which='minor', length=4, width=2, direction='in',color="white")
    axs.flat[i].yaxis.set_ticks_position('both')
    axs.flat[i].xaxis.set_ticks_position('both')
    plt.setp(axs.flat[i].get_xticklabels(), visible=False)
    plt.setp(axs.flat[i].get_yticklabels(), visible=False)
    anchored_text = AnchoredText(field[i],loc=2, borderpad=0.5,frameon=False,
                                prop={'family': 'serif', 'size': 16, 'fontweight': 'normal', 'color':'white'})
    axs.flat[i].add_artist(anchored_text)
pl.subplots_adjust(wspace=-0.5, hspace=0.00)
fig.suptitle('all DEC, vmin=5e-6, vmax =5e-4', fontsize=25)
#pl.tight_layout()

cbar = fig.colorbar(im, ax=axs.flat,pad=0.05,orientation='horizontal')
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_xlabel(r'Flux (Jy/beam)',fontsize=20)

pl.savefig(obs_name+'/PV_map_HI_cube_all_pix_'+obs_name+'.png',overwrite=True)

print('The whole cube selected!')

## HI cube in DEC vs RA
fig, axs = plt.subplots(3, 3, figsize=(20, 20), sharex=True, sharey=True,
                        gridspec_kw={'hspace': -0.05, 'wspace': -0.05})

for i in range(len(axs.flat)):
    x = np.flip(cube[i][250,:, :].array,0) #Inverting the declination with fliplr, channel fixed at 250
    im = axs.flat[i].imshow(x, vmin=5e-6, vmax=5e-4) #no transposition
    #all the pixels in the cube
    #vmin and vmax range of colors
    if i == 0 or i == 3:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
    elif i == 7 or i ==8:
        axs.flat[i].set_xlabel('RA(J2000)',fontsize=20,family='serif')
    elif i ==6:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
        axs.flat[i].set_xlabel('RA(J2000)',fontsize=20,family='serif')
    axs.flat[i].minorticks_on()
    axs.flat[i].tick_params(which='major', length=10, width=2, direction='in',color="white")
    axs.flat[i].tick_params(which='minor', length=4, width=2, direction='in',color="white")
    axs.flat[i].yaxis.set_ticks_position('both')
    axs.flat[i].xaxis.set_ticks_position('both')
    plt.setp(axs.flat[i].get_xticklabels(), visible=False)
    plt.setp(axs.flat[i].get_yticklabels(), visible=False)
    anchored_text = AnchoredText(field[i],loc=2, borderpad=0.5,frameon=False,
                                prop={'family': 'serif', 'size': 16, 'fontweight': 'normal', 'color':'white'})
    axs.flat[i].add_artist(anchored_text)
    #axs.flat[i].text(2000, 100, field[i], size=20, rotation=0., family='serif', color = 'white',
         #ha="right", va="top",
         #bbox={'boxstyle': 'square','facecolor': 'black', 'alpha': 0.5}
         #)
pl.subplots_adjust(wspace=-0.0, hspace=0.00)
fig.suptitle('HI map, vmin=5e-6, vmax = 5e-4', fontsize=25)
#pl.tight_layout()

cbar = fig.colorbar(im, ax=axs.flat,pad=0.05,orientation='horizontal')
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_xlabel(r'Flux (Jy/beam)',fontsize=20)

pl.savefig(obs_name+'/HI_map_'+obs_name+'.png',overwrite=True)

print('HI map plotted!')
    
## Continuum maps
cube = [cube_cont1,cube_cont2,cube_cont3,cube_cont4,cube_cont5,cube_cont6,cube_cont7,cube_cont8,cube_cont9]

fig, axs = plt.subplots(3, 3, figsize=(20, 20), sharex=True, sharey=True,
                        gridspec_kw={'hspace': -0.05, 'wspace': -0.05})

for i in range(len(axs.flat)):
    x = np.flip(cube[i][0,:, :].array,0) #Inverting the declination with fliplr, channel fixed
    im = axs.flat[i].imshow(x, vmin=-0.00005, vmax=0.00005) #no transposition
    #all the pixels in the cube
    #vmin and vmax range of colors
    if i == 0 or i == 3:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
    elif i == 7 or i ==8:
        axs.flat[i].set_xlabel('RA(J2000)',fontsize=20,family='serif')
    elif i ==6:
        axs.flat[i].set_ylabel('DEC(J2000)',fontsize=20,family='serif')
        axs.flat[i].set_xlabel('RA(J2000)',fontsize=20,family='serif')
    axs.flat[i].minorticks_on()
    axs.flat[i].tick_params(which='major', length=10, width=2, direction='in',color="white")
    axs.flat[i].tick_params(which='minor', length=4, width=2, direction='in',color="white")
    axs.flat[i].yaxis.set_ticks_position('both')
    axs.flat[i].xaxis.set_ticks_position('both')
    plt.setp(axs.flat[i].get_xticklabels(), visible=False)
    plt.setp(axs.flat[i].get_yticklabels(), visible=False)
    anchored_text = AnchoredText(field[i],loc=2, borderpad=0.5,frameon=False,
                                prop={'family': 'serif', 'size': 16, 'fontweight': 'normal', 'color':'black'})
    axs.flat[i].add_artist(anchored_text)
    #axs.flat[i].text(2000, 100, field[i], size=20, rotation=0., family='serif', color = 'white',
         #ha="right", va="top",
         #bbox={'boxstyle': 'square','facecolor': 'black', 'alpha': 0.5}
         #)
pl.subplots_adjust(wspace=-0.0, hspace=0.00)
fig.suptitle('Cont map, vmin=-5e-5, vmax = 5e-5', fontsize=25)
#pl.tight_layout()

cbar = fig.colorbar(im, ax=axs.flat,pad=0.05,orientation='horizontal')
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_xlabel(r'Flux (Jy/beam)',fontsize=20)

pl.savefig(obs_name+'/Continuum_map_'+obs_name+'.png',overwrite=True)

print('Code run succesfully!')