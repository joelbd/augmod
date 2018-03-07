# Make sure needed packages are imported
import os
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd


# Temporarily add includes to $PATH
# os.system("export PATH=$PATH:/includes")

# For more information on Galfit paramters see:
# https://users.obs.carnegiescience.edu/peng/work/galfit/README.pdf

# Take user input parameters and write a galfit feedme file.
# Input paramters as defined in the example above
def galfit(
    file_input, # The path to the fits image to model
    file_output, # The desired path for output files
    image_region, # The size of the image file
    plate_scale, # Given in format: "<float>  <float>"
    profile_type, # See galfit documentation for options
    pos_xy, # XY location for center of object to be modeled
    integrated_mag, # Estimated integrated magnitude
    eff_radius, # Estimated effective radius
    sersic_index, # Sersic index if using sersic profile
    axis_ratio, # Axis ration of the galaxy to be modeled

    #Optional parameters
    sigma_image = "none", # Not required for succesful modeling
    position_angle = 0,
    display_type = "regular",
    psf_image = "none", # Not required for successful modeling
    psf_sample_factor = 1,
    pixel_mask_image = "none", # Given as a fits or ascii file of coords
    constraints_file = "none", # Not typically used, exercise caution
    convolution_box = "50  50", # Size of convolution box.
    zero_point_mag = "26.563"
    ):

    # Confirm zero_point_mag is correctly formatted.
    try:
        float(zero_point_mag)
    except TypeError:
        print("Please enter a valid Zero Point Magnitude.")


    # Remove old feedme file
    if os.path.isfile('galfit.feedme'):
        os.remove('galfit.feedme')

    # Create feedme file
    feedMe = open('galfit.feedme', 'w')

    raw_file_name = file_input.split("/")[-1]

    # Fitting parameters
    fit_params = {
        "A" : f"A) {file_input}",
        "B" : f"B) {file_output}gal_{raw_file_name}",
        "C" : f"C) {sigma_image}",
        "D" : f"D) {psf_image}",
        "E" : f"E) {psf_sample_factor}",
        "F" : f"F) {pixel_mask_image}",
        "G" : f"G) {constraints_file}",
        "H" : f"H) {image_region}",
        "I" : f"I) {convolution_box}",
        "J" : f"J) {zero_point_mag}",
        "K" : f"K) {plate_scale}",
        "O" : f"O) {display_type}",
        "P" : "P) 0"
        }

    # Object parameters
    object_params = {
        "0"  : f"0) {profile_type}",
        "1"  : f"1) {pos_xy}  1",
        "3"  : f"3) {integrated_mag}  1",
        "4"  : f"4) {eff_radius}  1",
        "5"  : f"5) {sersic_index}  1",
        "9"  : f"9) {axis_ratio}  1",
        "10" : f"10) {position_angle}  1",
        "Z"  : f"Z) 0"
        }

    # Write fitting paramters to file
    for key,val in fit_params.items():
        feedMe.write(val + "\n")

    feedMe.write("\n")

    # Write object paramters to file
    for key,val in object_params.items():
        feedMe.write(val + "\n")

    feedMe.close()

    # Run galfit with feedme file
    os.system("galfit galfit.feedme")
    os.system("rm galfit.*") # Remove galfit temporary files



# Make a cutout
def cutOut(
    inputFile, # .fits file to be cropped
    ra,
    dec,
    co_size = 101, # Default cutout size
    output_location = "cropped/", # Default relatrive to current working directory
    fits_layer = 0, # Set to whatever layer you wish to crop
    silent = False
    ):

    # total box side of cutout in pixels
    co_size = 101.0

    # open the file
    filename = inputFile
    bigdrc = fits.open(filename)

    # define the cutout's central coordinate
    # convert coordinate to pixels
    co_ra     = ra
    co_dec    = dec
    co_c      = SkyCoord(co_ra + ' ' + co_dec, unit=(u.hour, u.deg), frame = 'icrs')
    w         = wcs.WCS(bigdrc[0].header)
    co_x,co_y = w.wcs_world2pix(co_c.ra.deg, co_c.dec.deg, 1, ra_dec_order = True)
    root_file_name   = (str(inputFile.split("/")[-1])).split(".fits")
    co_name = (root_file_name[0] + "_" + str(co_x).split(".")[0] + "_" + str(co_y).split(".")[0] + ".fits")

    # cut out a stamp around that x and y
    # note: sourced how-to from https://github.com/astropy/photutils/issues/338
    position = (co_x, co_y)
    size     = (co_size, co_size)
    cutout   = Cutout2D(bigdrc[0].data, position, size, wcs=w)
    newdata  = cutout.data

    # make a new header
    # note: don't use "cutout.wcs.to_header()" b/c it loses all the HST-related information
    newhdr = bigdrc[0].header
    # update the WCS
    newhdr['CRPIX1'] = cutout.wcs.wcs.crpix[0]
    newhdr['CRPIX2'] = cutout.wcs.wcs.crpix[1]

    # Delete an existing file (or else fatal error), and write out
    try:
        os.system('rm ' + output_location + co_name)
    except:
        pass

    fits.writeto(output_location + co_name, newdata, newhdr)

    if not silent:
        print('Made: ' + co_name)


# Make multiple cut outs given a regions file
def bulkCutOut(
    input_file,
    regions_file, # A two column comma separated ascii file with no headers
    co_size = 101, # Default cutout size
    output_location = "cropped/", # Default relatrive to current working directory
    fits_layer = 0, # Set to whatever layer you wish to crop
    silent = False
    ):

    coords = pd.read_csv(regions_file, delimiter = ",", header = None, names = ["ra", "dec"])

    for index,row in coords.iterrows():
        cutOut(input_file, row["ra"], row["dec"], co_size, output_location, fits_layer, silent)





# Take a galfit output image and pull pertinent metrics from the headers. Then output as a list.
def _getMetrics(path, fileName, source):
    fitsFile = path + fileName
    hdulist = fits.open(fitsFile)

    imgName = fileName
    r_eff, r_err = hdulist[2].header['1_RE'].split(' +/- ')
    sersicIndex, sersicErr = hdulist[2].header['1_N'].split(' +/- ')
    axisRatio, axisErr = hdulist[2].header['1_AR'].split(' +/- ')
    positionAngle, positionErr = hdulist[2].header['1_PA'].split(' +/- ')

    results = [imgName,
               r_eff,
               r_err,
               sersicIndex,
               sersicErr,
               axisRatio,
               axisErr,
               positionAngle,
               positionErr,
               source]
    return results

def getPixelMask(input_file, output_location = "cropped/", img_index = 0, std_dev = 3):
    hdulist = fits.open(output_location + input_file)
    pixData = hdulist[img_index].data
    meanPix = pixData.mean()
    pixData[pixData < (meanPix + std_dev * pixData.std())] = 0
    pixData[pixData > (meanPix + std_dev * pixData.std())] = 1
    mask_file = 'mask_' + input_file

    if os.path.isfile(output_location + mask_file):
        os.remove(output_location + mask_file)
    mask = fits.PrimaryHDU(pixData)

    try:
        hdulist[img_index].writeto(output_location + mask_file)
        print("Made " + mask_file)
    except:
        print("Mask failed.")

def _createTable():
    # Create an empty table with appropriate headers
    new_table = QTable(names = ('Image',
                                    'Effective Radius',
                                    'r_eff Err',
                                    'Sersic Index',
                                    'Sersic Err',
                                    'Axis Ratio',
                                    'Axis Err',
                                    'Position Angle',
                                    'Position Err',
                                    'Source'),
                          dtype = ['U25', 'U25', 'U25','U25', 'U25', 'U25', 'U25', 'U25', 'U25', 'U25'])
    return new_table