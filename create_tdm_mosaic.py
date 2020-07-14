#!/usr/bin/env python
#coding=utf-8

"""
Description : Create a mosaic TDX DEM from tiles stored locally. The tiles are merged, making sure that no horizontal shift is introduced, filtered by removing HEM > 0.5 pixels and setting the geoid to EGM96.

Author : Amaury Dehecq
Last modified : July 2020
"""

# Python libraries
from glob import glob
import numpy as np
import os, argparse
from subprocess import call, check_call

# Own libraries
import georaster as geor


## Set up arguments ##
parser = argparse.ArgumentParser(description='Create a mosaic TDX DEM from tiles stored locally. The tiles are merged, making sure that no horizontal shift is introduced, filtered by removing HEM > 0.5 pixels and setting the geoid to EGM96.')

# Positional arguments
parser.add_argument('dataDir', type=str, help='str, path to the folder containing all DEM files. The tree below contains [NS]??/[EW]???/TDM1_DEM__30_[NS]??[EW]???_V01_C/DEM/TDM1_DEM__30_[NS]??[EW]???_DEM.tif and [NS]??/[EW]???/TDM1_DEM__30_[NS]??[EW]???_V01_C/AUXFILES/TDM1_DEM__30_[NS]??[EW]???_HEM.tif.')
parser.add_argument('outfile', type=str, help='str, path to the output raster file')
    
# Optional arguments
parser.add_argument('-te', dest='extent', type=float, help='output DEM extent xmin ymin xmax ymax',nargs=4)
parser.add_argument('-r', dest='resampling', type=str, default='bilinear', help='resampling algorithm to use, check GDAL for a list of algorithms available (Default is bilibear)')
parser.add_argument('-overwrite', dest='overwrite', action='store_true', help='if set, will overwrite output file')
    
args = parser.parse_args()

# Colored print
def cprint(cmd):
    GREEN = '\033[32m'
    ENDC = '\033[0m'
    print(GREEN + cmd + ENDC)


## Sanity checks ##

# Check that output file doesn't already exist
if (not args.overwrite) & (os.path.exists(args.outfile)):
    raise SystemExit('ERROR: Output file already exists')

# Check that extent is not empty
if args.extent is None:
    raise SystemExit('ERROR: Extent must be specified')

# Check that ASP is loaded
try:
    FNULL = open(os.devnull, 'w')
    call('dem_geoid', stdout=FNULL)
except FileNotFoundError:
    raise SystemExit("ERROR: ASP's dem_geoid must be in your path")


# Find all DEM tiles
DEMfiles = np.sort(glob(args.dataDir + '/[NS]??/[EW]???/TDM1_DEM__*/DEM/TDM1_DEM__*_DEM.tif'))

# Create temporary output directory
outDir = os.path.dirname(args.outfile)
if outDir=='':  # in case output in current directory
    outDir = '.'

outDir += '/tmp/'
if os.path.exists(outDir):
    raise SystemExit('ERROR: Temporary output folder %s already exists' %outDir)
else:
    print("* Create temporary folder %s *" %outDir)
    os.makedirs(outDir)

## Find lat and lons coordinates of the tiles ##
lats, lons = [], []
for f in DEMfiles:
    basef = os.path.basename(f)
    tileID = basef.split('_')[4]
    lat = tileID[:3]
    lon = tileID[3:]
    if lat[0] == 'N':
        lat = int(lat[1:])
    elif lat[0] == 'S':
        lat = -int(lat[1:])
    else:
        print("File not understood")
    if lon[0] == 'E':
        lon = int(lon[1:])
    elif lon[0] == 'W':
        lon = -int(lon[1:])
    lats.append(lat)
    lons.append(lon)
lons = np.array(lons)
lats = np.array(lats)

# Select tiles overlapping with requested extent

xmin, ymin, xmax, ymax = args.extent

tiles2keep = []
for k in range(len(DEMfiles)):
    if (lons[k] >= xmin) & (lons[k] <= xmax+2):  #+2 because some tiles are 2 degrees wide
        if (lats[k] >= ymin) & (lats[k] < ymax):
            tiles2keep.append(k)

DEMfiles = DEMfiles[tiles2keep]


## First mask pixels with HEM > 0.5 in TDX DEM ##

print("* Mask TDX tiles *")
masked_DEM_files = []

for f in DEMfiles:

    # Filename info
    basef = os.path.basename(f)
    prefix, ext = os.path.splitext(basef)
    outPrefix = outDir + prefix

    # Open DEM and HEM files
    parent_folder = os.path.dirname(os.path.dirname(f))
    HEMfolder = parent_folder + '/AUXFILES/'
    HEMf = HEMfolder + basef.replace('DEM.tif','HEM.tif')
    hem = geor.SingleBandRaster(HEMf)
    dem = geor.SingleBandRaster(f)

    # Apply mask
    nodata = -32767
    dem.r[hem.r>0.5] = nodata

    # Save - Note: GTiff attribute 'AREA_OR_POINT=Point' is replaced by 'Area' removing issue with ASP geodiff
    maskedDEMf = outPrefix + '_masked.tif'
    dem.save_geotiff(maskedDEMf, dem.ds.GetRasterBand(1).DataType, nodata_value=nodata)
    masked_DEM_files.append(maskedDEMf)
    print("Saved file %s" %maskedDEMf)

# Save to text file
listfile = outDir + '/list_DEMs.txt'
np.savetxt(listfile,masked_DEM_files,fmt='%s') 


# Create merged VRT file
print("* Merge tiles *")
vrtfile = outDir + 'merged_DEM.vrt'
cmd = "gdalbuildvrt -r %s %s -input_file_list %s" %(args.resampling, vrtfile, listfile)
cprint(cmd); check_call(cmd.split())

# Crop mosaic to final extent
# Must use gdalwarp, as option -te in gdalbuildvrt creates horizontal shift...
print("* Crop DEM *")
mosaic_file = outDir + 'merged_dem.tif'
cmd = "gdalwarp -te %g %g %g %g -r %s %s %s" %(xmin, ymin, xmax, ymax, args.resampling, vrtfile, mosaic_file)
cprint(cmd); check_call(cmd.split())

# Correct geoid
print("* Convert to geoid EGM96 *")
outPrefix = outDir + 'merged_DEM'
corrDEMf = outPrefix + '-adj.tif'
cmd = "dem_geoid --geoid EGM96 %s -o %s" %(mosaic_file,outPrefix)
cprint(cmd); check_call(cmd.split())

# Mv to final DEM
print("* Move and delete temporary folder *")
cmd = "mv %s %s" %(corrDEMf,args.outfile)
cprint(cmd); check_call(cmd.split())

# Delete tmp folder
cmd = "rm -r %s" %outDir
cprint(cmd); check_call(cmd.split())
