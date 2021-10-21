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
from shapely.geometry.polygon import Polygon
import shapely.ops
import pyproj

# Own libraries
import georaster as geor


## Set up arguments ##
parser = argparse.ArgumentParser(description='Create a mosaic TDX DEM from tiles stored locally. The tiles are merged, making sure that no horizontal shift is introduced, filtered by removing HEM > 0.5 pixels and optionally setting the geoid to EGM96.')

# Positional arguments
parser.add_argument('dataDir', type=str, help='str, path to the folder containing all DEM files. The tree below contains [NS]??/[EW]???/TDM1_DEM__30_[NS]??[EW]???_V01_C/DEM/TDM1_DEM__30_[NS]??[EW]???_DEM.tif and [NS]??/[EW]???/TDM1_DEM__30_[NS]??[EW]???_V01_C/AUXFILES/TDM1_DEM__30_[NS]??[EW]???_HEM.tif.')
parser.add_argument('outfile', type=str, help='str, path to the output raster file')
    
# Optional arguments
parser.add_argument('-te', dest='extent', type=float, help='output DEM extent xmin ymin xmax ymax',nargs=4)
parser.add_argument('-tr', dest='tr', type=float, default=None, help='float, the output DEM resolution in t_srs units. Default is determined automatically by gdal from original tiles. WARNING: this can lead to strong resampling for high latitudes as TDX tiles are not square, but gdal assumes square tiles.')
parser.add_argument('-t_srs', dest='srs', type=str, default=None, help='str, a PROJ-4 string of the output projection. Default is same as input tiles, i.e. EPSG 4326.')
parser.add_argument('-r', dest='resampling', type=str, default='cubic', help='resampling algorithm to use, check GDAL for a list of algorithms available (Default is cubic)')
parser.add_argument('-ot', dest='ot', type=str, default='Int16', help='output data type. Original is Float32 but default is saved as Int16 to save space.')
parser.add_argument('-geoid', dest='geoid', action='store_true', help='if set, will express elevation with reference to the EGM96 geoid. TDX DEM is w.r.t. the WGS84 ellipsoid.')
parser.add_argument('-co', dest='co', type=str, default=['COMPRESS=LZW','TILED=YES','BLOCKXSIZE=256','BLOCKYSIZE=256','BIGTIFF=IF_SAFER'], help="GDAL creation options (Default is COMPRESS=LZW TILED=YES BLOCKXSIZE=256 BLOCKYSIZE=256 BIGTIFF=IF_SAFER", nargs='*')
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
if args.geoid:
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

## Select tiles overlapping with requested extent ##

# Bounding box of output
xmin, ymin, xmax, ymax = args.extent
bbox = Polygon([(xmin, ymin), (xmax, ymin),(xmax, ymax), (xmin, ymax)])

# Lat/lon projection
latlon_crs = pyproj.crs.CRS('epsg:4326')

if args.srs is not None:

    # Output CRS
    out_crs = pyproj.crs.CRS(args.srs)

    # Instance to reproject from lat/lon
    reproj = pyproj.Transformer.from_crs(latlon_crs, out_crs, always_xy=True, skip_equivalent=True).transform

else:
    reproj = pyproj.Transformer.from_crs(latlon_crs, latlon_crs, always_xy=True, skip_equivalent=True).transform

# Select tiles
tiles2keep = []
for k in range(len(DEMfiles)):

    # Convert tile bounds to output CRS
    tile = Polygon([(lons[k], lats[k]), (lons[k]+1, lats[k]), (lons[k]+1, lats[k]+1), (lons[k], lats[k]+1)])
    tile_proj = shapely.ops.transform(reproj, tile)

    # Keep if intersect with output bounds
    if tile_proj.intersects(bbox):
        tiles2keep.append(k)

DEMfiles = DEMfiles[tiles2keep]
print("Found %i tiles overlapping with set extent" %len(DEMfiles))

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
# if no resampling method specified, can lead to N/S oriented banding artifacts when tiles above 70N (with different grid cell) are included
print("* Merge tiles *")
vrtfile = outDir + 'merged_DEM.vrt'
cmd = "gdalbuildvrt -r %s %s -input_file_list %s -resolution highest" %(args.resampling, vrtfile, listfile)
cprint(cmd); check_call(cmd.split())

## Create final mosaic ##
# Must use gdalwarp, as option -te in gdalbuildvrt creates horizontal shift...
print("* Create final mosaic *")

mosaic_file = outDir + 'merged_dem.tif'

# gdalwarp options
gdal_opts = ''
if args.tr is not None:
    gdal_opts += " -tr %g %g" %(args.tr, args.tr)
if args.srs is not None:
    gdal_opts += " -t_srs '%s'" %args.srs
if len(args.co)>0:
    gdal_opts += ' -co ' + ' -co '.join(args.co)

# gdalwarp command
cmd = "gdalwarp -te %f %f %f %f -r %s %s %s -ot %s %s" %(xmin, ymin, xmax, ymax, args.resampling, vrtfile, mosaic_file, args.ot, gdal_opts)
cprint(cmd); check_call(cmd, shell=True)

# Correct geoid
if args.geoid:
    print("* Convert to geoid EGM96 *")
    outPrefix = outDir + 'merged_DEM'
    corrDEMf = outPrefix + '-adj.tif'
    cmd = "dem_geoid --geoid EGM96 %s -o %s" %(mosaic_file,outPrefix)
    cprint(cmd); check_call(cmd, shell=True)
else:
    corrDEMf = mosaic_file

# Mv to final DEM
print("* Move and delete temporary folder *")
cmd = "mv %s %s" %(corrDEMf,args.outfile)
cprint(cmd); check_call(cmd, shell=True)

# Delete tmp folder
cmd = "rm -r %s" %outDir
cprint(cmd); check_call(cmd, shell=True)
