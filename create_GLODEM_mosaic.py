#!/usr/bin/env python
#coding=utf-8

"""
Description : Create a mosaic of the copernicus 30m GLODEM from tiles stored locally. The tiles are merged, making sure that no horizontal shift is introduced.
GLO-30 DEM vertical reference is EGM2008. With the -ellips option, the DEM is converted to the WGS84 ellispoid for comparison with other DEMs.

Author : Amaury Dehecq
Last modified : Jan 2021
"""

# Python libraries
from glob import glob
import numpy as np
import os, argparse
from subprocess import call, check_call
from shapely.geometry.polygon import Polygon
import shapely.ops
import pyproj

## Set up arguments ##
parser = argparse.ArgumentParser(description='Create a mosaic of the copernicus 30m GLODEM from tiles stored locally. The tiles are merged, making sure that no horizontal shift is introduced.')

# Positional arguments
parser.add_argument('dataDir', type=str, help='str, path to the folder containing all DEM files. The tree below contains Copernicus_DSM_10_[NS]??_00_[EW]???_00/DEM/*_DEM.tif.')
parser.add_argument('outfile', type=str, help='str, path to the output raster file')
    
# Optional arguments
parser.add_argument('-te', dest='extent', type=float, help='output DEM extent xmin ymin xmax ymax',nargs=4)
parser.add_argument('-tr', dest='tr', type=float, default=None, help='float, the output DEM resolution in t_srs units. Default is determined automatically by gdal from original tiles. WARNING: this can lead to strong resampling for high latitudes as TDX tiles are not square, but gdal assumes square tiles.')
parser.add_argument('-t_srs', dest='srs', type=str, default=None, help='str, a PROJ-4 string of the output projection. Default is same as input tiles, i.e. EPSG 4326.')
parser.add_argument('-r', dest='resampling', type=str, default='cubic', help='resampling algorithm to use, check GDAL for a list of algorithms available (Default is cubic)')
parser.add_argument('-ot', dest='ot', type=str, default='Int16', help='output data type. Original is Float32 but default is saved as Int16 to save space.')
parser.add_argument('-ellips', dest='ellips', action='store_true', help='if set, will express elevation with reference to the ellispoid. GLO-30 DEM is w.r.t. the EGM2008 geoid.')
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
if args.ellips:
    try:
        FNULL = open(os.devnull, 'w')
        call('dem_geoid', stdout=FNULL)
    except FileNotFoundError:
        raise SystemExit("ERROR: ASP's dem_geoid must be in your path")


# Find all DEM tiles
DEMfiles = np.sort(glob(args.dataDir + '/Copernicus_DSM_10_[NS]??_00_[EW]???_00/DEM/*_DEM.tif'))

# Create temporary output directory
outDir = os.path.dirname(args.outfile)
if outDir=='':  # in case output in current directory
    outDir = './'


## Find lat and lons coordinates of the tiles ##
lats, lons = [], []
for f in DEMfiles:
    basef = os.path.basename(f)
    items = basef.split('_')
    lat = items[3]
    lon = items[5]
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

# Save to text file
listfile = outDir + '/list_DEMs.txt'
np.savetxt(listfile, DEMfiles, fmt='%s') 


# Create merged VRT file
# if no resampling method specified, can lead to N/S oriented banding artifacts when tiles above 70N (with different grid cell) are included
print("* Merge tiles *")
vrtfile = os.path.splitext(args.outfile)[0] + '.vrt'
cmd = "gdalbuildvrt -r %s %s -input_file_list %s -resolution highest" %(args.resampling, vrtfile, listfile)
cprint(cmd); check_call(cmd.split())

## Create final mosaic ##
# Must use gdalwarp, as option -te in gdalbuildvrt creates horizontal shift...
print("* Create final mosaic *")

# gdalwarp options
gdal_opts = ''
if args.tr is not None:
    gdal_opts += " -tr %g %g" %(args.tr, args.tr)
if args.srs is not None:
    gdal_opts += " -t_srs '%s'" %args.srs
if len(args.co)>0:
    gdal_opts += ' -co ' + ' -co '.join(args.co)

# gdalwarp command
cmd = "gdalwarp -te %f %f %f %f -r %s %s %s -ot %s %s" %(xmin, ymin, xmax, ymax, args.resampling, vrtfile, args.outfile, args.ot, gdal_opts)
if args.overwrite: cmd += ' -overwrite'
cprint(cmd); check_call(cmd, shell=True)

# Correct geoid
if args.ellips:
    print("* Convert to WGS84 ellispoid *")
    outPrefix = os.path.splitext(args.outfile)[0]
    corrDEMf = outPrefix + '-adj.tif'
    # Set nodata value to 0, which masks oceans
    cmd = "dem_geoid --geoid EGM2008 --reverse-adjustment %s -o %s --nodata_value 0" %(args.outfile,outPrefix)
    cprint(cmd); check_call(cmd, shell=True)

    # Copy itermediate file, temporarily
    cmd = "cp %s %s_egm2008.tif" %(args.outfile, outPrefix)
    cprint(cmd); check_call(cmd, shell=True)

    # Convert back to Int16, also set nodata to -32768 as 0 could be an actual value (e.g. cause issues with gdal_fillnodata.py)
    if args.ot != 'Float32':
        cmd = "gdalwarp -srcnodata 0 -dstnodata -32768 -ot %s -co COMPRESS=LZW -co TILED=YES -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 -co BIGTIFF=IF_SAFER %s %s -overwrite" %(args.ot, corrDEMf, args.outfile)
        cprint(cmd); check_call(cmd, shell=True)
        os.remove(corrDEMf)
    else:
        cmd = "gdalwarp -srcnodata 0 -dstnodata -32768 -co COMPRESS=LZW -co TILED=YES -co BLOCKXSIZE=256 -co BLOCKYSIZE=256 -co BIGTIFF=IF_SAFER %s %s -overwrite" %(corrDEMf, args.outfile)
        cprint(cmd); check_call(cmd, shell=True)
        os.remove(corrDEMf)


# Remove temporary files
os.remove(vrtfile)
os.remove(listfile)
