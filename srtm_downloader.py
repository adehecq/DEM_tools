#! /usr/bin/env python
#coding=utf-8
###############################################################################
#  scihub_opensearch.py
#
#  Purpose:  Tool to search, download and merge SRTM DEM tiles for a specific region
#
#  Author:   Amaury Dehecq
#  Created:  Mar 2016
#
###############################################################################

import argparse, sys, os
import numpy as np
from glob import glob

## add ASTER with link as follows :
## http://earthexplorer.usgs.gov/download/4220/ASTGDEMV2_0N33E078/STANDARD/EE
##
## SRTM
## http://earthexplorer.usgs.gov/download/8360/SRTM1N33E076V3/GEOTIFF/EE/
##
## Landsat
## http://earthexplorer.usgs.gov/download/3119/LT52030342010316MPS00/STANDARD/EE
##
## can use LANDSAT-Download/download_landsat_scene.downloadChunks for the last two, but not for ASTER since one has to agree with the term and conditions...


SERVER="http://e4ftl01.cr.usgs.gov//MODV6_Dal_D/SRTM/"  # found using http://reverb.echo.nasa.gov/reverb/
base_dirs = {
'srtm1':'SRTMGL1',
'srtm3':'SRTMGL3'
}
version='.003'
directory='2000.02.11/'
temp_folder='./tmp'

def srtm_downloader2(west,east,south,north,dem_type,outname):
    
    # download the data
    cmd="wget -c --no-check-certificate -O %s 'http://ot-data1.sdsc.edu:9090/otr/getdem?north=%f&south=%f&east=%f&west=%f&demtype=%s'" %(outname,north,south,east,west,base_dirs[dem_type])
    print(cmd); os.system(cmd)

    
def srtm_downloader(west,east,south,north,dem_type,outname,tile_size=5):

    # round up to 1°x1°
    westi = int(np.floor(west))
    southi = int(np.floor(south))
    easti = int(np.ceil(east))
    northi = int(np.ceil(north))

    # bounding box
    lats = np.arange(south,north,tile_size)
    lons = np.arange(west,east,tile_size)
    ntot = len(lats)*len(lons)
    
    print("\n## Downloading data for bounding box : ##")
    print("## %i, %i, %i, %i ##" %(west, east, south, north))
    print("## %i files to download ##" %ntot)

    lons = np.append(lons,east)
    lats = np.append(lats,north)
    
    # create a temporary directory where to download the data
    if os.path.exists(temp_folder):
        print("ERROR : ./tmp folder already exist")
        sys.exit(1)
    else:
        os.mkdir(temp_folder)
    os.chdir(temp_folder)

    # Loop on all tiles
    for i in range(len(lats)-1):

        lat1 = lats[i]
        lat2 = lats[i+1]
        
        # latitude string
        if lat1<0:
            lat_code = "S%02i" %np.abs(lat1)
        else:
            lat_code = "N%02i" %lat1

        for j in range(len(lons)-1):

            lon1 = lons[j]
            lon2 = lons[j+1]
            
            # longitude string
            if lon1<0:
                lon_code = "W%03i" %np.abs(lon1)
            elif lon1>180:
                lon_code = "W%03i" %np.abs(lon1-360)
            else:
                lon_code = "E%03i" %lon1
                
            print("\n## Downloading tile %s%s ##" %(lat_code,lon_code))
            tile_file="%s%s.%s.hgt" %(lat_code,lon_code,base_dirs[dem_type])

            # Download tile
            srtm_downloader2(lon1,lon2,lat1,lat2,dem_type,tile_file)

    # merge together the tile
    # without specifying the extent, the DEM extent should be about right, at the nearest pixel
    # setting -te or -ul_lr will shift the tiles wrt each other
    print("\n## Merge tiles together ##\n")
    cmd="gdal_merge.py -o ../%s *.hgt" %(outname)  
    #cmd="gdalwarp *.hgt ../%s -te %g %g %g %g" %(outname,west,south,east,north)
    print(cmd); os.system(cmd)
    
    # remove temporary files
    print("\n## Remove temporary files ##")
    os.chdir('..')
    #os.system('rm -rf %s' %temp_folder)
        

#### Set up arguments ####

parser = argparse.ArgumentParser(description='Tool to search, download and merge SRTM DEM tiles for a specific region')

parser.add_argument('west', type=float, help='float, bounding box west coordinates')
parser.add_argument('east', type=float, help='float, bounding box east coordinates')
parser.add_argument('south', type=float, help='float, bounding box south coordinates')
parser.add_argument('north', type=float, help='float, bounding box north coordinates')
parser.add_argument('outname', type=str, help='str, output file name')
parser.add_argument('dem_type', type=str, help='str, DEM to download (srtm1/srtm3)')

args = parser.parse_args()

# check that output file does not already exist
if os.path.exists(args.outname):
    print("ERROR : Output file %s already exists\n" %args.outname)
    sys.exit(1)

# check that coordinates are in the right order
if args.south>args.north:
    print("ERROR : south > north")
    sys.exit(1)
if args.west>args.east:
    print("ERROR : west > east")
    sys.exit(1)


srtm_downloader(args.west,args.east,args.south,args.north,args.dem_type,args.outname)
