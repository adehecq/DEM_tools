#!/usr/bin/env python
#coding=utf-8

"""
Description : Tool to bulk download the ArticDEM over a specified region.

Author : Amaury Dehecq
Last modified : Jul 2017
"""

import argparse, os, sys
import tempfile
import subprocess
from glob import glob
import numpy as np
from osgeo import gdal, ogr
from geoutils import geovector as vect
import georaster as raster

#Set up arguments
parser = argparse.ArgumentParser(description='Tool to bulk download the ArticDEM over a specified region.')


# Positional arguments
parser.add_argument('tiles_file', type=str, help='str, path to the shapefile containing the outlines of the ArcticDEM tiles')
parser.add_argument('outfile', type=str, help='str, path to the output file')


# Optional arguments
parser.add_argument('-shp', dest='shp', type=str, help='str, path to a shapefile containing RGI outlines. Only the tiles interesting with the glaciers will be downloaded (Default is None, but this or -te must be specified).', default=None)
parser.add_argument('-area_th', dest='area_th', type=str, help='float, glacier with area (as read from RGI attributes in km2) below this threshold will be excluded (Default is 5 km2).', default=5)
parser.add_argument('-te', dest='te', type=str, help='extent (lonmin, lonmax, latmin, latmax) of the output DEM.', nargs=4, default=None)
parser.add_argument('-tr', dest='tr', type=str, help='resolution (xres, yres) of the output DEM (Default is from the input DEMs).', nargs=2, default=None)
parser.add_argument('-t_srs', dest='t_srs', type=str, help='projection of the output DEM in PROJ4 format (Default is from the input DEMs).', default=None)
parser.add_argument('-overwrite', dest='overwrite', help='if set, will overwrite the output file if exists.', action='store_true')
parser.add_argument('-outdir', dest='outdir', type=str, help='str, path to output directory where to save the tiles (Default, create a temporary file that is deleted a the end).', default=None)

args = parser.parse_args()

## Hard-coded server parameters ##

# URL of the ArctiDEM server
URL = 'http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/'

# Version to be downloaded
version = 'v2.0'

## Check that all necessary input arguments are specified
if (args.shp==None) & (args.te==None):
      print "ERROR: At least one of -shp or -te must be specified"
      sys.exit(1)



## Functions to convert MultiPolygons to Polygons
def multipoly2poly(in_lyr, out_lyr):
      k=0
      for in_feat in in_lyr:
            geom = in_feat.GetGeometryRef()
            if geom.GetGeometryName() == 'MULTIPOLYGON':
                  geom_part=geom.GetGeometryRef(0)
                  addPolygon(geom_part.ExportToWkb(), out_lyr)
            else:
                  addPolygon(geom.ExportToWkb(), out_lyr)
            k+=1

      print "Processed %i features" %k

      
def addPolygon(simplePolygon, out_lyr):
      featureDefn = out_lyr.GetLayerDefn()
      polygon = ogr.CreateGeometryFromWkb(simplePolygon)
      out_feat = ogr.Feature(featureDefn)
      out_feat.SetGeometry(polygon)
      out_lyr.CreateFeature(out_feat)


## Read ArticDEM tiles shapefile ##
tiles = vect.SingleLayerVector(args.tiles_file)

if args.te!=None:
      lonmin, lonmax, latmin, latmax = args.te
      tiles.crop(lonmin,lonmax,latmin,latmax,latlon=True)

tiles.read()


## Read and simplify the external shapefile (RGI for now)

if args.shp!=None:
      
      print "*** Generate simplified geometry for shapefile ***"

      ## Read glacier outlines ##
      gl_outlines = vect.SingleLayerVector(args.shp)
      gl_outlines = gl_outlines.reproject(tiles.srs)
      gl_outlines.layer.SetAttributeFilter('Area>%s' %args.area_th)
      gl_outlines.read()
      #gl_outlines.update_extent()  # causes layer.GetNextFeature() to not return anything

      ## Replace MultiPolygons by Polygons ##
      nfeat = gl_outlines.FeatureCount()  # somehow necessary for the loop on features to work
      print "%i features to process" %nfeat
      out_ds = ogr.GetDriverByName('Memory').CreateDataSource('')
      out_lyr = out_ds.CreateLayer('poly', srs=gl_outlines.srs)
      multipoly2poly(gl_outlines.layer, out_lyr)
      
      ## Create a single geometry for simplicity ##
      union = ogr.Geometry(ogr.wkbMultiPolygon)
      for feat in out_lyr:
            union.AddGeometry(feat.GetGeometryRef())
      union=union.Simplify(0)

      # Check that geometry is valid (otherwise will fail later)
      if not union.IsValid():
            print "ERROR with geometry"
            sys.exit(1)


      ## Find polygons that intersect with glacier outlines (union)
      
      print "*** Find tiles intersecting with input outlines ***"
      inds = []
      for k in xrange(tiles.FeatureCount()):
            gdal.TermProgress_nocb(float(k)/tiles.FeatureCount())
    
            feat = tiles.features[k]
            if union.Intersect(feat.GetGeometryRef()):
                  inds.append(k)

      list_tiles = np.sort(np.unique(tiles.fields.values['tile'][inds]))

else:
      list_tiles = np.sort(np.unique(tiles.fields.values['tile']))


#list_tiles = ['32_34', '32_35', '32_36', '31_34','31_35', '31_36']



## Create output directory ##
if args.outdir!=None:
      if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
      outdir = args.outdir
else:
      # Generate a temporary dir to save downloaded files
      outdir = tempfile.mkdtemp(dir='/tmp')

## Download files ##

print "\n*** Download tiles in %s ***" %outdir
print "%i tiles to download" %len(list_tiles)

for t in list_tiles:

      print "\n ** Tile %s **" %t
      
      # wget options:
      # -r for recursive download
      # -N replace file locally if timestamp is older
      # -nd to save everything in one folder, remove original tree
      # -np does not include parent directories
      # -nv non verbose, display less information
      # -R to exclude files
      # -P destination folder
      # Don't forget the slash at the end of URL!
      wget_cmd = ['wget','-r','-N','-nd','-np','-nv','-R','index.html*','-R','robots.txt*','%s/%s/%s/' %(URL,version,t), '-P', '%s' %outdir]
      print ' CMD = ' + ' '.join(wget_cmd)
      out=subprocess.call(wget_cmd)
      if out!=0:
            print "Error in download !!"
            continue
      

## Untar all files ##

print "\n*** Extract archives ***"

# Get list of tar files to extract
tar_files = []
for t in list_tiles:
    l=glob(outdir + '/%s_?_?_5m_%s.tar*' %(t,version))
    tar_files.extend(l)

# From each archive, extract only the dem file
dem_files = []
for f in tar_files:
      if f[-7:]=='.tar.gz':
            print 'tar.gz', f 
            dem_file = f.replace('.tar.gz','_reg_dem.tif')
            cmd = ['tar','xzvf', '%s' %f,'-C',outdir,os.path.basename(dem_file)]
      else:
            print 'tar', f
            dem_file = f.replace('.tar','_reg_dem.tif')
            cmd = ['tar','xvf', '%s' %f,'-C',outdir,os.path.basename(dem_file)]
      print ' '.join(cmd)
      out=subprocess.call(cmd)
      if out!=0:
            print "Error extracting file %s !!" %f
            continue
      dem_files.append(dem_file)

# Save the list of dem files to file
list_file = outdir + '/list_files.txt'
np.savetxt(list_file,dem_files,fmt='%s')

# try:      
#       cmd = "for f in `ls %s/*.tar.gz`; do echo $f; tar xzvf $f -C %s -q '*dem.tif'; done" %(outdir,outdir)
#       print cmd; out=os.system(cmd)
#       if out!=0:
#             sys.exit()
# except SystemExit:
#       cmd = "for f in `ls %s/*.tar.gz`; do echo $f; tar xzvf $f -C %s --wildcards '*dem.tif'; done" %(outdir,outdir)
#       print cmd; out=os.system(cmd)
#       if out!=0:
#             print "Error in archive"
#             sys.exit()
      


## Generate output DEM ##

print "\n*** Generate final DEM ***"

#cmd = 'gdalwarp  -r average -co "COMPRESS=LZW" --optfile %s %s' %(list_file,args.outfile) # compression will fail if file > 4GB or BIGTIFF=Yes must be used
cmd = 'gdalwarp  -r average -ot Int16 --optfile %s %s' %(list_file,args.outfile)  

if args.te!=None:

      # reproject extent to DEM extent
      img = raster.SingleBandRaster(dem_files[0],load_data=False)
      x1, y1 = img.proj(lonmin,latmin)
      x2, y2 = img.proj(lonmin,latmax)
      x3, y3 = img.proj(lonmax,latmax)
      x4, y4 = img.proj(lonmax,latmin)
      xmin = min(x1,x2,x3,x4)
      xmax = max(x1,x2,x3,x4)
      ymin = min(y1,y2,y3,y4)
      ymax = max(y1,y2,y3,y4)
      cmd+= ' -te %f %f %f %f' %(xmin, ymin, xmax, ymax) 

if args.tr!=None:
      cmd+= ' -tr %s' %' '.join(args.tr)

if args.t_srs!=None:
      cmd+= ' -ts_srs %s' %args.t_srs

if args.overwrite!=False:
      cmd+= ' -overwrite'

print cmd; out=subprocess.call(cmd,shell=True)
if out!=0:
      print "Error reprojecting file"
      sys.exit(1)


## Remove temporary directory/files ##

if args.outdir==None:
      os.system('rm -r %s' %outdir)
else:
      for f in dem_files:
            os.remove(f)
      #os.remove(list_file)
      os.system('rm -r %s/robots.txt*' %outdir)
