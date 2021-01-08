# DEM_tools
Contains various tools to download and process global Digital Elevation Model datasets (SRTM, TanDEM-X...)

### DEM sources ###


#### ArcticDEM ####

Downloaded from https://www.pgc.umn.edu/data/arcticdem/. \
Vertical reference: WGS84 ellipsoid. \
Projection: Polar stereographic \

Script: `arcticDEM_downloader.py`

#### TanDEM-X 90 m ####

Infos: https://tandemx-90m.dlr.de/ \
Vertical reference: WGS84 ellipsoid \
Projection: lat/lon \

Script to create DEM mosaics from locally downloaded tiles: `create_tdm_mosaic.py`

#### Copernicus GLO-30 ####

Infos: https://spacedata.copernicus.eu/web/cscda/dataset-details?articleId=394198\
Vertical reference: EGM2008 geoid \
Projection: lat/lon \

Script to create DEM mosaics from locally downloaded tiles: `create_GLODEM_mosaic.py`


#### SRTM ####

Downloaded from OpenTopography (?) => link seem deprecated \
Vertical reference: unsure \
Projection: lat/lon

Script: `srtm_downloader.py`=> Deprecated
