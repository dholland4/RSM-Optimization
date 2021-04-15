# Created by Darren Holland
# Modified by Darren Holland 2020-11-02
#**********************************************************************
# File for creating the mask geometry
#**********************************************************************
# Input: python3 Geoinp.py $1 $2 $3 $4 $5 $6 $7 $8 $9 $10  
# $1 is the fin thickness
# $2 is the wall thickness
# $3 is the theta discretization
# $4 is the phi discretization
# $5 is the fin width
# $6 is the wall width
# $7 is the starting phi
# $8 is the ending phi
# $9 is the detector radius
# $10 is the detector (half) height

# Output is a design's node and element mapping in Abaqus format

# Load files to create matrix and node/element mapping
import CreateGeo,CreateMat,sys

# Input/output files
Userfile = 'EigMat.txt'     # Design matrix
MainDir=''                  # Directory to save files
dir_settings='Settings.cc'  # Geant settings file
Nodesfile='/nodes.inp'      # Output Node file name
Elemfile='/elem.inp'        # Output Element file name

# Geometry settings in cm (overwritten by SurrWrapper.sh or StartDakota.sh)
det_height=7.62             # Detector height
det_rad=3.81                # Detector radius
MaskMinThick=0.1            # Encasement thickness
sleeve_inner_rad=3.81       # Sleeve inner radius
sleeve_outer_rad=4.1275     # Sleeve outer radius
sleeve_height=5.0           # Sleeve top
sleeve_bottom=55            # Sleeve bottom
s_dist=86.36                # Source distance

start_cells=700

# Depreciated setting
RSMmaxsize=20

# Create matrix for geometry according to method in optimization paper.  
# If invalid geometry returns "1"
Bad_Geo = CreateMat.CreateMat(Userfile,float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),\
int(sys.argv[5]),int(sys.argv[6]),float(sys.argv[7]),float(sys.argv[8]),MaskMinThick)
print(Bad_Geo)

# If valid geometry then continue with the mask creation
if Bad_Geo == 0: 
    CreateGeo.CreateGeo(RSMmaxsize,float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[7]),float(sys.argv[8]),\
    det_height,det_rad,sleeve_inner_rad,sleeve_outer_rad,sleeve_height,sleeve_bottom,\
	start_cells,s_dist,Userfile,MainDir,Nodesfile,Elemfile,dir_settings,float(sys.argv[9]),float(sys.argv[10]))