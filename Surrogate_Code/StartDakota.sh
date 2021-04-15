#!/bin/bash
# Created by Darren Holland
# Modified by Darren Holland 2020-11-02
#**********************************************************************
# This code starts the Dakota optimization.  The desired parameters and 
# settings are passed into Dakota (written to MOGAMAStemplate.in), which 
# calls the Geant wrapper (GeantWrapper.sh).  The wrapper creates the
# geometry (CreateMat.py and CreateGeo.inp), compiles the Geant evaluation
# code, runs Geant, and calls the analysis script (Analyze.py). The analysis
# code evaluates the objective functions based on the Geant-produced 
# spectrum and returns the values to Dakota.
#**********************************************************************
# Output: Dakota optimization information including summary files
#**********************************************************************

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Running Code
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ./StartDakota $1
# where $1 is the folder/project name

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Settings (Other settings available in SurrWrapper.sh)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#max_iterations=1	# Maximum number of iteration to run
min_thick=0.81		# Minimum mask thickness (cm)
min_mask=3.0		# Minimum wall/fin thickness
max_mask=12.		# Maximum wall/fin thickness (cm)
n_obj=4			    # Number of MOGA objective functions
StartSource=35		# Source phi start limit
nnodes=1			# Number of nodes to run Dakota
numproc=1           # Number of processors to run Dakota
# Concurrancy > 1 will not work since it only uses 1 node (have to use mpiexec instead, but that only runs on 1 or 2 processors!)
#concurrancy=$nnodes		# Number of designs/distance to run in Dakota simultaneously
concurrancy=1		# Number of designs/distance to run in Dakota simultaneously

#**********************************************************************
# Original Spartan Geometry Settings
# Limit Dakota parameter settings to this geometry to evaluate this design
#**********************************************************************
#deltatheta=10		# Theta discretization for source position 
#deltaphi=10		# Phi discretization for source position
#finthick=5.0
#wallthick=7.5
#finwidth=2
#wallwidth=1
#material=0
#phifinal=170           # Phi final limit (Source position and geometry creation)
			# Must be <= 180 - deltaphi if set manually (Dakota changes deltaphi)
			# If using uncomment later phifinal sed command 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create Dakota script
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
homedir=$(pwd)
dakotafile="$1DakotaMOGAMAS.in"
cp "MOGAMAStemplate.in" $dakotafile

# Apply setting changes to Geant wrapper template
cp $homedir/SurrWrapper.sh $homedir/SurrEval.sh
sed -i -e  "s?StartSource=1?StartSource=$StartSource?g" "$homedir/SurrEval.sh"
sed -i -e  "s?MaskMinThick=0.4053?MaskMinThick=$min_thick?g" "$homedir/SurrEval.sh"
sed -i -e  "s?numproc=2?numproc=$numproc?g" "$homedir/SurrEval.sh"

# Set the number of threads for running Geant in parallel
geantthreads=$(awk -v dd=$nnodes -v ee=$numproc -v ff=$concurrancy 'BEGIN {printf "%d\n",int(dd*ee/ff)}');	
sed -i -e  "s?numproc=2?numproc=$numproc?g" "$homedir/SurrEval.sh"
sed -i -e  "s?threads=2?threads=$geantthreads?g" "$homedir/SurrEval.sh"
totalthreads=$(awk -v dd=$nnodes -v ee=$numproc 'BEGIN {printf "%d\n",dd*ee}');	
sed -i -e  "s?SpartOpt.dat?$1.dat?g" $dakotafile

# Set the maximum number of interations (for testing) and summary file name
#sed -i -e  "s?max_iterations = 5?max_iterations = $max_iterations?g" $dakotafile
sed -i -e  "s?Summary.dat?$1Summary.dat?g" $dakotafile

# Create Dakota parameters
if [[ -f "$1Dakota.template" ]]
then
	rm "$1Dakota.template"
fi
#**********************************************************************
# MOGA parameters
#**********************************************************************
# Uncomment lines (and comment corresponding lines) to run Original design
echo " " >> $dakotafile
echo "variables" >> $dakotafile
echo "  id_variables = 'V1'" >> $dakotafile
# Choose fin and wall thicknesses
echo "    continuous_design = 2" >> $dakotafile
#echo "        initial_point 5.0	 7.5" >> $dakotafile
#echo "        upper_bounds 5.0	7.5" >> $dakotafile
#echo "        lower_bounds 5.0	7.5" >> $dakotafile
echo "        initial_point $min_mask	$min_mask" >> $dakotafile
echo "        upper_bounds $max_mask	$max_mask" >> $dakotafile
echo "        lower_bounds $min_mask	$min_mask" >> $dakotafile
echo "        descriptors 'finthick'       'wallthick'" >> $dakotafile

# Choose angular discretizations
echo "    discrete_design_set real = 2" >> $dakotafile
#echo "        num_set_values 1 1" >> $dakotafile
#echo "        set_values 10  10" >> $dakotafile
echo "        num_set_values 8 8" >> $dakotafile
echo "        set_values 2 3 4 5 6 8 9 10 2 3 4 5 6 8 9 10" >> $dakotafile
echo "        descriptors 'deltatheta'       'deltaphi'" >> $dakotafile

# Choose fin and wall widths
echo "    discrete_design_range = 2" >> $dakotafile
#echo "        upper_bounds 2 1" >> $dakotafile
#echo "        lower_bounds 2 1" >> $dakotafile
echo "        upper_bounds 12 12" >> $dakotafile
echo "        lower_bounds 1 1" >> $dakotafile
echo "        descriptors 'finwidth'        'wallwidth'" >> $dakotafile

#**********************************************************************
# MAS parameters
#**********************************************************************
echo " " >> $dakotafile
echo "variables" >> $dakotafile
echo "  id_variables = 'V2'" >> $dakotafile
echo "    continuous_design = 2" >> $dakotafile
# Choose fin and wall thicknesses
echo "        initial_point $min_mask	$min_mask" >> $dakotafile
echo "        upper_bounds $max_mask	$max_mask" >> $dakotafile
echo "        lower_bounds $min_mask	$min_mask" >> $dakotafile
echo "        descriptors 'finthick'       'wallthick'" >> $dakotafile

# Choose angular discretizations
echo "    discrete_design_set real = 2" >> $dakotafile
echo "        num_set_values 8 8" >> $dakotafile
echo "        set_values 2 3 4 5 6 8 9 10 2 3 4 5 6 8 9 10" >> $dakotafile
echo "        descriptors 'deltatheta'       'deltaphi'" >> $dakotafile

# Choose fin and wall widths
echo "    discrete_design_range = 2" >> $dakotafile
echo "        upper_bounds 12 12" >> $dakotafile
echo "        lower_bounds 1 1" >> $dakotafile
echo "        descriptors 'finwidth'        'wallwidth'" >> $dakotafile

# Set input/output
sed -i -e  "s?objective_functions = 2?objective_functions = $n_obj?g" $dakotafile
sed -i -e  "s?evaluation_concurrency = 1?evaluation_concurrency = $concurrancy?g" $dakotafile
sed -i -e  "s?MOGAworkdir?$homedir/$1/MOGAbld?g" $dakotafile
sed -i -e  "s?MASworkdir?$homedir/$1/MASbld?g" $dakotafile
sed -i -e  "s?MOGAparameters.in?$1.MOGApara?g" $dakotafile
sed -i -e  "s?MASparameters.in?$1.MASpara?g" $dakotafile
sed -i -e  "s?MOGAresults.out?$1.MOGAres?g" $dakotafile
sed -i -e  "s?MASresults.out?$1.MASres?g" $dakotafile

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run Dakota script
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run on computer
dakota $dakotafile > $1.log
#**********************************************************************
# Run on cluster
# Create PBS job submission script
#cp submitGeant.pbs DakotarunTEMP.pbs
# Submit job
#./batchDakota.sh "$dakotafile" $nnodes $numproc $totalthreads $1
