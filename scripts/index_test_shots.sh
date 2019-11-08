###############################################################################
# Index test shots from 14-BM. This presupposes you have taken 5 images per   #
# scan                                                                        #
###############################################################################


datadir=`readlink -f $1`
workdir=/tmp/hekstra/$datadir
reference_geo=/home/userbmc/processing/reference_geometry.expt
space_group=19

mkdir -p $workdir
cd $workdir

dials.import invert_rotation_axis=True \
	reference_geometry=$reference_geo \
	image_range=1,5 \
	$datadir/*.cbf

dials.find_spots imported.* nproc=10
dials.index imported.expt strong.refl space_group=$space_group

dials.refine scan_varying=False indexed.* 
dials.integrate refined.* nproc=10

