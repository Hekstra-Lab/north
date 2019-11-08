#!/bin/bash

#########################################################################
# Run dials beam search to make a reference geometry. The input should  #
# ideally be a full, high quality diffraction data set from a single    #
# crystal like lysozyme or dhfr.                                        #
#########################################################################


datadir=/mnt/ChemData2/Data/BIOCARS/hekstra_201911/actual_dhfr_n23pp
space_group=19


mkdir .tmp
cd .tmp
dials.import invert_rotation_axis=True $datadir/*.cbf
dials.find_spots imported.* nproc=50
dials.search_beam_position strong.refl imported.expt n_macro_cycles=3 max_reflections=70000 nproc=50
dials.index space_group=$space_group optimised.expt strong.refl
dials.refine scan_varying=True indexed.*
mv refined.expt ../reference_geometry.expt
cd ..
rm -r .tmp
