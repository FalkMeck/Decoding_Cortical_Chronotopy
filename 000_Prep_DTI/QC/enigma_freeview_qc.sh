#!/bin/bash -p

STUDY_DIR=.../Freesurfer_data
subject=$1
SUB_DIR=$STUDY_DIR/$subject

source $FREESURFER_HOME/SetUpFreeSurfer.sh
export FS_LICENSE=$STUDY_DIR/license.txt 

#export LIBGL_ALWAYS_INDIRECT=1

# opens two freeview windows back to back
# Internal
freeview -v \
	$SUB_DIR/mri/orig.mgz \
	$SUB_DIR/mri/aparc+aseg.mgz:colormap=lut:opacity=0.4 \
	
# External
freeview -f \
	$SUB_DIR/surf/lh.pial:annot=aparc.annot:name=pial_aparc:visible=0 \
	$SUB_DIR/surf/rh.pial:annot=aparc.annot:name=pial_aparc:visible=0

# same images as  in the internal/extern QC websites
