#!/bin/bash

export SUBJECTS_DIR=.../Freesurfer_data/2_Anatomy_Reoriented/Day2/

fsqcDIR=.../Freesurfer_data/2_Anatomy_Reoriented/Day2/fsqc/

mkdir $fsqcDIR

cd $fsqcDIR

source .../fsqc.sh
