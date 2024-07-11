#!/bin/bash -p
export FREESURFER_HOME=.../freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export FS_LICENSE=$FREESURFER_HOME/license.txt


declare -a subjectNames=("HFG_121")

declare -a templates=("lausanne250")


for subj in "${subjectNames[@]}" 
do
	SUBJECTS_DIR=.../Analysis
	cd $SUBJECTS_DIR
	echo $subj

	for parc in "${templates[@]}" 
	do
		outDir=$SUBJECTS_DIR/$subj/FreeSurfer/$parc
		mkdir $outDir
		echo $parc
	
		Hem=lh
		mri_annotation2label --subject $subj/FreeSurfer \
		--hemi $Hem \
		--labelbase $outDir/$parc-$Hem \
		--annotation $parc

		Hem=rh
		mri_annotation2label --subject $subj/FreeSurfer \
		--hemi $Hem \
		--labelbase $outDir/$parc-$Hem \
		--annotation $parc
	
	done
done
