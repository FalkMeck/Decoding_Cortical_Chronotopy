#!/bin/bash
# author: Dominik Grotegerd
# institution: Translational Psychiarty / University of Muenster
# Version 1.1 / 2018

subDir=$(readlink -e $1) # uses first character string pasted after running shell from command as input for subDir
# necessery to loop through code via Matalb unix()-command
subBase=$(basename $subDir) # basename takes the last part of subDir (here name of participant) for subBase

# overwrite copy flag 
cpFlag="-i" # unclear what this is used for since it just defines a variable and 
cpFlag=""

# FSL configuration
#/usr/local/fsl/etc/fslconf/fsl.sh # maybe needs to be changed, but is just comment not code ???

# Naming Conventions # Defining Names and Folders for variables that get created in shell script
MagImName="$subBase"_MagnitudeImage 
phaImName="$subBase"_PhaseImage
revImName="$subBase"_ReversedImage
mergedImName="$subBase"_dwi_merged.nii.gz
mergedImNameRPE="$subBase"_dwi_topuped_merged.nii.gz
B0_ref="$subBase"_B0reference.nii.gz
outDir=$subDir/DTI/DWI_prepared 

bvName=$(echo $mergedImName | sed 's/.nii.gz//g') # taking the mergedImName(RPE) and cuting of the ".nii.gz"-Ending
bvNameRPE=$(echo $mergedImNameRPE | sed 's/.nii.gz//g') # to create new Names for the bval and bvec variables


# simple check if Subject Directory/Path exists in the specified form
if [ ! -d "$subDir"/DTI/ ]; then  
	echo $subDir DTI directory does not exist  
	exit 1
fi

mkdir -p $outDir # create the output directory $subDir/DTI/DWI_prepared 
cd $outDir # set current directory to that path


# mag image / phase image # define the magnitude and pahse image from the gradient fields
magIm="$subDir"/DTI/GreField_1/$(ls "$subDir"/DTI/GreField_1/ | head -1) 
phaIm="$subDir"/DTI/GreField_2/*.nii
revIm="$subDir"/DTI/ReversePhaseEncoding/*.nii

# get FieldMap Images
cp $cpFlag $phaIm $phaImName.nii # copy gradient niftis to new file (name) specified in Naming Conventions
cp $cpFlag $magIm $MagImName.nii
cp $cpFlag $revIm $revImName.nii

# merge simple dwi
$FSLDIR/bin/fslmerge -t "$mergedImName" "$subDir"/DTI/DTI_1/*.nii "$subDir"/DTI/DTI_2/*.nii "$subDir"/DTI/DTIB0/*nii
#/usr/local/fsl/bin/fslmerge -t "$mergedImName" "$subDir"/DTI/DTI_1/*.nii "$subDir"/DTI/DTI_2/*.nii "$subDir"/DTI/DTIB0/*nii
# merge the multiple DWI nifitis to one file

# merge bval / bvecs
paste   "$subDir"/DTI/DTI_1/*bvec "$subDir"/DTI/DTI_2/*bvec  <(echo 0 0 0; echo 0 0 0; echo 0 0 0)  | column -s $'\t' -t > "$bvName".bvec
paste   "$subDir"/DTI/DTI_1/*bval "$subDir"/DTI/DTI_2/*bval  <(echo 0 0 0) | column -s $'\t' -t > "$bvName".bval
# paste the two bvec and beval files into one with 0 0 0 in the third column
# so 3 columns 3(bvec) or 1(bval) lines
# 1st col = all 3/1 lines from DTI_1 files
# 2nd col = all 3/1 lines from DTI_2 files
# 3rd col = 0 0 0 the three B0 Images

# if there's Reverse Phase Encoding
if [ -d "$subDir"/DTI/ReversePhaseEncoding ]; then
	$FSLDIR/bin/fslmerge -t "$mergedImNameRPE"  "$mergedImName" "$subDir"/DTI/ReversePhaseEncoding/*.nii
	# /usr/local/fsl/bin/fslmerge -t "$mergedImNameRPE"  "$mergedImName" "$subDir"/DTI/ReversePhaseEncoding/*.nii
	paste "$bvName".bvec <(echo -0 -0 -0 -0; echo -0 -0 -0 -0; echo -0 -0 -0 -0)  | column -s $'\t' -t > "$bvNameRPE".bvec
	paste "$bvName".bval <(echo -0 -0 -0 -0 ) | column -s $'\t' -t > "$bvNameRPE".bval
fi
# If a ReversePhaseEncoding-Image exists then merge the merged files from DTI_1 and DTI_2 and DTIB0 with this image # and create new files in the same way which now have the filename "dwi_topuped_merged"
# for pasting add column with -0 -0 -0 -0 to represent 4 reverese(-) Images
