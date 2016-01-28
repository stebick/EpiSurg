#!/bin/sh

usage='\nUSAGE:\n  ct2mriBbreg.sh freesurferSubjectName ctNiiFile\n\nEXAMPLE:\n ct2mriBbreg.sh TWH014 /Users/dgroppe/Desktop/TWH_14_DICOMS/ct.nii.gz\n'

if [[ "$#" = 0 ]]; then
 echo -e $usage
 exit 2
fi


# ct2mri.sh
# 
# Registers a CT to an MRI using bbregister
#
# Created by David Groppe on 2/11/15.
# Questions? Email: david.m.groppe@gmail.com
# Copyright 2015 __MyCompanyName__. All rights reserved.

#echo 'Registering postopCT.nii.gz to T1.nii.gz'
#echo Subject $1

sub=$1
fsDir=$SUBJECTS_DIR$sub
if [ ! -d $fsDir ]; then
  echo
  echo "... ${fsDir} is not a directory."
  echo "...you have the wrong FreeSurfer codename for the subject or you have not yet run recon-all on this subject's MRI " 
  echo
  exit 2
fi

if [ ! -f  $2 ]; then
 echo
 echo "...File ${2} not found. Exit."
 echo
 exit 2
fi

elecReconPath=$SUBJECTS_DIR$sub/elec_recon
echo 'Creating directory ' $elecReconPath
mkdir $elecReconPath

echo 'Creating T1.nii.gz in elec_recon folder for coregistration.'
mriPath=$SUBJECTS_DIR$sub/mri
mri_convert $mriPath/T1.mgz $elecReconPath/T1.nii.gz

echo 'Creating brainmask.nii.gz in elec_recon folder for use with BioImageSuite later.'
mri_convert $mriPath/brainmask.mgz $elecReconPath/brainmask.nii.gz

echo 'Copying CT nii.gz file to elec_recon folder.'
cp $2 $elecReconPath/.

#bbregister --s $sub --mov $2 --reg $elecReconPath/ct2mri.dat --fslmat $elecReconPath/ct2mri.mat --init-fsl --bold
#flirt -in $2 -ref $elecReconPath/T1.nii.gz -out $elecReconPath/ctINt1.nii.gz -interp trilinear -init $elecReconPath/ct2mri.mat -applyxfm
#slices $elecReconPath/ctINt1.nii.gz $elecReconPath/T1.nii.gz
#slices $elecReconPath/T1.nii.gz  $elecReconPath/ctINt1.nii.gz  
echo 'Run this for interactive GUI'
echo 'fslview ' $elecReconPath/T1.nii.gz ' ' $elecReconPath/ctINt1.nii.gz  
