#! /bin/bash                                                                                                                                                                                               

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

setupATLAS
asetup 23.6.60,AthGeneration
source setupRivet

export RIVET_PLOT_PATH=$PWD
