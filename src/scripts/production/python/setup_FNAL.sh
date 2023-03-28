#!/bin/bash

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup ifdhc v2_6_6
export IFDH_CP_MAXRETRIES=0

setup root v6_12_06a -q e17:debug 
setup lhapdf v5_9_1k -q e17:debug 
setup log4cpp v1_1_3a -q e17:debug 
setup pdfsets v5_9_1b 
setup gdb v8_1
setup git v2_15_1
setup boost v1_66_0a -q e17:prof
setup cmake v3_14_3


