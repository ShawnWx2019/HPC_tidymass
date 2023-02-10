#!/bin/bash

##################################################
#       Prj: Tidymass 
#       Assignment: Convert raw data to mzXML
#       Author: Shawn
#       Date: Jul 1, 2022
#       Location: HENU
###################################################


echo -e "\033[31mStart convert .raw file to .mzXML file, This will take a while. \033[0m"


docker run  --rm -e WINEDEBUG=-all -v $1:/data -v $2:/outpath chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /data/*.raw --ignoreUnknownInstrumentError --mzXML --64 --zlib --filter "peakPicking cwt snr=0.1 peakSpace=0.1 msLevel=1-" --filter "scanTime [60,1080]" --filter "msLevel 1-2" --filter "zeroSamples removeExtra 1-" -o /outpath/

echo -e "\033[31mFinish. \033[0m"

