#!/bin/bash

############################################################
#       Prj: Tidymass pipeline
#       Assignment: initialization
#       Author: Shawn Wang
#       Date: Feb 7 2023
############################################################

set -e ## 报错打断，防止一直错下去

echo -e "\033[32mHPC_tidyMass initialization start....\033[0m"

echo -e "\033[32m-------------------------------\033[0m"

##> conf file

mkdir ~/.HPC_tidymass && cd ~/.HPC_tidymass

cp -r src ~/.HPC_tidymass/src

cp -r MS_db ~/.HPC_tidymass/MS_db

##> alias

touch ~/.bash_alias && echo "source ~/.bash_alias" >> ~/.bashrc

echo 'alias hpc-runTidymass="bash ~/.HPC_tidymass/src/shell/runTidymass.sh"' >> ~/.bash_alias
echo 'alias hpc-msConvert="bash ~/.HPC_tidymass/src/shell/01.msconvert.sh"' >> ~/.bash_alias
echo 'alias hpc-msConvert2="bash ~/.HPC_tidymass/src/shell/03.ms2convert.sh"' >> ~/.bash_alias
echo 'alias hpc-peakPicking="Rscript ~/.HPC_tidymass/src/R/02.PeakPicking.R"' >> ~/.bash_alias
echo 'alias hpc-dataCleaning="Rscript ~/.HPC_tidymass/src/R/04.DataCleaning.R"' >> ~/.bash_alias
echo 'alias hpc-annotation="~/Rscript .HPC_tidymass/src/R/05.Metabolomics_annotation.R"' >> ~/.bash_alias
echo 'alias hpc-annoFiltering="Rscript ~/.HPC_tidymass/src/R/06.Annotation_filtering.R"' >> ~/.bash_alias

source ~/.bash_alias

echo -e "\033[32mFinish!\033[0m"