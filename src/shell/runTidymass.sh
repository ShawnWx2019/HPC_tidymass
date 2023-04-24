#!/bin/bash

############################################################
#       Prj: Tidymass pipeline
#       Assignment: Convert ms format and peak detection
#       Author: Shawn Wang
#       Date: Feb 7 2023
############################################################

## getopts
set -e ## 报错打断，防止一直错下去

start_time=$(date +%s)

## 帮助内容
func(){
    echo -e "\033[32mTidymass pipeline part1. Date transform ,peak picking and annotation\033[0m"
    echo -e "\033[32mUsage:\033[0m"
    echo -e "\033[32m-------------------------------\033[0m"
    echo -e "\033[35mrunTidymass \033[32m[-i input] \033[33m[-t type] \033[31m[-c column]\033[0m"
    echo -e "\033[32m-------------------------------\033[0m"
    echo -e "\033[32mAuthor\033[0m Shawn Wang (shawnwang2016@126.com)"
    echo -e "\033[32m-------------------------------\n\033[0m"
    echo -e "\033[32mDate\033[0m Tus Feb 09, 2023"
    echo -e "\033[32m-------------------------------\n\033[0m"
    echo -e "\033[32mVersion\033[0m V.0.0.0.99 beta"
    echo -e "\033[32m-------------------------------\n\033[0m"
    echo -e "\033[32mDescription\033[0m"
    echo -e "\033[32m-------------------------------\n\033[0m"
    echo -e "\033[32m[-i]:input\033[0m,   The file path of .raw data\033[0m"
    echo -e "\033[33m[-t]:type\033[0m,  The type of ion model, 1: NEG+POS in one file, 2: NEG and POS in differnet files\033[0m"
    echo -e "\033[31m[-c]:column\033[0m,  rp or hilic"
    exit -1
}
## 是否备份原始文件，如果不写-b默认不备份

## 设置各个参数
while getopts 'h:i:t:c:' OPT;
do
    case $OPT in
        i) input=`echo "$OPTARG"`;;
        t) type=`echo "$OPTARG"`;;
        c) column=`echo "$OPTARG"`;;
        h) func;;
        ?) func;;
    esac
done
## 提示你的参数设置
echo -e "\033[32m====================================================\033[1m"
echo -e "\033[34mYour setting:\033[0m"
echo -e "\033[32mThe input file located in:\033[0m ${input}" 
echo -e "\033[32mThe input ion model is :\033[0m  ${type}"
echo -e "\033[32mThe column type is:\033[0m  ${column}"
echo -e "\033[32m====================================================\033[0m"

## step1. Project init

echo -e "\033[32m====================================================\nStep1 Project init....\n==================================================== \033[0m"

mkdir -p working_dir && cd working_dir && mkdir -p 01.data 02.progress 03.result

main_dir=$PWD

filecount=`ls 01.data | wc -w`

if [ $filecount == "0" ]
then
    echo -e "Copy raw file to working dir, please wait... "
    cp -r $input 01.data/rawdata && cd 01.data/rawdata
else
    cd 01.data/rawdata
fi  

## data organization

if [ -d "QC" ] || [ -d "Subject" ] || [ -d "POS/QC" ] || [ -d "NEG/QC" ] || [ -d "POS/Subject" ] || [ -d "NEG/Subject" ]
then
    echo -e "\033[34mData orgniazation have already fininshed! start next step.\033[0m"
else
    if [ $type == "1" ] 
    then
        ls -ltr $input > injection_order_pos.txt && mv injection_order_pos.txt ../ && cp ../injection_order_pos.txt ../injection_order_neg.txt
        mkdir -p QC Subject
        mv QC*.raw QC && mv *.raw Subject
        pos_dir_raw=$(realpath ./)
        neg_dir_raw=$(realpath ./)
        echo -e "\033[34mProcess: Data organization finishi waiting for data trans-format!\033[0m"
        echo -e "\033[34mThe rawdata were placed at\033[0m $pos_dir_raw"
        tree -d
    elif [ $type == "2" ]
    then
        ## 进行负谱数据整合 
        cd NEG && ls -ltr $input/NEG > injection_order_neg.txt && mv injection_order_neg.txt ../
        mkdir -p QC Subject 
        mv QC*.raw QC && mv *.raw Subject 
        neg_dir_raw=$(realpath ./)
        ## 进行正谱数据整合 
        cd ../POS && ls -ltr $input/POS > injection_order_pos.txt && mv injection_order_pos.txt ../
        mkdir -p QC Subject
        mv QC*.raw QC && mv *.raw Subject
        pos_dir_raw=$(realpath ./) && cd ../
        tree -d
        echo -e "\033[34mProcess: Data organization finishi waiting for data trans-format!\033[1m"
        echo -e "\033[34mThe rawdata (POS) were placed at\033[0m $pos_dir_raw"
        echo -e "\033[34mThe rawdata (NEG) were placed at\033[0m $neg_dir_raw"
    else
        echo -e "\033[41;37mERROR:Wrong type, please check parameter -t, only accept 1 or 2!\033[0m"
        exit
    fi
fi


## step2. convert .raw to mzXML and .mgf format

echo -e "\033[32m====================================================\nStep2 Convert .raw file to .mzXML and .mgf\n==================================================== \033[0m"

cd ../../02.progress && mkdir -p transform peak_picking && cd transform

filecount2=`ls ./ | wc -w`
out_put_dir=$(realpath ./)
if [ $filecount2 == "0" ]
then
    mkdir -p ${out_put_dir}/MS1/POS/QC ${out_put_dir}/MS1/POS/Subject ${out_put_dir}/MS1/NEG/QC ${out_put_dir}/MS1/NEG/Subject
    mkdir -p ${out_put_dir}/MS2/POS/QC ${out_put_dir}/MS2/POS/Subject ${out_put_dir}/MS2/NEG/QC ${out_put_dir}/MS2/NEG/Subject
    echo -e "\033[32m====================================================\nParalle run.\n==================================================== \033[0m"
    commands=(
        "bash ~/01.src/02.script/02.Tidymass/01.msconvert.sh ${pos_dir_raw}/QC  ${out_put_dir}/MS1/POS/QC"
        "bash ~/01.src/02.script/02.Tidymass/01.msconvert.sh ${pos_dir_raw}/Subject  ${out_put_dir}/MS1/POS/Subject"
        "bash ~/01.src/02.script/02.Tidymass/01.msconvert.sh ${neg_dir_raw}/QC  ${out_put_dir}/MS1/NEG/QC"
        "bash ~/01.src/02.script/02.Tidymass/01.msconvert.sh ${neg_dir_raw}/Subject  ${out_put_dir}/MS1/NEG/Subject"
        "bash ~/01.src/02.script/02.Tidymass/03.ms2convert.sh ${pos_dir_raw}/QC ${out_put_dir}/MS2/POS/QC"
        "bash ~/01.src/02.script/02.Tidymass/03.ms2convert.sh ${pos_dir_raw}/Subject ${out_put_dir}/MS2/POS/Subject"
        "bash ~/01.src/02.script/02.Tidymass/03.ms2convert.sh ${neg_dir_raw}/QC ${out_put_dir}/MS2/NEG/QC"
        "bash ~/01.src/02.script/02.Tidymass/03.ms2convert.sh ${neg_dir_raw}/Subject ${out_put_dir}/MS2/NEG/Subject"
    )
    parallel --jobs 0 ::: "${commands[@]}"
    echo -e "\033[32mTransform finish!\033[0m"
    
else
    echo -e "\033[42;37malreay finished! skip this step.\033[0m"
fi  


echo -e "\033[32m====================================================\nTransformat finish! run.\n==================================================== \033[0m"

echo -e "\033[32m====================================================\nStep3 peak picking\n==================================================== \033[0m"

cd ../peak_picking

## Step3 peak picking...

filecount3=`ls ./ | wc -w`

if [ $filecount3 == "0" ]
then
    echo -e "\033[32m====================================================\nPeak picking run.\n==================================================== \033[0m"
    /usr/local/bin/Rscript ~/01.src/02.script/02.Tidymass/02.PeakPicking.R \
                -p ../transform/MS1/POS/ \
                -n ../transform/MS1/NEG/ \
                -t 100 \
                -g "QC"
    mkdir -p NEG POS && mv ${out_put_dir}/MS1/NEG/Result NEG/Result && mv  ${out_put_dir}/MS1/POS/Result POS/Result
    cp NEG/Result/object object.neg && cp POS/Result/object object.pos 
    
else
    echo -e "\033[42;37malreay finished! skip this step.\033[0m"
fi  

echo -e "\033[32m====================================================\nPeak picking done!\n==================================================== \033[0m"

echo -e "\033[32m====================================================\nStep4 Data cleaning\n==================================================== \033[0m"

## Step4 Data cleaning

cd $main_dir/02.progress

if [ -d Data_cleaning ]
then
    echo -e "\033[34mData cleaning have already fininshed! start next step.\033[0m"
else
    echo -e "\033[34mStart data cleaing progress\033[0m"
    mkdir -p Data_cleaning && cd Data_cleaning 
    ln -s ../../01.data/rawdata/injection_order_* ./
    ln -s ../peak_picking/object* ./
    /usr/local/bin/Rscript ~/01.src/02.script/02.Tidymass/04.DataCleaning.R 
fi

echo -e "\033[32m====================================================\nData cleaning done!\n==================================================== \033[0m"
echo -e "\033[32m====================================================\nStep5. Feature annotation!\n==================================================== \033[0m"

cd $main_dir/02.progress

if [ -d Annotation ]
then
    echo -e "\033[34mFeature annotation have already fininshed! start next step.\033[0m"
else
    echo -e "\033[34mStart feature annotation progress\033[0m"
    mkdir -p Annotation && cd Annotation 
    ln -s ../transform/MS2/NEG/ ./
    ln -s ../transform/MS2/POS/ ./
    ln -s ../Data_cleaning/*.rds ./
    /usr/local/bin/Rscript ~/01.src/02.script/02.Tidymass/05.Metabolomics_annotation.R $column
fi

cd $main_dir/02.progress/Annotation

##> annotation filter

if [ -d Clean_annotation/Only_MS2 ] 
then 
    echo -e "033[34mAnnotation filter have already fininshed! start next step.\033[0m"
else
    echo -e "\033[34mStart feature annotation filtering progress\033[0m"
    /usr/local/bin/Rscript ~/01.src/02.script/02.Tidymass/06.Annotation_filtering.R
fi

echo -e "\033[32m====================================================\nFeature annotation done!\n==================================================== \033[0m"
echo -e "\033[32m====================================================\nStep6. Down Stream analysis!\n==================================================== \033[0m"


##> check running time
end_time=$(date +%s)
run_time=$((end_time - start_time))


hours=$((run_time / 3600))
minutes=$(((run_time / 60) % 60))
seconds=$((run_time % 60))

printf "\033[34mThe script took %d hours, %d minutes, and %d seconds to run.\033[0m\n" $hours $minutes $seconds