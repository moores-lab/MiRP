#!/usr/bin/env csh 

#Joe Atherton 30/01/19#

#TAKES RELION .MRCS STACKS AND .STAR FILES (IN EXTRACT FOLDER) WITH SEVERAL MTS AND THEIR SEGMENTS, DECONSTRUCTS INTO INVDIVDUAL MTS THEN MAKES SEGMENT AVERAGES FOR EACH MT, THEN RESTACKS INTO ONE SEGMENT AVERAGE STACK PER MT

#Run as: source preprocess_segment_averages.csh 108, where 108 is the boxsize


##############################################################
#Will need to source imod, bsoft and relion/v3.0/beta4. Modify the below accordingly for your institute's computing set up.

module use -a /s/emib/s/modules
module load eman2
module load imod
module load bsoft
module load relion/v3.0/beta4

##############################################################

##############################################################

#calculate background radius box for normalisation

set box_size=$1

set background_box_radius = `echo $box_size | awk '{print 0.75*$1/2}'`

##############################################################

#set up micrograph loop

foreach micrograph_stack (`ls -1 ./*.mrcs`)
echo 'working on' $micrograph_stack

##############################################################

#work out number of segments in stack from .star file.

set total_number_of_segments = `grep -o '.mrcs' $micrograph_stack:r_extract.star | wc -l`

#works out number of segments in each microtubule from micrograph .star file. Only for MAX 30 MTs so far. Dangerous using grep, as depends on column separation.

set number_of_segments_MT1 = `grep -o '            1.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT2 = `grep -o '            2.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT3 = `grep -o '            3.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT4 = `grep -o '            4.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT5 = `grep -o '            5.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT6 = `grep -o '            6.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT7 = `grep -o '            7.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT8 = `grep -o '            8.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT9 = `grep -o '            9.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT10 = `grep -o '           10.*.mrcs' $micrograph_stack:r_extract.star | wc -l`

set number_of_segments_MT11 = `grep -o '           11.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT12 = `grep -o '           12.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT13 = `grep -o '           13.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT14 = `grep -o '           14.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT15 = `grep -o '           15.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT16 = `grep -o '           16.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT17 = `grep -o '           17.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT18 = `grep -o '           18.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT19 = `grep -o '           19.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT20 = `grep -o '           20.*.mrcs' $micrograph_stack:r_extract.star | wc -l`

set number_of_segments_MT21 = `grep -o '           21.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT22 = `grep -o '           22.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT23 = `grep -o '           23.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT24 = `grep -o '           24.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT25 = `grep -o '           25.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT26 = `grep -o '           26.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT27 = `grep -o '           27.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT28 = `grep -o '           28.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT29 = `grep -o '           29.*.mrcs' $micrograph_stack:r_extract.star | wc -l`
set number_of_segments_MT30 = `grep -o '           30.*.mrcs' $micrograph_stack:r_extract.star | wc -l`


##############################################################

#determine segment ranges for each MT. DEPENDS ON MTS BEING NUMBERED SEQUENTIALLY IN .STAR FILE (DOES NOT WORK IF HELICAL ID SKIPS FROM 1 TO 3). Only for MAX 30 MTs so far.

set MT1_start_section = 1
set MT1_end_section = $number_of_segments_MT1

set MT2_start_section = `echo $MT1_end_section | awk '{print $1+1}'`
set MT2_end_section = `echo $MT1_end_section $number_of_segments_MT2 | awk '{print $1+$2}'`

set MT3_start_section = `echo $MT2_end_section | awk '{print $1+1}'`
set MT3_end_section = `echo $MT2_end_section $number_of_segments_MT3 | awk '{print $1+$2}'`

set MT4_start_section = `echo $MT3_end_section | awk '{print $1+1}'`
set MT4_end_section = `echo $MT3_end_section $number_of_segments_MT4 | awk '{print $1+$2}'`

set MT5_start_section = `echo $MT4_end_section | awk '{print $1+1}'`
set MT5_end_section = `echo $MT4_end_section $number_of_segments_MT5 | awk '{print $1+$2}'`

set MT6_start_section = `echo $MT5_end_section | awk '{print $1+1}'`
set MT6_end_section = `echo $MT5_end_section $number_of_segments_MT6 | awk '{print $1+$2}'`

set MT7_start_section = `echo $MT6_end_section | awk '{print $1+1}'`
set MT7_end_section = `echo $MT6_end_section $number_of_segments_MT7 | awk '{print $1+$2}'`

set MT8_start_section = `echo $MT7_end_section | awk '{print $1+1}'`
set MT8_end_section = `echo $MT7_end_section $number_of_segments_MT8 | awk '{print $1+$2}'`

set MT9_start_section = `echo $MT8_end_section | awk '{print $1+1}'`
set MT9_end_section = `echo $MT8_end_section $number_of_segments_MT9 | awk '{print $1+$2}'`

set MT10_start_section = `echo $MT9_end_section | awk '{print $1+1}'`
set MT10_end_section = `echo $MT9_end_section $number_of_segments_MT10 | awk '{print $1+$2}'`

#####

set MT11_start_section = `echo $MT10_end_section | awk '{print $1+1}'`
set MT11_end_section = `echo $MT10_end_section $number_of_segments_MT11 | awk '{print $1+$2}'`

set MT12_start_section = `echo $MT11_end_section | awk '{print $1+1}'`
set MT12_end_section = `echo $MT11_end_section $number_of_segments_MT12 | awk '{print $1+$2}'`

set MT13_start_section = `echo $MT12_end_section | awk '{print $1+1}'`
set MT13_end_section = `echo $MT12_end_section $number_of_segments_MT13 | awk '{print $1+$2}'`

set MT14_start_section = `echo $MT13_end_section | awk '{print $1+1}'`
set MT14_end_section = `echo $MT13_end_section $number_of_segments_MT14 | awk '{print $1+$2}'`

set MT15_start_section = `echo $MT14_end_section | awk '{print $1+1}'`
set MT15_end_section = `echo $MT14_end_section $number_of_segments_MT15 | awk '{print $1+$2}'`

set MT16_start_section = `echo $MT15_end_section | awk '{print $1+1}'`
set MT16_end_section = `echo $MT15_end_section $number_of_segments_MT16 | awk '{print $1+$2}'`

set MT17_start_section = `echo $MT16_end_section | awk '{print $1+1}'`
set MT17_end_section = `echo $MT16_end_section $number_of_segments_MT17 | awk '{print $1+$2}'`

set MT18_start_section = `echo $MT17_end_section | awk '{print $1+1}'`
set MT18_end_section = `echo $MT17_end_section $number_of_segments_MT18 | awk '{print $1+$2}'`

set MT19_start_section = `echo $MT18_end_section | awk '{print $1+1}'`
set MT19_end_section = `echo $MT18_end_section $number_of_segments_MT19 | awk '{print $1+$2}'`

set MT20_start_section = `echo $MT19_end_section | awk '{print $1+1}'`
set MT20_end_section = `echo $MT19_end_section $number_of_segments_MT20 | awk '{print $1+$2}'`

#####

set MT21_start_section = `echo $MT20_end_section | awk '{print $1+1}'`
set MT21_end_section = `echo $MT20_end_section $number_of_segments_MT21 | awk '{print $1+$2}'`

set MT22_start_section = `echo $MT21_end_section | awk '{print $1+1}'`
set MT22_end_section = `echo $MT21_end_section $number_of_segments_MT22 | awk '{print $1+$2}'`

set MT23_start_section = `echo $MT22_end_section | awk '{print $1+1}'`
set MT23_end_section = `echo $MT22_end_section $number_of_segments_MT23 | awk '{print $1+$2}'`

set MT24_start_section = `echo $MT23_end_section | awk '{print $1+1}'`
set MT24_end_section = `echo $MT23_end_section $number_of_segments_MT24 | awk '{print $1+$2}'`

set MT25_start_section = `echo $MT24_end_section | awk '{print $1+1}'`
set MT25_end_section = `echo $MT24_end_section $number_of_segments_MT25 | awk '{print $1+$2}'`

set MT26_start_section = `echo $MT25_end_section | awk '{print $1+1}'`
set MT26_end_section = `echo $MT25_end_section $number_of_segments_MT26 | awk '{print $1+$2}'`

set MT27_start_section = `echo $MT26_end_section | awk '{print $1+1}'`
set MT27_end_section = `echo $MT26_end_section $number_of_segments_MT27 | awk '{print $1+$2}'`

set MT28_start_section = `echo $MT27_end_section | awk '{print $1+1}'`
set MT28_end_section = `echo $MT27_end_section $number_of_segments_MT28 | awk '{print $1+$2}'`

set MT29_start_section = `echo $MT28_end_section | awk '{print $1+1}'`
set MT29_end_section = `echo $MT28_end_section $number_of_segments_MT29 | awk '{print $1+$2}'`

set MT30_start_section = `echo $MT29_end_section | awk '{print $1+1}'`
set MT30_end_section = `echo $MT29_end_section $number_of_segments_MT30 | awk '{print $1+$2}'`



##############################################################


#Splits original .mrcs segments stack for each micrograph in to individual microtubules.IN PROGRESS. Only for MAX 30 MTs so far.

if ($number_of_segments_MT1 >1) then
newstack -NumberedFromOne -secs $MT1_start_section-$MT1_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT01.mrcs
echo $number_of_segments_MT1 > $micrograph_stack:r_MT01.data
endif

if ($number_of_segments_MT2 >1) then
newstack -NumberedFromOne -secs $MT2_start_section-$MT2_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT02.mrcs
echo $number_of_segments_MT2 > $micrograph_stack:r_MT02.data
endif

if ($number_of_segments_MT3 >1) then
newstack -NumberedFromOne -secs $MT3_start_section-$MT3_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT03.mrcs
echo $number_of_segments_MT3 > $micrograph_stack:r_MT03.data
endif

if ($number_of_segments_MT4 >1) then
newstack -NumberedFromOne -secs $MT4_start_section-$MT4_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT04.mrcs
echo $number_of_segments_MT4 > $micrograph_stack:r_MT04.data
endif

if ($number_of_segments_MT5 >1) then
newstack -NumberedFromOne -secs $MT5_start_section-$MT5_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT05.mrcs
echo $number_of_segments_MT5 > $micrograph_stack:r_MT05.data
endif

if ($number_of_segments_MT6 >1) then
newstack -NumberedFromOne -secs $MT6_start_section-$MT6_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT06.mrcs
echo $number_of_segments_MT6 > $micrograph_stack:r_MT06.data
endif

if ($number_of_segments_MT7 >1) then
newstack -NumberedFromOne -secs $MT7_start_section-$MT7_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT07.mrcs
echo $number_of_segments_MT7 > $micrograph_stack:r_MT07.data
endif

if ($number_of_segments_MT8 >1) then
newstack -NumberedFromOne -secs $MT8_start_section-$MT8_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT08.mrcs
echo $number_of_segments_MT8 > $micrograph_stack:r_MT08.data
endif

if ($number_of_segments_MT9 >1) then
newstack -NumberedFromOne -secs $MT9_start_section-$MT9_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT09.mrcs
echo $number_of_segments_MT9 > $micrograph_stack:r_MT09.data
endif

if ($number_of_segments_MT10 >1) then
newstack -NumberedFromOne -secs $MT10_start_section-$MT10_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT10.mrcs
echo $number_of_segments_MT10 > $micrograph_stack:r_MT10.data
endif

###

if ($number_of_segments_MT11 >1) then
newstack -NumberedFromOne -secs $MT11_start_section-$MT11_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT11.mrcs
echo $number_of_segments_MT11 > $micrograph_stack:r_MT11.data
endif

if ($number_of_segments_MT12 >1) then
newstack -NumberedFromOne -secs $MT12_start_section-$MT12_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT12.mrcs
echo $number_of_segments_MT12 > $micrograph_stack:r_MT12.data
endif

if ($number_of_segments_MT13 >1) then
newstack -NumberedFromOne -secs $MT13_start_section-$MT13_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT13.mrcs
echo $number_of_segments_MT13 > $micrograph_stack:r_MT13.data
endif

if ($number_of_segments_MT14 >1) then
newstack -NumberedFromOne -secs $MT14_start_section-$MT14_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT14.mrcs
echo $number_of_segments_MT14 > $micrograph_stack:r_MT14.data
endif

if ($number_of_segments_MT15 >1) then
newstack -NumberedFromOne -secs $MT15_start_section-$MT15_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT15.mrcs
echo $number_of_segments_MT15 > $micrograph_stack:r_MT15.data
endif

if ($number_of_segments_MT16 >1) then
newstack -NumberedFromOne -secs $MT16_start_section-$MT16_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT16.mrcs
echo $number_of_segments_MT16 > $micrograph_stack:r_MT16.data
endif

if ($number_of_segments_MT17 >1) then
newstack -NumberedFromOne -secs $MT7_start_section-$MT17_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT17.mrcs
echo $number_of_segments_MT17 > $micrograph_stack:r_MT17.data
endif

if ($number_of_segments_MT18 >1) then
newstack -NumberedFromOne -secs $MT18_start_section-$MT18_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT18.mrcs
echo $number_of_segments_MT18 > $micrograph_stack:r_MT18.data
endif

if ($number_of_segments_MT19 >1) then
newstack -NumberedFromOne -secs $MT19_start_section-$MT19_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT19.mrcs
echo $number_of_segments_MT19 > $micrograph_stack:r_MT19.data
endif

if ($number_of_segments_MT20 >1) then
newstack -NumberedFromOne -secs $MT20_start_section-$MT20_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT20.mrcs
echo $number_of_segments_MT20 > $micrograph_stack:r_MT20.data
endif

###

if ($number_of_segments_MT21 >1) then
newstack -NumberedFromOne -secs $MT21_start_section-$MT21_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT21.mrcs
echo $number_of_segments_MT21 > $micrograph_stack:r_MT21.data
endif

if ($number_of_segments_MT22 >1) then
newstack -NumberedFromOne -secs $MT22_start_section-$MT22_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT22.mrcs
echo $number_of_segments_MT22 > $micrograph_stack:r_MT22.data
endif

if ($number_of_segments_MT23 >1) then
newstack -NumberedFromOne -secs $MT23_start_section-$MT23_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT23.mrcs
echo $number_of_segments_MT23 > $micrograph_stack:r_MT23.data
endif

if ($number_of_segments_MT24 >1) then
newstack -NumberedFromOne -secs $MT24_start_section-$MT24_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT24.mrcs
echo $number_of_segments_MT24 > $micrograph_stack:r_MT24.data
endif

if ($number_of_segments_MT25 >1) then
newstack -NumberedFromOne -secs $MT25_start_section-$MT25_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT25.mrcs
echo $number_of_segments_MT25 > $micrograph_stack:r_MT25.data
endif

if ($number_of_segments_MT26 >1) then
newstack -NumberedFromOne -secs $MT26_start_section-$MT26_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT26.mrcs
echo $number_of_segments_MT26 > $micrograph_stack:r_MT26.data
endif

if ($number_of_segments_MT27 >1) then
newstack -NumberedFromOne -secs $MT27_start_section-$MT27_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT27.mrcs
echo $number_of_segments_MT27 > $micrograph_stack:r_MT27.data
endif

if ($number_of_segments_MT28 >1) then
newstack -NumberedFromOne -secs $MT28_start_section-$MT28_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT28.mrcs
echo $number_of_segments_MT28 > $micrograph_stack:r_MT28.data
endif

if ($number_of_segments_MT29 >1) then
newstack -NumberedFromOne -secs $MT29_start_section-$MT29_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT29.mrcs
echo $number_of_segments_MT29 > $micrograph_stack:r_MT29.data
endif

if ($number_of_segments_MT30 >1) then
newstack -NumberedFromOne -secs $MT30_start_section-$MT30_end_section $micrograph_stack:r.mrcs $micrograph_stack:r_MT30.mrcs
echo $number_of_segments_MT30 > $micrograph_stack:r_MT30.data
endif

end

#end of first loop, now working on averaging individual MT stacks

###############################################################

#average_MT_sections for all MTs (MT segment averages are an average of 7 around a central original particle (central particle + 3 either side). ONLY PROCESSES 90 boxes SO FAR.
#Wont be many more (consider diagonal MTs) than only ~100 images max per MT per micrograph (4000pix K2 field of view. To find diagonal (max distance) 4000 x SQRT of 2 = 5656 pixels.7863A at a pixel size of 1.39A/pix. So equivalent of 5560/82 is 95 x 82A repeats.

foreach MT_stack (`ls -1 ./*MT??.mrcs`)
echo 'working on' $MT_stack

set number_of_segments = `grep '' $MT_stack:r.data`
echo 'number of segments in MT =' $number_of_segments

##########

if ($number_of_segments >0) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_00.mrc
0,3
EOF
endif

if ($number_of_segments >1) then 
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_01.mrc
0,4
EOF
endif

if ($number_of_segments >2) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_02.mrc
0,5
EOF
endif

if ($number_of_segments >3) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_03.mrc
0,6
EOF
endif

if ($number_of_segments >4) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_04.mrc
1,7
EOF
endif

if ($number_of_segments >5) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_05.mrc
2,8
EOF
endif

if ($number_of_segments >6) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_06.mrc
3,9
EOF
endif

if ($number_of_segments >7) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_07.mrc
4,10
EOF
endif

if ($number_of_segments >8) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_08.mrc
5,11
EOF
endif

if ($number_of_segments >9) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_09.mrc
6,12
EOF
endif

if ($number_of_segments >10) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_10.mrc
7,13
EOF
endif

###### up to 10

if ($number_of_segments >11) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_11.mrc
8,14
EOF
endif

if ($number_of_segments >12) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_12.mrc
9,15
EOF
endif

if ($number_of_segments >13) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_13.mrc
10,16
EOF
endif

if ($number_of_segments >14) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_14.mrc
11,17
EOF
endif

if ($number_of_segments >15) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_15.mrc
12,18
EOF
endif

if ($number_of_segments >16) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_16.mrc
13,19
EOF
endif

if ($number_of_segments >17) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_17.mrc
14,20
EOF
endif

if ($number_of_segments >18) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_18.mrc
15,21
EOF
endif

if ($number_of_segments >19) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_19.mrc
16,22
EOF
endif

if ($number_of_segments >20) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_20.mrc
17,23
EOF
endif

###### up to 20

if ($number_of_segments >21) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_21.mrc
18,24
EOF
endif

if ($number_of_segments >22) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_22.mrc
19,25
EOF
endif

if ($number_of_segments >23) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_23.mrc
20,26
EOF
endif

if ($number_of_segments >24) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_24.mrc
21,27
EOF
endif

if ($number_of_segments >25) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_25.mrc
22,28
EOF
endif

if ($number_of_segments >26) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_26.mrc
23,29
EOF
endif

if ($number_of_segments >27) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_27.mrc
24,30
EOF
endif

if ($number_of_segments >28) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_28.mrc
25,31
EOF
endif

if ($number_of_segments >29) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_29.mrc
26,32
EOF
endif

if ($number_of_segments >30) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_30.mrc
27,33
EOF
endif

###### up to 30

if ($number_of_segments >31) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_31.mrc
28,34
EOF
endif

if ($number_of_segments >32) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_32.mrc
29,35
EOF
endif

if ($number_of_segments >33) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_33.mrc
30,36
EOF
endif

if ($number_of_segments >34) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_34.mrc
31,37
EOF
endif

if ($number_of_segments >35) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_35.mrc
32,38
EOF
endif

if ($number_of_segments >36) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_36.mrc
33,39
EOF
endif

if ($number_of_segments >37) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_37.mrc
34,40
EOF
endif

if ($number_of_segments >38) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_38.mrc
35,41
EOF
endif

if ($number_of_segments >39) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_39.mrc
36,42
EOF
endif

if ($number_of_segments >40) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_40.mrc
37,43
EOF
endif

###### up to 40

if ($number_of_segments >41) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_41.mrc
38,44
EOF
endif

if ($number_of_segments >42) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_42.mrc
39,45
EOF
endif

if ($number_of_segments >43) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_43.mrc
40,46
EOF
endif

if ($number_of_segments >44) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_44.mrc
41,47
EOF
endif

if ($number_of_segments >45) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_45.mrc
42,48
EOF
endif

if ($number_of_segments >46) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_46.mrc
43,49
EOF
endif

if ($number_of_segments >47) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_47.mrc
44,50
EOF
endif

if ($number_of_segments >48) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_48.mrc
45,51
EOF
endif

if ($number_of_segments >49) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_49.mrc
46,52
EOF
endif

if ($number_of_segments >50) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_50.mrc
47,53
EOF
endif

###### up to 50

if ($number_of_segments >51) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_51.mrc
48,54
EOF
endif

if ($number_of_segments >52) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_52.mrc
49,55
EOF
endif

if ($number_of_segments >53) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_53.mrc
50,56
EOF
endif

if ($number_of_segments >54) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_54.mrc
51,57
EOF
endif

if ($number_of_segments >55) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_55.mrc
52,58
EOF
endif

if ($number_of_segments >56) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_56.mrc
53,59
EOF
endif

if ($number_of_segments >57) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_57.mrc
54,60
EOF
endif

if ($number_of_segments >58) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_58.mrc
55,61
EOF
endif

if ($number_of_segments >59) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_59.mrc
56,62
EOF
endif

if ($number_of_segments >60) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_60.mrc
57,63
EOF
endif

####### up to 60

if ($number_of_segments >61) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_61.mrc
58,64
EOF
endif

if ($number_of_segments >62) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_62.mrc
59,65
EOF
endif

if ($number_of_segments >63) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_63.mrc
60,66
EOF
endif

if ($number_of_segments >64) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_64.mrc
61,67
EOF
endif

if ($number_of_segments >65) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_65.mrc
62,68
EOF
endif

if ($number_of_segments >66) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_66.mrc
63,69
EOF
endif

if ($number_of_segments >67) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_67.mrc
64,70
EOF
endif

if ($number_of_segments >68) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_68.mrc
65,71
EOF
endif

if ($number_of_segments >69) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_69.mrc
66,72
EOF
endif

if ($number_of_segments >70) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_70.mrc
67,73
EOF
endif

####### up to 70

if ($number_of_segments >71) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_71.mrc
68,74
EOF
endif

if ($number_of_segments >72) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_72.mrc
69,75
EOF
endif

if ($number_of_segments >73) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_73.mrc
70,76
EOF
endif

if ($number_of_segments >74) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_74.mrc
71,77
EOF
endif

if ($number_of_segments >75) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_75.mrc
72,78
EOF
endif

if ($number_of_segments >76) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_76.mrc
73,79
EOF
endif

if ($number_of_segments >77) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_77.mrc
74,80
EOF
endif

if ($number_of_segments >78) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_78.mrc
75,81
EOF
endif

if ($number_of_segments >79) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_79.mrc
76,82
EOF
endif

if ($number_of_segments >80) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_80.mrc
77,83
EOF
endif

####### up to 80

if ($number_of_segments >81) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_81.mrc
78,84
EOF
endif

if ($number_of_segments >82) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_82.mrc
79,85
EOF
endif

if ($number_of_segments >83) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_83.mrc
80,86
EOF
endif

if ($number_of_segments >84) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_84.mrc
81,87
EOF
endif

if ($number_of_segments >85) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_85.mrc
82,88
EOF
endif

if ($number_of_segments >86) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_86.mrc
83,89
EOF
endif

if ($number_of_segments >87) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_87.mrc
84,90
EOF
endif

if ($number_of_segments >88) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_88.mrc
85,91
EOF
endif

if ($number_of_segments >89) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_89.mrc
86,92
EOF
endif

if ($number_of_segments >90) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_90.mrc
87,93
EOF
endif

####### up to 90

if ($number_of_segments >91) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_91.mrc
88,94
EOF
endif

if ($number_of_segments >92) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_92.mrc
89,95
EOF
endif

if ($number_of_segments >93) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_93.mrc
90,96
EOF
endif

if ($number_of_segments >94) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_94.mrc
91,97
EOF
endif

if ($number_of_segments >95) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_95.mrc
92,98
EOF
endif

if ($number_of_segments >96) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_96.mrc
93,99
EOF
endif

if ($number_of_segments >97) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_97.mrc
94,100
EOF
endif

if ($number_of_segments >98) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_98.mrc
95,101
EOF
endif

if ($number_of_segments >99) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_99.mrc
96,102
EOF
endif

if ($number_of_segments >100) then
avgstack << EOF
$MT_stack:r.mrcs
$MT_stack:r_SA_100.mrc
97,103
EOF
endif

####### up to 100
echo 'finished making MT segment averages for' $number_of_segments 'MT segments in' $MT_stack

###############################################################################

###STACK SEGMENT AVERAGES INTO MTS

newstack $MT_stack:r_SA_*.mrc $MT_stack:r_SA_stack.mrc

###NORMALIZE??? 

bimg -images -rescale 0,1 $MT_stack:r_SA_stack.mrc $MT_stack:r_SA_stack_norm.mrc
mv $MT_stack:r_SA_stack_norm.mrc $MT_stack:r_SA_stack_norm.st

#CLEAN UP TEMPORARY FILES. BIT DANGEROUS BUT SHOULD WORK. COMMENT OUT IF YOU WANT ALL INTERMEDIATE FILES.

rm -rf $MT_stack
rm -rf $MT_stack:r.data
rm -rf $MT_stack:r_SA_*.mrc

end

#############################################################################
######stack all MTs into segment average stacks for each micrograph

foreach micrograph (*.mrcs)
echo 'making micrograph SAs stack for' $micrograph

newstack $micrograph:r_MT*_SA_stack_norm.st $micrograph:r_SAs.mrcs

end

#############################################################################

#############################################################################
#Normalise final stacks with relion

foreach final_stack (*_SAs.mrcs)

echo 'normalising in relion' $final_stack

relion_preprocess --norm true --bg_radius $background_box_radius --operate_on $final_stack --operate_out $final_stack:r_norm.mrcs

rm -rf $final_stack
rm -rf $final_stack:r_norm.star

end

rename _SAs_norm.mrcs.mrcs _SAs.mrcs *_SAs_norm.mrcs.mrcs


#############################################################################

#CLEAN UP TEMPORARY FILES. BIT DANGEROUS BUT SHOULD WORK. COMMENT OUT IF YOU WANT ALL INTERMEDIATE FILES.

rm -rf *.st

###############################################################################

