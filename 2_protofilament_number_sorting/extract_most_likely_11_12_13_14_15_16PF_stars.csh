#!/bin/csh -f

#Joe Atherton 30/01/19#

#run as source extract_most_likely_11_12_13_14_15_16PF_stars.csh run_it000_data.star run_it001_data.star run_it001_data_corrected01.star where run_it000_data.star is the star file input into 3D class, run_it001_data.star is the output from the 3D class and run_it001_data_corrected01.star is the output file, with classes corrected to be the most common for each MT
#MUST HAVE A COLUMN $34 WHICH IS NEW CORRECTED CLASS NUMBER.
#extracts individual star files for different PF numbers, and resets all columns except 3D class. NormCorrection LogLikeliContribution MaxValueProbDistribution and NrOfSignificant samples.

set star_file=$1
set class_3D_star_file=$2
set corrected_class_3D_star_file=$3

#Merge original star file and new star file;

echo | grep '.mrc' $star_file > $star_file:r_temp1.star
echo | grep '.mrc' $corrected_class_3D_star_file > $corrected_class_3D_star_file:r_temp1.star

#take all columns from unrefined original .star file except 3D class;
echo | awk '{printf("%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28)}' $star_file:r_temp1.star >> $star_file:r_temp2.star

#take last 4 columns from refined new .star file, 3D class and norm correlation;
echo | awk '{printf(" %d\t%.6f\t%.6f\t%.6f\t%d\n", $34, $30, $31, $32, $33)}' $corrected_class_3D_star_file:r_temp1.star >> $class_3D_star_file:r_temp2.star

#paste the columns together;
echo | paste $star_file:r_temp2.star $class_3D_star_file:r_temp2.star >> $class_3D_star_file:r_temp3.star



#################


#copy headers to new 11pf class file
echo ' ' > $class_3D_star_file:r_11pf_class.star
echo 'data_images' >> $class_3D_star_file:r_11pf_class.star
echo ' ' >> $class_3D_star_file:r_11pf_class.star
echo 'loop_' >> $class_3D_star_file:r_11pf_class.star
#copy column headers to new 11pf class files
echo | grep '_rln*' $class_3D_star_file >> $class_3D_star_file:r_11pf_class.star
#add extra column of unique MT number

#copy headers to new 12pf class file
echo ' ' > $class_3D_star_file:r_12pf_class.star
echo 'data_images' >> $class_3D_star_file:r_12pf_class.star
echo ' ' >> $class_3D_star_file:r_12pf_class.star
echo 'loop_' >> $class_3D_star_file:r_12pf_class.star
#copy column headers to new 12pf class files
echo | grep '_rln*' $class_3D_star_file >> $class_3D_star_file:r_12pf_class.star

#copy headers to new 13pf class file
echo ' ' > $class_3D_star_file:r_13pf_class.star
echo 'data_images' >> $class_3D_star_file:r_13pf_class.star
echo ' ' >> $class_3D_star_file:r_13pf_class.star
echo 'loop_' >> $class_3D_star_file:r_13pf_class.star
#copy column headers to new 13pf class files
echo | grep '_rln*' $class_3D_star_file >> $class_3D_star_file:r_13pf_class.star

#copy headers to new 14pf class file
echo ' ' > $class_3D_star_file:r_14pf_class.star
echo 'data_images' >> $class_3D_star_file:r_14pf_class.star
echo ' ' >> $class_3D_star_file:r_14pf_class.star
echo 'loop_' >> $class_3D_star_file:r_14pf_class.star
#copy column headers to new 14pf class files
echo | grep '_rln*' $class_3D_star_file >> $class_3D_star_file:r_14pf_class.star

#copy headers to new 15pf class file
echo ' ' > $class_3D_star_file:r_15pf_class.star
echo 'data_images' >> $class_3D_star_file:r_15pf_class.star
echo ' ' >> $class_3D_star_file:r_15pf_class.star
echo 'loop_' >> $class_3D_star_file:r_15pf_class.star
#copy column headers to new 15pf class files
echo | grep '_rln*' $class_3D_star_file >> $class_3D_star_file:r_15pf_class.star

#copy headers to new 16pf class file
echo ' ' > $class_3D_star_file:r_16pf_class.star
echo 'data_images' >> $class_3D_star_file:r_16pf_class.star
echo ' ' >> $class_3D_star_file:r_16pf_class.star
echo 'loop_' >> $class_3D_star_file:r_16pf_class.star
#copy column headers to new 16pf class files
echo | grep '_rln*' $class_3D_star_file >> $class_3D_star_file:r_16pf_class.star

################################################

#extract column data for PF classes to temp files

echo | grep '	 1	0*' $class_3D_star_file:r_temp3.star >> $class_3D_star_file:r_11pf_class.star
echo | grep '	 2	0*' $class_3D_star_file:r_temp3.star >> $class_3D_star_file:r_12pf_class.star
echo | grep '	 3	0*' $class_3D_star_file:r_temp3.star >> $class_3D_star_file:r_13pf_class.star
echo | grep '	 4	0*' $class_3D_star_file:r_temp3.star >> $class_3D_star_file:r_14pf_class.star
echo | grep '	 5	0*' $class_3D_star_file:r_temp3.star >> $class_3D_star_file:r_15pf_class.star
echo | grep '	 6	0*' $class_3D_star_file:r_temp3.star >> $class_3D_star_file:r_16pf_class.star

echo 'star files generated for different PF numbers;'
echo $class_3D_star_file:r_11pf_class.star
echo | grep '.mrc' $class_3D_star_file:r_11pf_class.star | wc -l
echo 'particles/segments'
echo $class_3D_star_file:r_12pf_class.star
echo | grep '.mrc' $class_3D_star_file:r_12pf_class.star | wc -l
echo 'particles/segments'
echo $class_3D_star_file:r_13pf_class.star
echo | grep '.mrc' $class_3D_star_file:r_13pf_class.star | wc -l
echo 'particles/segments'
echo $class_3D_star_file:r_14pf_class.star
echo | grep '.mrc' $class_3D_star_file:r_14pf_class.star | wc -l
echo 'particles/segments'
echo $class_3D_star_file:r_15pf_class.star
echo | grep '.mrc' $class_3D_star_file:r_15pf_class.star | wc -l
echo 'particles/segments'
echo $class_3D_star_file:r_16pf_class.star
echo | grep '.mrc' $class_3D_star_file:r_16pf_class.star | wc -l
echo 'particles/segments'

##################################################################################

#remove temopory files. comment this out for diagnostics on temporary files.

rm -rf *temp*.star

##################################################################################
