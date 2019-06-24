#!/bin/csh -f

#run as source extract_unified_seam_classes_star_relionv3.csh run_it000_data.star run_it001_data.star run_it001_data_corrected01.star where run_it000_data.star is the star file input into 3D class, run_it001_data.star is the output from the 3D class and run_it001_data_corrected01.star is the output file, with classes corrected to be the most common for each MT
#MUST HAVE A COLUMN $35 WHICH IS NEW CORRECTED CLASS NUMBER.
#extracts individual star files for seam locations

set star_file=$1
set class_3D_star_file=$2
set corrected_class_3D_star_file=$3

#Merge original star file and new star file;

echo | grep '.mrc' $star_file > $star_file:r_temp1.star
echo | grep '.mrc' $corrected_class_3D_star_file > $corrected_class_3D_star_file:r_temp1.star

#take all columns from unrefined original .star file except 3D class;
echo | awk '{printf("%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28)}' $star_file:r_temp1.star >> $star_file:r_temp2.star

#take last 5 columns from refined new .star file, 3D class and norm correlation;
echo | awk '{printf(" %d\t%.6f\t%.6f\t%.6f\t%d\n", $35, $30, $31, $32, $33)}' $corrected_class_3D_star_file:r_temp1.star >> $class_3D_star_file:r_temp2.star

#paste the columns together;
echo | paste $star_file:r_temp2.star $class_3D_star_file:r_temp2.star >> $class_3D_star_file:r_temp3.star


#################


#copy headers to new unified class file
echo ' ' > $class_3D_star_file:r_unified_class.star
echo 'data_images' >> $class_3D_star_file:r_unified_class.star
echo ' ' >> $class_3D_star_file:r_unified_class.star
echo 'loop_' >> $class_3D_star_file:r_unified_class.star
#copy column headers to new unified class file
echo | grep '_rln*' $class_3D_star_file >> $class_3D_star_file:r_unified_class.star
#add extra column of unique MT number


################################################

#extract column data for SEAM classes to files

echo | grep '.mrc' $class_3D_star_file:r_temp3.star >> $class_3D_star_file:r_unified_class.star

##################################################################################

#remove temopory files. comment this out for diagnostics on temporary files.

rm -rf *temp*.star

##################################################################################
