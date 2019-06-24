#!/bin/csh -f

#Joe Atherton 30/01/19#
#rescales helical track length to correct value.
#run as source scale_helical_track_length.csh particles.star 0.25 where particles.star is the star file to be converted and 0.25 is the upscale factor for helical track length (from binx1 to binx4 data).

set star_file=$1
set upscale_factor=$2

#copy headers to new file
echo ' ' > $star_file:r_scaled_helical_track_length.star
echo 'data_images' >> $star_file:r_scaled_helical_track_length.star
echo ' ' >> $star_file:r_scaled_helical_track_length.star
echo 'loop_' >> $star_file:r_scaled_helical_track_length.star

#copy column headers to new file
echo | grep '_rln*' $star_file >> $star_file:r_scaled_helical_track_length.star

#extract column data to temp file
echo | grep '.mrc' $star_file > $star_file:r_temp1.star


#print all data to new file with modifications
echo | awk '{printf("%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n", $1, $2, $3, $4, $5, $6*"'$upscale_factor'", $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22)}' $star_file:r_temp1.star >> $star_file:r_scaled_helical_track_length.star

rm -rf $star_file:r_temp1.star

















