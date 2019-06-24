# MiRP
Microtubule Image Processing in RELION Pipeline

#### INSTRUCTIONS FOR MiRP v2.2 PROCESSING OF MTS #### Joe Atherton/Alex Cook 2019 #

N.B YOU WILL FIND ALL REQUIRED SCRIPTS IN THE âSCRIPTSâ FOLDER AND EXAMPLE REFERENCES IN THE âREFERENCESâ FOLDER.

ONLY TESTED USING RELION v3.0/BETA4- MAY VARY WITH VERSION AS SENSITIVE TO COLUMN NUMBERS.

--------------------------------------------------------------------------------------------------------------------

PREPROCESSING

--------------------------------------------------------------------------------------------------------------------

1. Import micrographs (FULL DOSE)

2. Run gCTF using GPU on FULL DOSE micrographs. Adjust parameters according to your dataset. After this point make sure you are using dose-weighted micrographs for all
subsequent steps.

3. Select best micrographs; remove those with poor power spectra (indicating perhaps poor motion correction) or ice contamination.

4. Use relion's manual picker (picking with start-end coordinates for helicies) to pick straight and undamaged regions of MTs.
		
5. Extract 4x binned particles from these coordinates in helical mode;
		Input Coordinates: *coords_suffix_manualpick.star file in manual pick job*
		OR re-extract refined particles? No
		Refined particles star file: N/A
		Re-center refined coordinates: No
		Particle box size (pix) = ~600A (432 pixels if 1.39A/pixel)
		Rescale particles? Yes
		Re-scaled size (pixels): *Original particle box size/4*
		
		Extract helical segments = yes
		Tube Diameter (A) = 400
		Use bimodal angular priors? yes
		Coordinates are start end only? yes
		Cut helical tubes into segments? Yes
		Number of asymmetrical units: 1
		Helical rise (A): 82
		
N.B Try using multiple MPI without GPU- e.g 30 MPI 1 Thread.
		
6. Scale helical track length; RELION at present does not scale helical track length column properly when changing particle binning. Therefore you need to run the
script scale_helical_track_length.csh on the binnedx4 extracted particles.star file to scale the track length (source scale_helical_track_length.csh particles.star 0.25).
	
7. Make rough binnedx4 MT segment averages;

	- Go to directory in the ./Extract/*micrographs* where there are individual stacks (.mrcs) and corresponding extract.star files.
	- Run the script preprocess_segment_averages.csh (source preprocess_segment_averages.csh)
	- Change your extracted particles.star to have links to the segment averages rather than regular partilces. Most easily done with opening the star file in nedit and using 'find and replace' function. Be careful to only change text you want to change!

------------------------------------------------------------------------------------------------------------------------

PROTOFILAMENT NUMBER SORTING

------------------------------------------------------------------------------------------------------------------------

	
8. Run supervised 3D class on binx4 segment averages to sort PF numbers. You will need 3D references in a directory for different PF numbers, and a star file with a list of these
references and their paths (see example in âreferencesâ directory). No need for a reference mask as one iteration only;
		
		I/O tab:
		
		Input images STAR file: *4x binned extracted segment averages star file with scaled track length* (MAKE SURE .star file particle paths are to the segment averages)
		Reference map: *path to references .star file*
		Reference mask: None
		
		Reference tab:
		
		Reference map is on absolute greyscale? No
		Initial low pass filter (A): 0 (current references are already 15A low pass filtered- if your references are not already low pass filtered this can be set to ~12-15A)
		Symmetry: C1
		
		CTF tab:
		
		Do CTF-correction? Yes
		Has reference been CTF-corrected? Yes
		Have the data been phase flipped: No
		Ignore CTFs until first peak: No
		
		Optimisation Tab:
		
		Number of classes: 6 (*or equal to number of references*)
		Reg param T: 4
		Number of iterations: 1 (output reconstuctions not good enough to be references for next iteration and 1 iteration is enough)
		Use fast subsets (for large datasets): NO (we want to use all the data as we only do 1 iteration)
		Mask diameter (A): 580 (keeping big for time being, as most signal required for classification and seam finding)
		Mask individual particles with zeros? Yes
		Limit resolution to E-step to (A): 12A (binned data Nyquist is ~11A anyway and default references are 15A low-pass filtered).
		
		Sampling tab:
		
		Perform image alignment? Yes
		Angular sampling interval: 1.8 deg (quick and accurate enough for time being)
		Offset search range (pix): *~70A* (manually picked boxes can be quite far off)
		Offset search step (pix) 1
		Perform local angular searches? No
		
		Helix tab:
		
		Do helical reconstruction? Yes
		Tube diameter - inner, outer (A): -1 400 (no inner mask; doesnt matter at this stage because only 1 iteration in 3D using synthetic references with no noise).
		Angular search range - tilt, psi (deg): 15 10 (fine as priors of tilt (90) and psi (from start and end coords) are close enough).
		Keep tilt-prior fixed: No (we want to test both polarities)
		Range factor of local averaging: -1
		Apply helical symmetry: No
		Do local searches of symmetry: No
		Additional arguments: --dont_check_norm
		
		N.B Fastest way run 3D classification with alignment is with GPU. Also see notes 5/6 in the âData Optimisationâ section below.
		
9. Unify the 3D classes, so that the most common (mode) 3D class found for particles from a particular MT are assigned to the whole MT.
		- Go to folder with the 3D classification iteration 1 results.
		- Run the perl script to reformat the .star file for R; (e.g ./reformat.pl<./run_it001_data.star>./run_it001_data_corrected00.star)
		- Type 'R' to start R (must have R available on the computer).
		- Type the following line by line;
			x= read.table ("run_it001_data_corrected00.star", header= T)
			library(dplyr)
			getmode <- function(v){
			uniqv <- unique(v)
			uniqv [which.max (table(match( v,uniqv))) ]
			}
			modeDataset <- x %>% group_by(GroupNumber, HelicalTubeID)%>% summarise(ClassNumber2=getmode(ClassNumber))
			y<- left_join(x,modeDataset, by=c("GroupNumber","HelicalTubeID"))
			write.table(y,file="run_it001_data_corrected01.star", quote=F,row=F,sep="\t")
			
		- Type 'q()' to quit R.
		- Output file is called run_it001_data_corrected01.star
		- Run the script 'extract_most_likely_11_12_13_14_15_16PF_stars.csh' to extract	star files for each PF number (11-16pfs). E.g source extract_most_likely_11_12_13_14_15_16PF_stars.csh run_it000_data.star run_it001_data.star run_it001_data_corrected01.star. 
		This resets all shifts and angles to priors.
		
p.s if loading library fails, you may need to install the library; If so type;
install.packages("dplyr")
			
From this point, you will be working with .star files belonging to a particlar 3D class. The example helical parameters used below are for a 13pf MT.

------------------------------------------------------------------------------------------------------------------------

GLOBAL SEARCH

------------------------------------------------------------------------------------------------------------------------


10. 1st Refine3D, roughly align particles to get PSI and THETA (TILT);
	
		I/O tab:
		
		Input images STAR file: *output 3D classification .star file belonging to particlar PF number class*
		Reference map: *4x binned decorated synthetic reference of corresponding PF number*
		Reference mask: None
		
		Reference tab:
		
		Reference map is on absolute greyscale? No
		Initial low pass filter (A): 0 (current references are already 15A low pass filtered- if your references are not already low pass filtered this can be set to ~12-15A)
		Symmetry: C1
		
		CTF tab:
		
		Do CTF-correction? Yes
		Has reference been CTF-corrected? Yes
		Have the data been phase flipped: No
		Ignore CTFs until first peak: No
		
		Optimisation Tab:
		
		Mask diameter (A): 580 (keeping big for time being, as most signal required for classification and seam finding)
		Mask individual particles with zeros? Yes
		Use solvent-flatened FSCs? Yes
		
		Sampling tab:
		
		Angular sampling interval: 1.8 deg (quick and accurate enough for time being)
		Offset search range (pix): *~45A* (adjust for pixel size) 
		Offset search step (pix) 1
		Local searches from auto-sampling: 0.9 deg 
		
		Helix tab:
		
		Do helical reconstruction? Yes
		Tube diameter - inner, outer (A): -1 400 (no inner mask; doesnt matter at this stage because only 1 iteration in 3D using synthetic references with no noise).
		Angular search range - tilt, psi (deg): 15 10 (fine as priors of tilt (90) and psi (from start and end coords) are close enough).
		Keep tilt-prior fixed: No
		Range factor of local averaging: 2 (helps smooth angles over an MT)
		Apply helical symmetry: No
		Number of asymmetrical units: 13 (i.e for 13pf kinesin-decorated reconstructions, 12 AUs for proteins that don't bind at the seam like DCX and CAMSAP)
		Intial twist (deg). rise (A): -27.6689 9.46346 (*example is 13pf MT roughly, modify accordingly*)
		Central Z length (%): 30
		Do local searches of symmetry: No
		Twist search - Min,Max,Step (deg): -27 -28 0.1 (*example is 13pf MT roughly, modify accordingly*)
	   	Rise search - Min,Max,Step (A): 9.4 9.7 0.1 (*example is 13pf MT roughly, modify accordingly*)
		
		Additional arguments: --dont_check_norm	--ignore_helical_symmetry --iter 1 (I think the option of 1 iteration does not work in auto-refine, so you will have to stop the job manually for the mo after iteration 1).
	
	
11. Reset the PSI and THETA (TILT) angles in the resulting iteration 1 star file to the priors and reset the PHI (ROT) and the shifts to 0, using the 'relion_reset_angles_to_priors_and_shifts_to_0_relion3.csh' script (E.g source relion_reset_angles_to_priors_and_shifts_to_0_relion3.csh run_it001_data.star)			

	
12. 2nd Refine3D, finer align particles and get initial PHI (Rot);

	
		I/O tab:
		
		Input images STAR file: *output .star file created above with reset angles and shifts* (still with links to segment averages)
		Reference map: *4x binned decorated synthetic reference of corresponding PF number*
		Reference mask: None
		
		Reference tab:
		
		Reference map is on absolute greyscale? No
		Initial low pass filter (A): 0 (current references are already 15A low pass filtered- if your references are not already low pass filtered this can be set to ~12-15A)
		Symmetry: C1
		
		CTF tab:
		
		Do CTF-correction? Yes
		Has reference been CTF-corrected? Yes
		Have the data been phase flipped: No
		Ignore CTFs until first peak: No
		
		Optimisation Tab:
		
		Mask diameter (A): 580 (keeping big for time being, as most signal required for classification and seam finding)
		Mask individual particles with zeros? Yes
		Use solvent-flatened FSCs? Yes
		
		Sampling tab:
		
		Angular sampling interval: 0.9 deg (finer now)
		Offset search range (pix): *~45A* (adjust for pixel size)
		Offset search step (pix) 0.5 (finer works a bit better)
		Local searches from auto-sampling: 0.5 deg 
		
		Helix tab:
		
		Do helical reconstruction? Yes
		Tube diameter - inner, outer (A): -1 400 (no inner mask; doesnt matter at this stage because only 1 iteration in 3D using synthetic references with no noise).
		Angular search range - tilt, psi (deg): 3 3 (small search range as we are close enough from last refinement round).
		Keep tilt-prior fixed: No
		Range factor of local averaging: 2 (helps smooth angles over an MT)
		Apply helical symmetry: No
		Number of asymmetrical units: 13 (i.e for 13pf kinesin-decorated reconstructions, 12 AUs for proteins that don't bind at the seam like DCX and CAMSAP)
		Intial twist (deg). rise (A): -27.6689 9.46346 (*example is 13pf MT roughly, modify accordingly*)
		Central Z length (%): 30
		Do local searches of symmetry: No
		Twist search - Min,Max,Step (deg): -27 -28 0.1 (*example is 13pf MT roughly, modify accordingly*)
	   	Rise search - Min,Max,Step (A): 9.4 9.7 0.1 (*example is 13pf MT roughly, modify accordingly*)
		
		Additional arguments: --dont_check_norm	--ignore_helical_symmetry --iter 1 (I think the option of 1 iteration does not work in auto-refine, so you will have to stop the job manually for the mo after iteration 1).
	
		N.B Best to run this on GPU but as not a massive search range you could try multiple MPIs also (e.g 20 MPIs 1 Thread on diana, minerva or sulis).		

------------------------------------------------------------------------------------------------------------------------

INITIAL SEAM ASSIGNMENT

------------------------------------------------------------------------------------------------------------------------

13. Unify PHI angles. 
		- Go to the directory with the resulting iteration 1 star file.
		- First load the right software: module load anaconda/v2-5.2.0 (or equivalent in your institution).
		- RUN THE SCRIPT mirpy.py on ouput star file from 2nd 3Drefine iteration 1 (python mirpy.py -ang run_it001_data.star -id rlnAngleRot)
		- Reset shifts to 0 and PSI/THETA angles to priors, keeping the unified PHI (ROT) angles using the script relion_keep_PHI_reset_angles_to_priors_and_shifts_to_0_relion3.csh (e.g source relion_keep_PHI_reset_angles_to_priors_and_shifts_to_0_relion3.csh fitted_angles.star).

The next step requires raw particles, not segment averages. Therefore copy the output star file above *.star to a new name *_NOT_SAs.star, then nedit the new .star file to have links to original particles rather than segment averages (can use replace funtion in
nedit, find '_SAs.mrcs' and replace with '.mrcs'.

14. 3rd Refine3D, get fine shifts whilst keeping the PHI (Rot) fixed around the previous integer minimum;
	
		I/O tab:
		
		Input images STAR file: *output .star file created above with links to raw particles (NOT segment averages)*
		Reference map: *4x binned decorated synthetic reference of corresponding PF number*
		Reference mask: None
		
		Reference tab:
		
		Reference map is on absolute greyscale? No
		Initial low pass filter (A): 0 (current references are already 15A low pass filtered- if your references are not already low pass filtered this can be set to ~12-15A)
		Symmetry: C1
		
		CTF tab:
		
		Do CTF-correction? Yes
		Has reference been CTF-corrected? Yes
		Have the data been phase flipped: No
		Ignore CTFs until first peak: No
		
		Optimisation Tab:
		
		Mask diameter (A): 580 (keeping big for time being, as most signal required for classification and seam finding)
		Mask individual particles with zeros? Yes
		Use solvent-flatened FSCs? Yes
		
		Sampling tab:
		
		Angular sampling interval: 0.9 deg (finer now)
		Offset search range (pix): *~45A* (adjust for pixel size)
		Offset search step (pix) 0.25 (IMPORTANT: FINE SAMPLING IS REQUIRED FOR NEXT STEP.. FOR THE NEXT SCRIPT TO WORK PROPERLY)
		Local searches from auto-sampling: 0.5 deg 
		
		Helix tab:
		
		Do helical reconstruction? Yes
		Tube diameter - inner, outer (A): -1 400 (no inner mask; doesnt matter at this stage because only 1 iteration in 3D using synthetic references with no noise).
		Angular search range - tilt, psi (deg): 3 3 (small search range as we are close enough from last refinement round).
		Keep tilt-prior fixed: No
		Range factor of local averaging: 2 (helps smooth angles over an MT)
		Apply helical symmetry: No
		Number of asymmetrical units: 13 (i.e for 13pf kinesin-decorated reconstructions, 12 AUs for proteins that don't bind at the seam like DCX and CAMSAP)
		Intial twist (deg). rise (A): -27.6689 9.46346 (*example is 13pf MT roughly, modify accordingly*)
		Central Z length (%): 30
		Do local searches of symmetry: No
		Twist search - Min,Max,Step (deg): -27 -28 0.1 (*example is 13pf MT roughly, modify accordingly*)
	   	Rise search - Min,Max,Step (A): 9.4 9.7 0.1 (*example is 13pf MT roughly, modify accordingly*)
		
		Additional arguments: --dont_check_norm --iter 1  --sigma_rot 3 --ignore_helical_symmetry (IMPORTANT THAT SIGMA ROT IS 3, RESTRAINS PHI ANGLES TO REMAIN IN INTEGER
		MINIMUM, also I think the option of 1 iteration does not work in auto-refine, so you will have to stop the job manually for the mo after iteration 1).
	
		N.B Best to run this on GPU but as not a massive search range; you could try multiple CPU MPIs also (e.g 20 MPIs 1 Thread). Also see note 5 in the âData Optimisationâ section below.	
		
------------------------------------------------------------------------------------------------------------------------

X/Y SHIFT SMOOTHING

------------------------------------------------------------------------------------------------------------------------
		
15. Smooth XY shifts to remove overlapping boxes and keep individual boxes ~AU apart.

- Go to the directory with the resulting iteration 1 star file from 3rd Refine3D.
		- First load the right software: module load anaconda/v2-5.2.0 (or equivalent in your institution).
		- RUN THE SCRIPT mirpy.py on ouput star file from 3rd 3Drefine iteration 1 (python mirpy.py -xy run_it001_data.star). Don't worry about the error messages they don't seem to be a problem. The resulting file will be called uniXY_data.star.
		
16. 4th Refine3D, get final shifts whilst keeping the PHI (Rot) fixed to the previous integer minimum;
	
		I/O tab:
		
		Input images STAR file: output .star file created above (uniXY_data.star) (NOT segment averages)*
		Reference map: *4x binned decorated synthetic reference of corresponding PF number*
		Reference mask: None
		
		Reference tab:
		
		Reference map is on absolute greyscale? No
		Initial low pass filter (A): 0 (current references are already 15A low pass filtered- if your references are not already low pass filtered this can be set to ~12-15A)
		Symmetry: C1
		
		CTF tab:
		
		Do CTF-correction? Yes
		Has reference been CTF-corrected? Yes
		Have the data been phase flipped: No
		Ignore CTFs until first peak: No
		
		Optimisation Tab:
		
		Mask diameter (A): 580 (keeping big for time being, as most signal required for classification and seam finding)
		Mask individual particles with zeros? Yes
		Use solvent-flatened FSCs? Yes
		
		Sampling tab:
		
		Angular sampling interval: 0.9 deg
		Offset search range (pix): *~35A* (less than monomer repeat distance to stp boxes jumping back into each other, adjust for pixel size).
		Offset search step (pix) 0.5
		Local searches from auto-sampling: 0.5 deg 
		
		Helix tab:
		
		Do helical reconstruction? Yes
		Tube diameter - inner, outer (A): -1 400 (no inner mask; doesnt matter at this stage because only 1 iteration in 3D using synthetic references with no noise).
		Angular search range - tilt, psi (deg): 2 2 (smaller search range as we are close enough from last refinement round).
		Keep tilt-prior fixed: No
		Range factor of local averaging: 2 (helps smooth angles over an MT)
		Apply helical symmetry: No
		Number of asymmetrical units: 13 (i.e for 13pf kinesin-decorated reconstructions, 12 AUs for proteins that don't bind at the seam like DCX and CAMSAP)
		Intial twist (deg). rise (A): -27.6689 9.46346 (*example is 13pf MT roughly, modify accordingly*)
		Central Z length (%): 30
		Do local searches of symmetry: No
		Twist search - Min,Max,Step (deg): -27 -28 0.1 (*example is 13pf MT roughly, modify accordingly*)
	   	Rise search - Min,Max,Step (A): 9.4 9.7 0.1 (*example is 13pf MT roughly, modify accordingly*)
		
		Additional arguments: --dont_check_norm --iter 1  --sigma_rot 3 --ignore_helical_symmetry (IMPORTANT THAT SIGMA ROT IS 3, RESTRAINS PHI ANGLES TO REMAIN IN INTEGER
		MINIMUM, also I think the option of 1 iteration does not work in auto-refine, so you will have to stop the job manually for the mo after iteration 1)
	
		N.B Best to run this on GPU but as not a massive search range you could try multiple MPI CPUs also (e.g 20 MPIs 1 Thread).

------------------------------------------------------------------------------------------------------------------------

REFINED SEGMENT AVERAGES

------------------------------------------------------------------------------------------------------------------------

17.	Rextract centered 4x binned particles from these coordinates in helical mode;
		Input Coordinates: N/A (leave empty)
		OR re-extract refined particles? YES
		Refined particles star file: *your refined particle star file created in the last step*
		Reset the refined offsets to 0? No
		Re-center refined coordinates: Yes
		Recenter on- X,Y,Z (pix): 0,0,0
		Particle box size (pix) = ~600A (432 pixels if 1.39A/pixel)
		Rescale particles? Yes
		Re-scaled size (pixels): 108 (*Original box size/4*)
		
		Extract helical segments = yes
		Tube Diameter (A) = 400
		Use bimodal angular priors? yes
		Coordinates are start end only? No
		Cut helical tubes into segments? Yes
		Number of asymmetrical units: 1
		Helical rise (A): 82

N.B I've found the fastest way to run extraction jobs is through multiple MPI on a CPU computer- e.g 30 MPI 1 Thread.

18.	Make centered binnedx4 MT segment averages;

	- Go to directory in the ./Extract/*micrographs* where there are new centered individual stacks (.mrcs) and corresponding extract.star files.
	- Run the script divide_stacks_and_make_MT_segment_average_averages.csh (source divide_stacks_and_make_MT_segment_average_averages.csh)
	- Will need to now go back to segment averages for the coming 3D classification seam finding, so edit the extraction job particles.star file to have links to segment averages (_SAs.mrcs) rather than raw particles (.mrcs). You should not have to rescale helical track length as with the last extraction.

------------------------------------------------------------------------------------------------------------------------

SEAM CHECK

------------------------------------------------------------------------------------------------------------------------


19.  Seam finding by 3D classification with references with different seam positions. You will need 3D references of the appropriate PF number with different seam positions and their 41A shifted positions in a directory (26 for 13pfs, 28 for
14pfs), and a star file with a list of these references and their paths. IT IS IMPORTANT THAT THE SEAM POSITIONS IN YOUR .STAR FILE AND CORRESPONDING REFERENCES MATCH THOSE IN THE EXAMPLE .STAR FILE AND REFERENCES (LATER SEAM POSITION CORRECTION WILL NOT WORK OTHERWISE).
Examples can be found for 13pf and 14pf CKK only MTs in the directories 'seam_check/13pf_refs/' and 'seam_check/14pf_refs/' respectively. At this stage I found that using 'decorator' only references (e.g CKK or kinesin only) works better than tubulin+decorator references. No need for a reference mask as one iteration only.

		I/O tab:
		
		Input images STAR file: *4x binned extract.star file from above centered 2nd extraction* (MAKE SURE .star file particle paths are to the segment averages)
		Reference map: *path to references .star file* 
		Reference mask: None
		
		Reference tab:
		
		Reference map is on absolute greyscale? No
		Initial low pass filter (A): 0 (current references are already 15A low pass filtered- if your references are not already low pass filtered this can be set to ~12-15A)
		Symmetry: C1
		
		CTF tab:
		
		Do CTF-correction? Yes
		Has reference been CTF-corrected? Yes
		Have the data been phase flipped: No
		Ignore CTFs until first peak: No
		
		Optimisation Tab:
		
		Number of classes: * equal to number of references* (26 for 13pf, 28 for 14pf)
		Reg param T: 4
		Number of iterations: 1 (output reconstuctions not good enough to be references for next iteration and 1 iteration is enough)
		Use fast subsets (for large datasets): NO (we want to use all the data in the only iteration)
		Mask diameter (A): 580 (keeping big for time being, as most signal required for classification and seam finding)
		Mask individual particles with zeros? Yes
		Limit resolution to E-step to (A): 12A (binned data Nyquist is ~11A anyway and default references are 15A low-pass filtered).
		
		Sampling Tab
		
		Perform image alignment? NO (IMPORTANT)
		Angular sampling interval: N/A
		Offset search range (pix): N/A
		Offset search step (pix) N/A
		Perform local angular searches? N/A
		
		Helix tab:
		
		Do helical reconstruction? Yes
		Tube diameter - inner, outer (A): -1 400 (no inner mask; doesnt matter at this stage because only 1 iteration in 3D using synthetic references with no noise).
		Angular search range - tilt, psi (deg): 3 3 
		Keep tilt-prior fixed: No
		Range factor of local averaging: -1
		Apply helical symmetry: No
		Do local searches of symmetry: No
		Additional arguments: --dont_check_norm --ignore_helical_symmetry
		
		N.B Try 20 threads on a CPU machine. Also see notes 5/6 in the âData Optimisationâ section below.
		
20. Unify seam location for each MT. 
	
		- Go to folder with the 3D classification iteration 1 results.
		- Run the perl script to reformat the .star file for R; (e.g ./reformat.pl<./run_it001_data.star>./run_it001_data_corrected00.star)
		- Type 'R' to start R (must have R installed on computer).
		- Type the following line by line;
			x= read.table ("run_it001_data_corrected00.star", header= T)
			library(dplyr)
			getmode <- function(v){
			uniqv <- unique(v)
			uniqv [which.max (table(match( v,uniqv))) ]
			}
			modeDataset <- x %>% group_by(GroupNumber, HelicalTubeID)%>% summarise(ClassNumber2=getmode(ClassNumber))
			y<- left_join(x,modeDataset, by=c("GroupNumber","HelicalTubeID"))
			write.table(y,file="run_it001_data_corrected01.star", quote=F,row=F,sep="\t")
			
		- Type 'q()' to quit R.
		- Output file is called run_it001_data_corrected01.star
		
		- Run the script extract_unified_seam_classes_star_relionv3.csh on the .star file created above (e.g source extract_unified_seam_classes_star_relionv3.csh run_it000_data.star run_it001_data.star run_it001_data_corrected01.star 
		  where run_it000_data.star is the star file input into 3D class, run_it001_data.star is the output from the 3D class and run_it001_data_corrected01.star is the file created by R above, with classes corrected to be the most common for each MT). The resulting output file should be called
		  run_it001_data_unified_class.star.
			
21. Correct the seam location for each MT according to 3D class assignments.

		- If you have a 13pf MT run the script 13pf_correct_unified_seam_classes_star_relionv3.csh on the star file created above (e.g source 13pf_correct_unified_seam_classes_star_relionv3.csh run_it001_data_unified_class.star 13 5.56  where
		run_it001_data_unify_seam_class.star is the star file to be converted, 13 is the pf number and 5.56 is the binned pixel size). N.B IMPORTANT THAT THE MT HAS THE RIGHT POLARITY. IF YOU ARE USING A REFERENCE WITH OPPOSITE POLARITY TO THE EXAMPLE CAMSAP REFERENES USE THE SCRIPT 13pf_correct_unified_seam_classes_star_relionv3_REVERSED.csh INSTEAD.
		
		- If you have a 14pf MT run the script 14pf_correct_unified_seam_classes_star_relionv3.csh on the star file created above (e.g source 14pf_correct_unified_seam_classes_star_relionv3.csh run_it001_data_unified_class.star 14 5.56  where
		run_it001_data_unify_seam_class.star is the star file to be converted, 14 is the pf number and 5.56 is the binned pixel size). N.B IMPORTANT THAT THE MT HAS THE RIGHT POLARITY. IF YOU ARE USING A REFERENCE WITH OPPOSITE POLARITY TO THE EXAMPLE CAMSAP REFERENES USE THE SCRIPT 14pf_correct_unified_seam_classes_star_relionv3_REVERSED.csh INSTEAD.
		


------------------------------------------------------------------------------------------------------------------------

HIGH RESOLUTION RECONSTRUCTION

------------------------------------------------------------------------------------------------------------------------	
	  
22.  Rextract centred seam corrected unbinned particles for final reconstruction in helical mode;

		Input Coordinates: N/A (leave empty)
		OR re-extract refined particles? YES
		Refined particles star file: *your refined/seam corrected particle star file created in the last step*
		Reset the refined offsets to 0? No
		Re-center refined coordinates: Yes
		Recenter on- X,Y,Z (pix): 0,0,0
		Particle box size (pix) = ~600A (432 pixels if 1.39A/pixel)
		Rescale particles? No
		Re-scaled size (pixels): N/A
		
		Extract helical segments = yes
		Tube Diameter (A) = 400
		Use bimodal angular priors? yes
		Coordinates are start end only? No
		Cut helical tubes into segments? Yes
		Number of asymmetrical units: 1
		Helical rise (A): 82


23. Scale helical track length; RELION at present does not scale helical track length column properly when changing particle binning.  Therefore you need to run the
script scale_helical_track_length_2.csh (MAKE SURE YOU RUN scale_helical_track_length_2.csh AND NOT scale_helical_track_length.csh) on the new binx1 extracted particles.star file to scale the track length (source scale_helical_track_length_2.csh particles.star 4).
     
	
24.  5th Refine3D, FINAL fine local refinement. Can either run as a run with or without symmetry. With symmetry applied only the PF opposite the seam will be good;
	
		I/O tab:
		
		Input images STAR file: *binx1 output .star file with scaled helical track length created above* (NOT segment averages)
		Reference map: *unbinned decorated synthetic reference of corresponding PF number OR rescaled output reconstruction from previous round* (see examples in
		'references/high_resolution_refinement/' directory)
		Reference mask: 1xbinned soft edged mask encapsulating all non-noise density created with MaskCreate to cover only the central 30% of the MT.
		
		Reference tab:
		
		Reference map is on absolute greyscale? No
		Initial low pass filter (A): 0 (current references are already 15A low pass filtered- if your references are not already low pass filtered this can be set to ~12-15A)
		Symmetry: C1
		
		CTF tab:
		
		Do CTF-correction? Yes
		Has reference been CTF-corrected? Yes
		Have the data been phase flipped: No
		Ignore CTFs until first peak: No
		
		Optimisation Tab:
		
		Mask diameter (A): 460 (REDUCED MASk SIZE TO FOCUS FINE ALIGNMENT ON CENTRAL REGION)
		Mask individual particles with zeros? Yes
		Use solvent-flatened FSCs? Yes
		
		Sampling tab:
		
		Angular sampling interval: 0.9 deg
		Offset search range (pix): ~20A (less than monomer repeat distance to stop boxes jumping back into each other).
		Offset search step (pix) 0.5
		Local searches from auto-sampling: 0.9 deg (restricts to local alignment)
		
		Helix tab:
		
		Do helical reconstruction? Yes
		Tube diameter - inner, outer (A): 120 400 (NOW AN INNER MASK AS YOU WILL BE GENERATING REAL DATA 3D REFS AT EACH ITERATION).
		Angular search range - tilt, psi (deg): 2 2 (smaller search range as we are close enough from last refinement round).
		Keep tilt-prior fixed: No
		Range factor of local averaging: -1 (think this works slightly better at this stage as MTs are rather flexible)
		Apply helical symmetry: No (CHANGE TO YES FOR SYMMETRIC RUN)
		Number of asymmetrical units: 13 (i.e for 13pf kinesin-decorated reconstructions, 12 AUs for proteins that don't bind at the seam like DCX and CAMSAP)
		Intial twist (deg). rise (A): -27.6689 9.46346 (*example is 13pf MT roughly, modify accordingly*)
		Central Z length (%): 30
		Do local searches of symmetry: No (CHANGE TO YES FOR SYMMETRIC RUN)
		Twist search - Min,Max,Step (deg): -27 -28 0.1 (*example is 13pf MT roughly, modify accordingly*)
	   	Rise search - Min,Max,Step (A): 9.4 9.7 0.1 (*example is 13pf MT roughly, modify accordingly*)
		
		Additional arguments: --dont_check_norm --sigma_rot 2 --ignore_helical_symmetry true (REMOVE --IGNORE HELICAL SYM ARGUMENT IF DOING SYMMETRIC RUN, IMPORTANT
		THAT SIGMA ROT IS 2 or 3, RESTRAINS PHI ANGLES TO REMAIN IN INTEGER MINIMUM, also I think the option of 1 iteration does not work in auto-refine, so you will have to stop the job manually for the mo after iteration 1).
	
		N.B Best to run this on GPUs, but as not a massive search range you could try multiple MPIs also (e.g 20 MPIs 1 Thread).
		
------------------------------------------------------------------------------------------------------------------------

DATA OPTIMISATION

------------------------------------------------------------------------------------------------------------------------

The following can be used to try and improve your data once you have your rough alignment parameters and you have moved to 1xbinned data, either before or between iterations of the high-resolution refinement step.
The use of these steps depends hugely on your sample/data.

1. You can try running an additional supervised 3D classification of 2 references: one with decorating protein and one without (just tubulin). 
In my hands this gives you a class with good decoration and one with poor decoration. Do a further refinement on the class with good decoration. This will probably give you a better decorated map and perhaps a better seam but perhaps at the cost of resolution if you don't have many particles.

2. Try iterative rounds of Bayesian Polishing,CTF Refinement and local refinement- can improve reconstruction resolution.

3. Try running a 3D classification with no alignment with 2 classes with your output reconstruction as a reference you will get poorly aligned particles coming out in one of the 2 classes (class will look poor)- then you can keep the good particles.

4. You can manually check particles by plotting ROT angles against particle numbers and checking you have straight lines within particle from each microtubule. This is a good way to check the protocol has worked well for your dataset. To make the plots, you could either use a standard data analysis program like Excel or Prism, or use the mirpy.py script (python mirpy.py -p *_data.star).

The following should be used at earlier stages of the process;

5. You may want to remove short microtubules from your dataset if you are unconvinced by the angular and translational assignments, as a) shorter microtubules give less statistical certainty for pf number/seam allocations and b) particles towards the ends of MTs include less averaged particles in their segment averages and thus have lower signal to noise. The mirpy.py script can be used for this (e.g for minimum length 10 particles; python mirpy.py -s *_data.star -lw 10)

6. You can pick individual MTs based on how internally consistent (particles within a single MT) their protofilament number sorting or seam check class allocations were. This will probably give you a better seam but sometimes at the cost of resolution if you are low on particles. You can use the mirpy.py script can be used for this, to show confidence plots and/or provide a cutoff (e.g to make a confidence plot; python mirpy.py -c run_it001_data.star, to cutoff MTs below 50% confidence of a class assignment; python mirpy.py -c run_it001_data.star -lw 50).
