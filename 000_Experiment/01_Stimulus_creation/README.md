MATALB Code
create_all_block_orders
- Creates n (n = 50) different orders of blocks for al 4 possible starting blocks
- while checking that all the 50 order are indeed differnt
- block order are based on DeBruijn sequences, transitions between blocks are equally distributed 
 Change to normal DeBruijn sequences: the same block cannot repeat
- requires the functions create_block_order
create_block_order3
- Creates a sequence of two DeBruijn sequences, appends them and removes same numbers back to back, so the same blocks cannot repeat
- requires the input of a number (1-4) and generates a DeBruijn seuqenece with teh respective startint number
- requires the function debruijn_generator (W. Owen Brimijoin (2024). de Bruijn sequence generator (https://www.mathworks.com/matlabcentral/fileexchange/28015-de-bruijn-sequence-generator), MATLAB Central File Exchange. Abgerufen 9. Juli 2024.)
- as sequences are circular with the starting number, all blocks have the same frequency 
- also last number ~= first number 
- also first 5 Blocks are A B C D A, so first 4 blocks can be used fro trainign al 4 different blocs

DeBruijn_generator_test
- creates DeBruijn sequecnes Level1/Single 
- 4 sequences à 36 = 144 -> 6* to fill entire experiment
- checks, if 1-3 and 2-4 dont have the same sequences of the length 3 - 9 und that 1 und 4 not exactly the same
- tries 1.000.000 times
- does not necessarily find a solution for equivalents of the length 3 - 9 and also no circularity
- save sequence as mat and text file 
- theoretically a loop for all 6 order, but its too slow as a loop
- renames number for easier translation to other code

DeBruijn_Level2_test
- same as DeBruijn_generator_test but for Level2/triplet
- creates one DeBruijne sequence for 1,2,3,4  is repalces with 123,321,567,765  and later with 456,654,789,987
- controlls for circularity and triplet to nonet repetitions
- code creates 3 sequences (16*3*3 = 144) and saves them

DeBruijn_Level3_test
- same as DeBruijn_generator_test and DeBruijn_Level2_test but for Level3/Nonet
- less controll necessary, only one sequnce per 144 digits, only cirkularity test
-  renames Benennt die 1,2,3,4 dann um in 123123123, 321321321, 567567567, 765765765, and than later in Presentation calls them  zu 456456456, 654654654, 789789789, 987987987 

Define_violations
- loads input lists for the experiment (rpl und ins) and defines for every inserted violation if it is single only, triplet only, nonet only or concerns a compelte transition
- 1 = single, 2 = triplet, 3 = nonet, 4 = full  always ¼ of kind of transistion are violations, in total 36
- Single/Level 1: 36/144 single only
- Triplet/Level 2: 24/96 single only, 12/48 triplet only
- Nonet/Level 3: 24/96 single only, 8/32 triplet only, 4/16 nonet only
- Complete/Level 4: 24/96 single only, 8/32 triplet only, 3/12 nonet only, 1/4 full
- Creates text files that mirror input files, will be read in parallel in Presentation

RRP_generator_v3
- generates input lists
- defines perfect Transition Matrix for all 4 Levels 
- generates transistion positions for all Kinds of transistions, draws violated transisiton 1-4 times 
	- Level 1: every Transition 2 times, violated once
	- Level 2: single only are 8 different transistion, everyone is violated 1-2 times, 16 differnt triplet 		transsition 6 violations 
	- Level 3: 1-2 time 8 single only, 1 times most common triplet transistions, 2 of the 12 nonet only 	transsitions
	- Level 4: similar to Level 3, 4 full transsition, once violated 
- calcultes transsition Matrix from violated sequence
- compares with perfect transistion Matrix and sums distances: thresholds: Level 1 = 20, Level 2+3 = 10, Level 4 = 0
	- a bit arbitrary but otherwise it does not converge	
- 18 transistion per block -
- uniform distribution of violations across all lists
Txt_Exprt_rpl_and_ins
- genereates Input files 

# -------------------------------------------------------------------------------------------------------------#

R Code

find_isi.R
- creates one set for interstimulus intervals with added noise to better estimate the BOLD function
- are matched to the same blocks to have same distance fro TR for all subject to have synchrony for inter-subejct correlation analysis

WuS2013Similarity.R (for the necessary files please contact f_meck01@uni-muenster.de)
- finds the font with the lowest Pixel overlap for numbers in Windows
- choose the best one that is still readable 
Wong, B., & Szücs, D. (2013). Single-digit Arabic numbers do not automatically activate magnitude representations in adults or in children: Evidence from the symbolic same–different task. Acta Psychologica, 144(3), 488–498. https://doi.org/10.1016/j.actpsy.2013.08.006

- 
- 
