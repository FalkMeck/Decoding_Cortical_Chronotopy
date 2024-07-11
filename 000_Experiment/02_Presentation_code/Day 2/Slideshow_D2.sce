scenario = "Slideshow_D2.sce";
pcl_file = "Slideshow_D2.pcl";

default_background_color = 40,40,40;

default_font_size = 8;
default_font = "Microsoft JhengHei Light";	
default_text_color = 230, 230, 230;
default_optimize = false;

response_matching = simple_matching;	
response_logging = log_all;			

active_buttons = 1;      
button_codes = 1;

#Fontgroessen
$font_fixcross = 135;
$font_delimiter = 68;
$font_text = 60;

# Koordinaten
$x_coord_zentral = 0;
$y_coord_zentral = 0; 

# Header End

#scenario_type = fMRI;  # für Scanner ##############
scenario_type = fMRI_emulation;		################ # for testing, not connected to fMRI
#scenario_type = fMRI; 						##############
scan_period = 1000;						##############
pulses_per_scan = 1; 					##############
pulse_code = 55; 							

begin; 

## Zeichen
text {caption = "+"; system_memory = true; font_size = $font_fixcross;} fixcross_text;

text {caption = " Willkommen zum letzten Teil des Experiments.

Bleiben Sie bitte einfach entspannt liegen, 
bewegen Sie sich nicht und 
betrachten Sie das Fixationskreuz."; font_size = 60;} Start_text; # change that it starts automatically?

text {caption = "Super!\n Sie haben es geschafft!"; font_size = $font_text;} End_text;

# define instruction, fixation, picture and video
bitmap {filename = "2009-07-Kollosseum.JPG"; scale_factor = 1.25;} my_bitmap;

picture {bitmap my_bitmap; x = 0 ; y = 0;} pic;

picture {} default;

picture {text Start_text; x = 0; y = 0;} Start_pic;

picture { text fixcross_text; x = $x_coord_zentral; y = $y_coord_zentral;} fixcross_pic;

picture {text End_text; x = 0; y = 0;} End_pic;


# Slideshow trial
trial {
	trial_type = fixed;
	trial_duration = 4000;
	stimulus_event {	
		picture pic;
	} pic_event;
} pic_trial;


# empty trial
trial{
	trial_type = fixed;
	trial_duration = 2000;
	picture default;}mt_trial; 

# Welcome and Ready trial
trial {
	trial_duration = forever;
	trial_type = specific_response;
		terminator_button = 1;
	picture Start_pic;
} Start_trial;

# Fixcross Trial
trial {
   trial_type = specific_response;
	trial_duration = forever;
	terminator_button = 1;
	picture fixcross_pic;
} fixcross_trial;

#EndeTrial
trial { 
			trial_type= specific_response;
			trial_duration = forever;
			terminator_button = 1;
			picture End_pic;} End_trial;
			
			

###################################################
###				   fMRT Standard			   	   ### #????		################
###################################################

# for write_log subroutine	#das ist ein Standard trial - IMMER SO UEBERNEHMEN BEI MR-MESSUNGEN
#trial {	
#stimulus_event {
#nothing {};
#time = 0;
#code = "nothing";	#diesen code machen wir uns nacher zu nutze
#}evt_log;
#} param_log;

