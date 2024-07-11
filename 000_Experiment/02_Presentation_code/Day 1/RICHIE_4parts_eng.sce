# Header
scenario = "RICHIE_4parts_eng.sce";
pcl_file = "RICHIE_4parts_eng.pcl";

default_font_size = 8;
default_font = "Microsoft JhengHei Light";	
default_background_color = 40,40,40;###oder:150,150,150
default_text_color = 230, 230, 230;###oder:255,255,255

screen_width=1920;
screen_height=1080;
screen_bit_depth = 32;

response_matching = simple_matching;	################
response_logging = log_all;				################

active_buttons = 2;      
button_codes = 1,2;

#Fontgroessen
$font_fixcross = 135;
$font_delimiter = 68;
$font_text = 60;
$font_numbers = 180;

# Koordinaten
$x_coord_zentral = 0;
$y_coord_zentral = 0; 

# Header End

scenario_type = fMRI;  # für Scanner
#scenario_type = fMRI_emulation;  # für Scanner
scan_period = 1000;						################
pulses_per_scan = 1; 					################
pulse_code = 55; 							################




begin;

##########################################################################
#		
#		Text
#
##########################################################################

## Zeichen
text {caption = "+"; system_memory = true; font_size = $font_fixcross; } fixcross_text;
text {caption = " "; system_memory = true; font_size = $font_delimiter; } delimiter_text;

text {caption = "Welcome to this test!"; font_size = $font_text;} Start_text_test;

text {caption = "Welcome!

In the following experiment you will see blocks
of number sequences. These sequences are structured 
in regular segments of vayring length depending on the block.

Independent of the structure your tasks is to press the 
answer button as soon as possible 
when you see the number 0. 

When you are ready,
please press the answer button."; font_size = 45;} Start_text; # change that it starts automatically?

text {caption = "Nice!\n You are done!

Please remain calm for a little while longer."; font_size = $font_text;} End_text;

text {caption = "Next block:\n\n"; font_size = $font_text;} block_break_text;

text {caption = "Block_name"; font_size = $font_text;}next_block_text;

text {caption = "No."; font_size = $font_numbers; font = "Microsoft JhengHei Light"; max_text_height = 1080; max_text_width = 1920;
#font = "Ink Free";
} trial_number; 



##########################################################################
#		
#		Pictures
#
##########################################################################

picture {} default;

# define instruction, fixation, picture and video
picture { text Start_text; x = 0; y = 0; } Start_pic;
picture { text fixcross_text; x = $x_coord_zentral; y = $y_coord_zentral; } fixcross_pic;

picture { text block_break_text; x = 30; y = 0; text next_block_text; x = -30; y = 0;} block_break_pic;

picture { text delimiter_text; x = $x_coord_zentral; y = $y_coord_zentral;} empty_pic;

picture { text trial_number; x = $x_coord_zentral; y = $y_coord_zentral;} trial_number_pic;

picture { text End_text; x=0; y=0;} End_pic;


##########################################################################
#		
#		Trials (also includes some unused ones)
#
##########################################################################

# Welcome and Ready trial

trial {

	trial_duration = forever;
	trial_type = specific_response;
		terminator_button = 1;
	picture Start_pic;
} Start_trial;


#Number Trial
trial {
   trial_type = fixed;
	trial_duration = 700;
	#trial_duration = 1500;
	   picture default;
	   stimulus_event {
		picture trial_number_pic;
		time = 0;
		duration = 400;
		#duration = 1500;
		code = "Num";
		target_button = 1;
   } number_event;
  }number_trial; 

#Fix Blink trial in same timing as Numbers
trial {
   trial_type = fixed;
	trial_duration = 700;
	   picture default;
	   stimulus_event {
		picture trial_number_pic;
		time = 0;
		duration = 400;
		code = "Fix_blink";
		#target_button = 1;
   } fix_blink_event;
  }fix_blink_trial; 


# Fixcross Trial
trial {
   trial_type = fixed;
	trial_duration = 2000;
		 stimulus_event {	
			picture fixcross_pic;
			code="Fix";
			time=0;  
			duration = 1700; 
		} fix_event;
} fixcross_trial;

# Break
trial{trial_type = fixed;
		trial_duration = 3000;
		picture block_break_pic;}break_trial;
		
#EndeTrial
trial { 
			trial_type= specific_response;
			trial_duration = forever;
			terminator_button = 2;
			picture End_pic;} End_trial;
			
			
# progress bar
# borders
box { height = 2; width = 1500; color = 230,230,230; } horiz;
box { height = 100; width = 2; color = 230,230,230; } vert;
box { height = 100; width = 1500; color = 230,230,230;} bar; 


trial {
	trial_duration = 100;
	picture{
		text{caption = "caption"; font_size = $font_text;}short_break; x = 0; y = 0;
      box horiz; x = 0; y = -350;
      box horiz; x = 0; y = -450;
      box vert; x = -750; y = -400;
      box vert; x = 750; y = -400;
		box bar; x = 0; y = -400;}bar_pic;
}progress_bar;

trial {trial_type = specific_response; terminator_button = 2; #wird extern vom Kontrollraum weiter gedrückt mit E
		trial_duration = forever;
		picture{text{caption = "caption"; font_size = $font_text;}short_break_cont; x = 0; y = 0;}short_break_pic;
		code = "break";
} br_trial;



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

