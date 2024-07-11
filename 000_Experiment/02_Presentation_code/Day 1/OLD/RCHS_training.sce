# Header
scenario = "RCHS_training.sce";
pcl_file = "RCHS_training.pcl";

default_font_size = 8;
default_font = "Microsoft JhengHei Light";	
default_background_color = 155,155,155;
default_text_color = 255, 255, 255;

active_buttons = 1;      
button_codes = 1;


#Fontgroessen
$font_fixcross = 90;
$font_delimiter = 45;
$font_text = 30;

# Koordinaten
$x_coord_zentral = 0;
$y_coord_zentral = 0; 

# Header End

begin;

##########################################################################
#		
#		Text
#
##########################################################################

## Zeichen
text {caption = "+"; system_memory = true; font_size = $font_fixcross; } fixcross_text;
text {caption = " "; system_memory = true; font_size = $font_delimiter; } delimiter_text;

text {caption = "Willkommen!

Im folgenden Training werden Ihnen vier Blöcke von 
Zahlenfolgen gezeigt. 

Ihre Aufgabe ist es dabei so schnell wie möglich 
auf die LEERTASTE zu drücken, 
wenn die Zahl 0 präsentiert wird.

Wenn Sie bereit sind,
drücken Sie bitte die LEERSTASTE"; font_size = 30;} Start_text_ohne;

text {caption = "Willkommen!

Im folgenden Experiment werden Ihnen mehrere Blöcke von 
Zahlenfolgen gezeigt. Diese Zahlenfolgen sind je nach Block 
in regelmäßige Abschnitte unterschiedlicher Länge strukturiert.

Ihre Aufgabe ist es dabei unabhängig von der Strukturierung, 
so schnell wie möglich auf die Antwort-Taste zu drücken, 
wenn die Zahl 0 präsentiert wird.

Wenn Sie bereit sind,
drücken Sie bitte auf die Antwort-Taste"; font_size = 30;} Start_text; # change that it starts automatically?

text {caption = "Super!\n Sie haben es geschafft!"; font_size = $font_text;} End_text;

text {caption = "Nächster Block:\n\n"; font_size = $font_text;} block_break_text;

text {caption = "Block_name"; font_size = $font_text;}next_block_text;

text {caption = "No."; font_size = 90; font_color = 255, 255, 255; font = "Microsoft JhengHei Light";
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
picture { text Start_text_ohne; x = 0; y = 0; } Start_pic_ohne;

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

trial {

	trial_duration = forever;
	trial_type = specific_response;
		terminator_button = 1;
	picture Start_pic_ohne;
} Start_trial_ohne;

#Number Trial
trial {
   trial_type = fixed;
	trial_duration = 700;
	   picture default;
	   stimulus_event {
		picture trial_number_pic;
		time = 0;
		duration = 400;
		code = "Num";
		target_button = 1;
   } number_event;
  }number_trial; 


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
trial{trial_type = specific_response;
		trial_duration = forever;
		terminator_button = 1;
		picture block_break_pic;}break_trial;
		
#EndeTrial
trial { 
			trial_type= specific_response;
			trial_duration = forever;
			terminator_button = 1;
			picture End_pic;} End_trial;


