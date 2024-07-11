scenario = "Slideshow_D1.sce";
pcl_file = "Slideshow_D1.pcl";

default_background_color = 40,40,40;
default_text_color = 230, 230, 230;
default_optimize = false;

response_matching = simple_matching;	
response_logging = log_all;			

active_buttons = 1;      
button_codes = 1;

begin; 

# define instruction, fixation, picture and video
bitmap {filename = "2009-07-Kollosseum.JPG"; scale_factor = 1.25;} my_bitmap;

picture {bitmap my_bitmap; x = 0 ; y = 0;} pic;

# Slideshow trial
trial {
	trial_type = fixed;
	trial_duration = 4000;
	stimulus_event {	
		picture pic;
	} pic_event;
} pic_trial;


