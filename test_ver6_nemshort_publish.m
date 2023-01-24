clear
clc 
close all 

%choose file to analyze
[name,path] = uigetfile;
filename=fullfile(path,name);
load(filename);

%restart from existing incomplete analysis 
if exist('newdata','var')
    data=newdata;
else
data= data';
end

%set a dinput= the minimum diameter particle you wish to detect/measure,
%the threshold for detection will be slightly below this diameter to ensure
%all particles with diameter=dinput are detected
dinput=2.5;

%% read all, search for pulses
    
    [ym, yasls,cornercontext,cornerindex,ydetrend,pulses,pulsesforprint] = mNPS_ver6_nemshort_publish(data, sampleRate,name,dinput);
  
        
