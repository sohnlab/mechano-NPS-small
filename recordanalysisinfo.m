function [pulses]=recordanalysisinfo(pulses,dinput,info,noise)
%the purpose of this function is to populate the excel table with the
%conditions of this experimental run, most of these conditions are stored
%in the filename, but will be input by the user, and should be passed along
%if this is pt>1 

pulses(:,27)=dinput; 
pulses(:,23)=info(1,1);
pulses(:,24)=info(1,2);
pulses(:,25)=info(1,3);
pulses(:,26)=info(1,4);
pulses(:,28)=noise;
