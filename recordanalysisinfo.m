function [pulses]=recordanalysisinfo(pulses,dinput,info,noise)
%the purpose of this function is to populate the excel table with the
%conditions of this experimental run, most of these conditions are stored
%in the filename, but will be input by the user, and should be passed along
%if this is pt>1 


pulses(:,12)=info(1,1);%pressure
pulses(:,13)=info(1,2);%voltage 
pulses(:,14)=dinput; 
pulses(:,15)=noise;
