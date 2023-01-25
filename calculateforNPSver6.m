function [pulses]=calculateforNPSver6(indices,ym,yasls,samplerate,N,Length,Deff,szlength,sqlength,hchannel,pulses,sqchan)
%indices should be matrix with 2+numsq*2 columns: 1=start sz 2=endsz 3=start sq 4=
%end sq... etc. 

%%pulses columns: 1) filename, 2) sizing start, 3) size end, 4) Iavg sz, 5)
%Iavgbaseline, 6)diameter, 7)uflow (sizing), 8) sqstart1, 9) sqend1, 10-13)
%same as 8,9 but for sq2,3, 14-16) vratio1-3, 17-19) strain1-3, 20) tau,
%21) Rsquared, 22) numpast sizing, 23) pressure, 24) geometry, 25) device
%num, 26) voltage, 27) dinput 


%pulses columns: 1) filename, 2) sizing start, 3) size end, 4) Iavg sz, 5)
%Iavgbaseline, 6)diameter, 7)uflow (sizing), 8) sqstart1, 9) sqend1, 10)wcdi 11)strain 
%12) pressure, 13) voltage, 14) dinput, 15) noise level used in asls2 

[m,~]=size(indices);
pulsesforprint=num2cell(pulses);
sizeptsexist=0;

      
for i=1:m 
   
    if indices(i,1) ~=0 && indices(i,2) ~=0 %are there sizing points recorded (this shoudl always be true) 
        sizeptsexist=1;
        szstart=indices(i,1);
        szend=indices(i,2);
        
        pulses(i,2)=szstart;
        pulses(i,3)=szend;
        pulses(i,4)=mean(ym(szstart:szend,1));%Ipore
        pulses(i,5)=mean(yasls(indices(i,1):indices(i,2),1));%Ibaseline
        %calculate cell diameter
        deltaI=pulses(i,5)-pulses(i,4);
        deltaIoverI=deltaI/pulses(i,5);%deltaI/I
        dtop=deltaIoverI*Length*(Deff^2);
        dbottom=1+(0.8*deltaIoverI*Length/Deff);
        pulses(i,6)=(dtop/dbottom)^(1/3); %calculated cell diameter
 
        deltaT=(szend-szstart)*N/samplerate;
        pulses(i,7)=szlength/deltaT; %flowspeed
        pulses(i,11)=(pulses(i,6)-sqchan)/pulses(i,6); %strain 
    end
    
    %calculate wCDI:
    if indices(i,3) ~=0 && indices(i,4)~=0 % is sq 1 is recorded
        sqstart=indices(i,3);
        sqend=indices(i,4);

        pulses(i,8)=sqstart;
        pulses(i,9)=sqend;
        deltaTsq=(sqend-sqstart)*N/samplerate; %time in squeeze channel (sec) 
        vsq=sqlength/deltaTsq; %velocity of particle in squeeze channel (microns/sec)
        if sizeptsexist==1
            pulses(i,10)=(vsq/pulses(i,7))*(pulses(i,6)/hchannel); %wCDI
        end
    end
    sizeptsexist=0;

end
