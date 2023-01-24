 function [ ym, yasls,cornercontext,cornerindex,ydetrend,pulses,pulsesforprint ] = mNPS_ver6_nemshort_publish( data, sampleRate,name,dinput)

saveto=strcat('analyzed_',name); %enables saving as you analyze



%% parameters specific to device 
sqchan=2.9; %microns, width of squeeze channel
sizingchan=5.5; %microns, width of sizing channel 
Length=1550; %microns, overall length of device (including entrance node and exit node) 
szlength = 700; %microns, length of sizing channel 
sqlength = 700; %microns, length of the squeeze channel 
hchannel = 9.9; %microns 
Deffective=10.19; %determined by measuring beads of known size in your device 

% calculate average deltaI/I for different values of interest
X=(dinput^3)/((Deffective^2)*Length)*(1/(1-0.8*((dinput/Deffective)^3))); %for dinput= minimum diameter we are interested in measuring 
Q=((sqchan)^3)/((Deffective^2)*Length)*(1/(1-0.8*((sqchan/Deffective)^3)));%for when diameter=sqchannel width--> threshold for 0% strain 
Y=(5^3)/((Deffective^2)*Length)*(1/(1-0.8*((5/Deffective)^3))); %CONSIDER DELETING
T=(sizingchan^3)/((Deffective^2)*Length)*(1/(1-0.8*((sizingchan/Deffective)^3))); %for when diameter=sizing channel width --> threshold above which we can't accurately measure diameter

% parameters specific to device and pressure
%   - entire pulse length on avg= 63 (73), with highest around 72 (77)
%   - avg pore length = 25 (31)
%   - avg length between pre-pore-drop and start of squeeze = (49)
%   - avg length between pre-pore-drop and end of squeeze = (67)
period = 123; % therefore encompasses the entire pulse length + 2 before start index
S1=49;
S2=55;
S3=67;
S4=linspace(30,40,11)'; %for use with 11 microns/ fast pulses
pulseperiod = 9; % average pulse length is 25, therefore it allows you
% to average over a significant portion- smallest pulse reported in
% original batch was 18 
wp=698; %what you have to add to encompass whole post 
interval=40; %less than this and there is no pulse (subject to change) (was 30 but then the p100 in 10L) (was 20- then 196 in 10)
numrecov=0; %number of recovery pulses
numpul=2; %number of total pulses in a cell event 
numsq=1;

%% option to continue where you left off last time 
%do we want to continue analyzing from a file we left off on? 
contanalysis=input('do you want to continue analysis on this dataset where you left off last time? 1=yes');
if contanalysis==1
    numanalysis=input('what part of analysis are we on now?');
    load(saveto,'stderror');
    if numanalysis == 2 
        load(saveto,'pulses')
        saveto=saveto(1:end-4);
        saveto=strcat(saveto,'pt',num2str(numanalysis));
    else
        saveto=saveto(1:end-4);
        saveto=strcat(saveto,'pt',num2str(numanalysis-1));
        load(saveto,'pulses');
        saveto=saveto(1:end-1);
        saveto=strcat(saveto,num2str(numanalysis));
    end
    info=pulses(1,23:26);
    %cutoff data after the last pulse analyzed 
    [numpulse,~]=size(pulses);
    lastpulse=pulses(numpulse,2);
    clear pulses
else
    %ask user for info 
    info=zeros(1,4);
    userinfop=input(strcat(name,'what is the pressure?'));
    if isfloat(userinfop)==0
        userinfop=input(strcat(name,'what is the pressure?'));
    end
    info(1,1)=userinfop;
    userinfod=input(strcat(name,'what is the device number?'));
    if isfloat(userinfod)==0
        userinfod=input(strcat(name,'what is the device number?'));
    end
    info(1,3)=userinfod;
    userinfod=round(userinfod);
    switch userinfod
        case {1,4,5}
            info(1,2)=789;
        case {6,10}
            info(1,2)=978;
        case {7,8,11}
            info(1,2)=897;
        case {3,9}
            info(1,2)=14;
        case{2}
            info(1,2)=987;
    end
    
    userinfov=input(strcat(name,'what is the voltage?'));
    if isfloat(userinfov)==0
        userinfov=input(strcat(name,'what is the voltage?'));
    end
    info(1,4)=userinfov;
end

%% data pre-processing
N = 20; %for downsampling
smoothingparam= 100; %for smoothing 
ysmoothed = fastsmooth(data,smoothingparam,1,1);  % perform rectangular smoothing
ym = downsample(ysmoothed,N); %downsample by period N


%calculate baseline (yasls), and 2 more filtered signals (yasls2, yasls3)
lambda = 1e8; %larger means smoother background 
p=0; % less than 0.5 means negative is strongly supressed %0.02
maxiter=20; % maximum iterations 
noise=1.6e-4; %found in first data run analyzed (specific to experimental setup) 
aslsparam1=struct('lambda', lambda, 'p', p, 'max_iter', maxiter, 'noise_margin', noise);
aslsparam2=struct('lambda', 5, 'p', 0.5, 'max_iter', maxiter, 'noise_margin', noise);
aslsparam3=struct('lambda',1e4, 'p', 0.5, 'max_iter', maxiter, 'noise_margin', noise);

yasls = -1*ASLS2(-1*ym,aslsparam1); %baseline
yasls2=-1*ASLS2(-1*ym,aslsparam2); %smoothed signal
yasls3=-1*ASLS2(-1*ym,aslsparam3); %smoothest signal 

ydetrend = ym-yasls; %subtract the baseline from the signal, bringing it flat and to zero
expectedpore=yasls.*(1-X); %this is the threshold for the actual pulse
diameter5 = yasls.*(1-Y); %CONSIDER DELETING 
upperlimit= yasls-((yasls-expectedpore).*0.9); % this is the initial threshold


sqlimit= yasls.*(1-Q);


ydiff=diff(ydetrend);
ydiff2=[0;ydiff];
ydiff2=(ydiff2.^3)*(10^10); %larger differences are enhanced exponentially
yasdet=yasls2-yasls;
yasdiff=diff(yasls2-yasls);%this was yasls3 until there was noise and it made it bad 
yasdiff=[0;yasdiff];
yasdiff=(yasdiff.^3)*(10^10);
save(saveto,'ym','yasls','ydetrend','ydiff2','data','yasls2');

%% plotting post processing 
    xtitle='Sample Number';
    ytitle='Current';
    % for plotting post downsampling 
    lengthpost = length(ym);
    
    f1=figure('Name','Ym with Baseline'); %figure 1 
    a1=axes;
    xaxispost = linspace(1,lengthpost,lengthpost)';
    plot(xaxispost,ym)
    hold on 
    plot(xaxispost,yasls)
    
    %for plotting expected pore pulse drop  
    plot(xaxispost,expectedpore)
    plot(xaxispost,upperlimit)
    plot(xaxispost,diameter5) %line denoting cells with a diameter of 5 microns, CONSIDER DELETING 
    plot(xaxispost,yasls2) 
    plot(xaxispost,yasls3)
    xlabel(xtitle);
    ylabel(ytitle);
    sz=10;
    mkr='.';
    hold(a1,'on');
    legend('ym','yasls','expectedpore','upperlimit', 'diameter5','asls2','asls3') %CONSIDER DELETING DIAMETER 5
    
    f2=figure('Name','Data Detrended'); 
    a2=axes;
    plot(a2,xaxispost,ydetrend,'r');
    hold on
    scatter(a2,xaxispost,ydetrend,sz,'r');
    xlabel(xtitle);
    ylabel(ytitle);
    hold on 
%     scatter(a2,xaxispost,ydiff2.*.001,sz,'b')
%     plot(a2,xaxispost,ydiff2.*.001,'b')
    plot(a2,xaxispost,upperlimit-yasls,'c')
    scatter(a2,xaxispost,yasdiff.*.001,sz,'m')
    plot(a2,xaxispost,yasdiff.*.001,'m')
    plot(a2,xaxispost,yasdet,'y')
    plot(a2,xaxispost,yasls3-yasls,'c')
    
    
    %pick region for stdev calculation --> uncomment when you are doing the
    if exist('stderror','var')==0
        d=datacursormode(f1);
          d.Enable='on';
         d.DisplayStyle='window';
         startstdi=input('click where you want region to start, then press enter');
         if isempty(startstdi)==1
             vals=getCursorInfo(d);
             startstd=vals.Position(1,1);
         end
         d.Enable='off';
         d.Enable='on';
         endstdi=input('click where you want region to end, then press enter');
         if isempty(endstdi)==1
             vals=getCursorInfo(d);
             endstd=vals.Position(1,1);
         end
         stderror=std(ydetrend(startstd:endstd))
              %stderror=        0.00036256116737692 %this is for training purposes
    end
    
    
    %back to plotting
    f3=figure('Name','Difference');
    a3=axes;
          scatter(a3,xaxispost,ydiff2,sz,mkr);
          xlabel(xtitle);
          ylabel(ytitle);
          hold(a3, 'on');

%% classic thresholding
threshpts = zeros(2,2); %list of pts in ym that cross the threshold 
%second column is whether its trending down or up (down=0, up =1, 3=error?) 
% it seems like noise is cross the threshold before the pulse actually
% finishes- do this crossing of threshold with a moving average (yasls3)
k=1;
i=100+wp;
if exist('lastpulse','var')==1
    i=lastpulse;
end
recovstop=0;
[threshpts]=classicthresholding(yasls3,upperlimit,threshpts,i,k,lengthpost,recovstop);
% while i<lengthpost-2
%     if yasls3(i)==upperlimit(i) % if a data point== the threshold
%         threshpts(k,1)=i;
% %         xline(a2,i,'k');
%         %determine whether trending down or up? 
%         if yasls3(i-1)>upperlimit(i) && yasls3(i+1)<upperlimit(i+1) %trending down 
%             threshpts(k,2)=0;
%         elseif yasls3(i-1)<upperlimit(i) && yasls3(i+1)>upperlimit(i+1) %trending up 
%             threshpts(k,2)=1;
%         else
%             threshpts(k,2)=3;
%         end
%             
%         k=k+1;
%     elseif yasls3(i)>upperlimit(i) && yasls3(i+1)<upperlimit(i+1) && yasls3(i+2)<upperlimit(i+2) % if btwn two points going down
%         threshpts(k,1)=i;
%         threshpts(k,2)=0;
% %         xline(a2,i,'k');
%         k=k+1;
%     elseif yasls3(i)<upperlimit(i) && yasls3(i+1)>upperlimit(i+1) && yasls3(i-1)<upperlimit(i-1) % if btwn two points going up 
%         threshpts(k,1)=i+1;
%         threshpts(k,2)=1;
% %         xline(a2,i+1,'k');
%         k=k+1;
%     end
%     i=i+1;
% end

save(saveto,'stderror','-append');
%% remove threshpts that are just error

[rowthresh,~]=size(threshpts);
i=1;
[threshpts,rowthresh]=fixthreshptserror(threshpts,interval,i,rowthresh);

% if threshpts(i,2)==1
% %     xline(a2,threshpts(i),'b');
%     threshpts(i,:)=[];
%     rowthresh=rowthresh-1;
% end
% while i<rowthresh
%     if threshpts(i+1,1)-threshpts(i,1)<interval
%         if threshpts(i,2)==0 && threshpts(i+1,2)==1
% %             xline(a2,threshpts(i),'b');
% %             xline(a2,threshpts(i+1),'b');
%             threshpts(i,:)=[];
%             threshpts(i,:)=[];
%             rowthresh=rowthresh-2;
%         elseif threshpts(i,2)==1     %idk what this section is?? it might be bad start
% %             xline(a2,threshpts(i),'b');
%             threshpts(i,:)=[];
%             rowthresh=rowthresh-1;
%         end
%     else
%         i=i+2;
%     end
%     
% end

for i=1:rowthresh
    xline(a2,threshpts(i),'k');
end




%% look at the windows 
f4= figure('Name','Window of Interest');
a4=subplot(2,1,1);
a5=subplot(2,1,2);
hstep= 15; %horizontal line step size (og 10)
hstepog=hstep;
vstep=3;
vstepr=4; %rising and falling lines step size (og 3)
steps=[hstep hstepog vstep vstepr];
cornercontext=zeros(1,23);
cornerindex=zeros(1,23);
cornerdiff=zeros(1,23);
guesscheck=zeros(1,4); %col1= startguess, col2=startactual, col3=endguess, col4=endactual
guesscheck(:,:)=-1;
i=1;
useranswer=1;
k=1; % the number of corner *sets* logged

%corner matrices, each row is different pulse segment (sizing, cont, recov 1-10), set of 13 rows can be from 1 cell 
%col: 1= pulse type; 2=actual start pulse; 3=actual endpulse; 4:13= the
%points surrounding the starting predicted corner (5 on each side, if it falls
%between two points); 14:23= the points surrounding the ending predicted
%corner 

while i<rowthresh && useranswer==1 && threshpts(i+1,1)+wp<lengthpost
    a=threshpts(i,1)-3;
    b=threshpts(i+1,1)+3;
    cla(a4)
    reset(a4)
    cla(a5)
    reset(a5)
        plot(a4,xaxispost(a:b),ydetrend(a:b),'r');
        plot(a5,xaxispost(a-wp:b+wp),ym(a-wp:b+wp),'r');
        
        hold(a4,'on')
        hold(a5,'on')
        
        plot(a5,xaxispost(a-wp:b+wp),yasls(a-wp:b+wp))
        plot(a5,xaxispost(a-wp:b+wp),upperlimit(a-wp:b+wp),'c')
        plot(a5,xaxispost(a-wp:b+wp),expectedpore(a-wp:b+wp),'k')
        plot(a5,xaxispost(a-wp:b+wp),yasls3(a-wp:b+wp),'b')
        plot(a5,xaxispost(a-wp:b+wp),yasls2(a-wp:b+wp),'y')
        plot(a5,xaxispost(a-wp:b+wp),sqlimit(a-wp:b+wp),'g')
        xlabel(a5,strcat('expected pore, in black, =',num2str(dinput)));
        legend(a5,{'ym','yasls','upperlimit','expectedpore','yasls3','yasls2','sqlimit'});
        
        
        scatter(a4,xaxispost(a:b),ydetrend(a:b),sz,'r');
        xlabel(a4,xtitle);
        ylabel(a4,ytitle); 
        %scatter(a4,xaxispost(a:b),ydiff2(a:b).*.001,sz,'b')
        %plot(a4,xaxispost(a:b),ydiff2(a:b).*.001,'b')
        plot(a4,xaxispost(a:b),upperlimit(a:b)-yasls(a:b),'c')
        scatter(a4,xaxispost(a:b),yasdiff(a:b).*.001,sz,'m')
        plot(a4,xaxispost(a:b),yasdiff(a:b).*.001,'m')
        plot(a4,xaxispost(a:b),yasdet(a:b),'y')
        
        
        xline(a4,threshpts(i,1),'b')
        xline(a5,threshpts(i,1),'b')
        xline(a4,threshpts(i+1,1),'b')
        xline(a5,threshpts(i+1,1),'b')
   
    
   
    
    
    
    useranswer0=input('bad start=0, next=1, analyze=2,skip sq and recov pulses=3,skip entire cell event=4, go back 1=5, go back 2 pulses=6, stop=anythingelse');
    
    if useranswer0==0
        i=i+1;
    elseif useranswer0==2
        % categorize as sizing, squeeze, or recovery for training purposes 
        
        
            useranswer2=input('is it sizing(1),squeeze (2), or recovery (3)');
            if isempty(useranswer2)
                 useranswer2=input('is it sizing(1),squeeze (2), or recovery (3)');
            end
            if useranswer2 ~= 1 && useranswer2 ~=2 && useranswer2~=3
                useranswer2=input('is it sizing(1),squeeze (2), or recovery (3)');
            end
            if useranswer2 ~= 1 && k==1
                 useranswer2=input('is it sizing(1),squeeze (2), or recovery (3)');
            end % CHECK TO MAKE SURE THIS WORKS 
            switch useranswer2
                case 1
                    cornercontext(k,1)=1;
                    cornerindex(k,1)=1;
                    cornerdiff(k,2)=1;
                    hstep=hstep+40;
                case 2
                    cornercontext(k,1)=2;
                    cornerindex(k,1)=2;
                    cornerdiff(k,2)=2;
                    hstep=hstep+10;
                case 3
                    cornercontext(k,1)=3;
                    cornerindex(k,1)=3;
                    cornerdiff(k,2)=3;
            end
            if k<3 || cornercontext(k,1)==1 || cornercontext(k,1)==2
            
                 %make horizontal line
                [~,minasls]=min(yasdet(a+3:b-3));
                minasls=minasls+a-1;
                xline(a4,minasls,'r')
                    %ensure that hstep does not exceed the threshpts
                    if minasls-hstep<threshpts(i,1) || minasls+hstep>threshpts(i+1,1)
                        if threshpts(i+1,1)-minasls > minasls-threshpts(i,1)
                            hstep=minasls-threshpts(i,1);
                        else 
                            hstep=threshpts(i+1,1)-minasls;
                        end
                    end

                    %
                horizregion=ydetrend(minasls-hstep:minasls+hstep);

                [minregion,~]=min(horizregion);
                xlabel(a4,strcat('hstep=',num2str(hstep)));

                %plot horiz line higher with error previously calculated
                horizline=minregion+stderror*2;%is normally *2 but with 60 Hz noise you need to to up it 
                yline(a4,horizline,'r')
                [~,mindiff]=min(yasdiff(a:minasls));
                mindiff=mindiff+a-1;
                [~,maxdiff]=max(yasdiff(minasls:b));
                maxdiff=maxdiff+minasls-1;    %all of it used to be vstepr+1
                p1=polyfit(xaxispost(mindiff-vstepr+2:mindiff+2),ydetrend(mindiff-vstepr+2:mindiff+2),1);
                x1=xaxispost(mindiff-vstepr+2:mindiff+2+hstep);
                y1=polyval(p1,x1);
                plot(a4,x1,y1,'k')
                p2=polyfit(xaxispost(maxdiff-vstepr+2:maxdiff+2),ydetrend(maxdiff-vstepr+2:maxdiff+2),1);
                x2=xaxispost(maxdiff+vstepr-hstep-2:maxdiff+2);
                y2=polyval(p2,x2);
                plot(a4,x2,y2,'k')

                %find possible corner point
                %start corner
                startxe=(horizline-p1(2))/p1(1);
                xline(a4,startxe,'b')
                %end corner
                endxe=(horizline-p2(2))/p2(1);
                xline(a4,endxe,'b') 




                % get the four points surrounding the startxe and endxe & +/- vertstep
                % before and after for context 
                for j=a:b+1
                    if xaxispost(j)<startxe && xaxispost(j+1)>startxe
                        %these are the two points in which the corner line crosses 
                        %the 4 pts of interest are j-1,j,j+1,j+2 
                        cornercontext(k,4:13)=ydetrend((j-1)-vstep:(j+2)+vstep);
                        cornerindex(k,4:13)=xaxispost(j-1-vstep:j+2+vstep);
                        cornerdiff(k,4:13)=ydiff2(j-1-vstep:j+2+vstep);
                    elseif xaxispost(j)<endxe && xaxispost(j+1)>endxe
                        %these are the two points in which the corner line crosses 
                        %the 4 pts of interest are j-1,j,j+1,j+2 
                        cornercontext(k,14:23)=ydetrend((j-1)-vstep:(j+2)+vstep);
                        cornerindex(k,14:23)=xaxispost(j-1-vstep:j+2+vstep);
                        cornerdiff(k,14:23)=ydiff2(j-1-vstep:j+2+vstep);
                    end
                end
                 %what do you think start is? 
                    %we know its either 7,8,9,10
                    if cornerdiff(k,8)>-0.1
                         xline(a4,cornerindex(k,7),'k')
                         guesscheck(k,1)=1;
                    elseif cornerdiff(k,9)>-0.1
                        xline(a4,cornerindex(k,8),'k')
                        guesscheck(k,1)=2;
                    elseif cornerdiff(k,10)>-0.1
                        xline(a4,cornerindex(k,9),'k')
                        guesscheck(k,1)=3;
                    elseif cornerdiff(k,11)>-0.1
                        xline(a4,cornerindex(k,10),'k')
                        guesscheck(k,1)=4;
                    else
                        xline(a4,cornerindex(k,8),'k')
                        guesscheck(k,1)=2;
                        noguess='noguess start'
                    end   
                 %end start guessing

                 %what do you think end is? 
                    %we know its either 17,18,19,20
                    if cornerdiff(k,18)-cornerdiff(k,17)>0.3 && cornerdiff(k,18)>0.1
                        xline(a4,cornerindex(k,17),'k')
                         guesscheck(k,3)=1;
                    elseif cornerdiff(k,19)-cornerdiff(k,18)>0.3 && cornerdiff(k,19)>0.1
                        xline(a4,cornerindex(k,18),'k')
                        guesscheck(k,3)=2;
                    elseif cornerdiff(k,20)-cornerdiff(k,19)>0.3 && cornerdiff(k,20)>0.1
                        xline(a4,cornerindex(k,19),'k')
                        guesscheck(k,3)=3;
                    elseif cornerdiff(k,21)-cornerdiff(k,20)>0.3 && cornerdiff(k,21)>0.1
                        xline(a4,cornerindex(k,20),'k')
                        guesscheck(k,3)=4;
                    else
                        xline(a4,cornerindex(k,18),'k')
                        guesscheck(k,3)=2;
                        noguess='noguess end'
                    end   
                 %end start guessing
                 %label the titles to give user information
                 title(a5,strcat('b1=',num2str(guesscheck(k,1)),' & b2=',num2str(guesscheck(k,3))));
                 x=0;
                 for y=1:k
                     if cornercontext(y,1)==1
                         x=x+1;
                     end
                 end
                 title(a4,strcat('number of cells analyzed=',num2str(x)));
                 
                 
                 %end title labeling
                 [mindetrend,~]=min(ydetrend(a:b));
                 [maxdetrend,~]=max(ydetrend(a:b));
                ylim(a4,[mindetrend,maxdetrend]);
                useranswer3=input('if startguess is right- press enter, otherwise: what position is the real startxe in? (1,2,3,4)');

                if isempty(useranswer3)==0
                    guesscheck(k,2)=useranswer3;
                    cornercontext(k,2)=useranswer3;
                    cornerindex(k,2)=useranswer3;
                else
                    guesscheck(k,2)=-1;
                    cornercontext(k,2)=guesscheck(k,1);
                    cornerindex(k,2)=guesscheck(k,1);
                end
                useranswer4=input('if endguess is right- press enter, otherwise: what position is the real endxe in? (1,2,3,4)');

                if isempty(useranswer4)==0
                    guesscheck(k,4)=useranswer4;
                    cornercontext(k,3)=useranswer4;
                    cornerindex(k,3)=useranswer4;
                else
                    guesscheck(k,4)=-1;
                    cornercontext(k,3)=guesscheck(k,3);
                    cornerindex(k,3)=guesscheck(k,3);
                end
                % end pt gathering
                k=k+1;
                i=i+2; 
                hstep=hstepog;
            else %its a recovery pulse
                %adjust threshold to allow for better recovery pulse
                %identification 
                %make threshold 1/2 of original drop, as long as that is
                %below the original threshold 
                sizing=k-1-numsq;
                squeezing=k-1;
                magsizing=mean([ym(cornerindex(sizing,13));ym(cornerindex(sizing,14))]); %avg value in sizing pulse(between latest start corner adn earliest end corner)
                magsizing=yasls(cornerindex(sizing,14))-magsizing;
                newupper=yasls-(magsizing/2);% setting a new upperlimit for all of the recovery pulses (this should be <the old upperlimit to allow for 
                plot(a5,xaxispost(a-wp:b+wp),newupper(a-wp:b+wp)-yasls(a-wp:b+wp),'b')
                if newupper(cornerindex(sizing,14))<upperlimit(cornerindex(sizing,14)) 
                %if the new upperlimit< old upperlimit then we need to 
                %recalculate where the pulse crosses the new threshold
                    indexofdata=cornerindex(squeezing,14); %look for threshold crossings after squeeze pulse
                    recovstop=1;
                    newthreshpts=zeros(2,2);
                    [newthreshpts]=classicthresholding(yasls3,newupper,newthreshpts,indexofdata,1,lengthpost,recovstop,numrecov);
                    %check if any of the next 10 pulses dont match up with the
                    %new ones
                    match=1;
                    m=2;
                    while m<(numrecov*2)+2 && match==1
                        if abs(threshpts(i+m-2,1)-newthreshpts(m,1))<100
                            match=1;
                            m=m+1;
                        else
                            match=0;
                        end
                    end
                    if m==(numrecov*2)+2
                        threshpts(i:i+numrecov*2-1,:)=newthreshpts(2:(numrecov*2)+1,:);
                    else
                        threshpts(i:i+m-3,:)=newthreshpts(2:m-1,:);
                        threshpts=[threshpts(1:i+m-3,:);newthreshpts(m:(numrecov*2)+1,:);threshpts(i+m-2:end,:)];
                        
                    end
                end   
                [rowthresh,~]=size(threshpts)
                %its time to analyze recovery! 
                [cornercontext,cornerindex,cornerdiff,guesscheck,i,k,hstep]=analyzerecov(a,threshpts,cornercontext,cornerindex,cornerdiff,guesscheck,k,i,yasdet,stderror,yasdiff,xaxispost,ydetrend,ydiff2,a4,a5,numrecov,steps);
            end
    elseif useranswer0==1 
        i=i+2;
    elseif useranswer0==3
        i=i+(numpul*2)-2; %skips the squeeze and recovery pulses
    elseif useranswer0==4
        i=i+(numpul*2); %skips the entire cell event 
    elseif useranswer0==5
        i=i-2;
    elseif useranswer0==6
        i=i-4;
    else
        useranswer=0;
    end
    save(saveto,'cornercontext','cornerindex','cornerdiff','guesscheck','-append');
end
%% post checking
[numr ~]=size(cornerindex);
rcornerindex=zeros(numr,23);
rcornercontext=zeros(numr,23);
rcornerdiff=zeros(numr,23);
for i=1:numr
    %start column indices in the corner matricies are 7,8,9,10 --> 1,2,3,4
    %end column indices in the corner matricies are 17,18,19,20 --> 1,2,3,4
    switch cornerindex(i,2)
        case 0 
            rcornerindex(i,4:5)=0;
            rcornercontext(i,4:5)=0;
            rcornerdiff(i,4:5)=0;
            rcornerindex(i,6:11)=cornerindex(i,4:9);
            rcornercontext(i,6:11)=cornercontext(i,4:9);
            rcornerdiff(i,6:11)=cornercontext(i,4:9);
        case 1
            rcornerindex(i,4)=0;
            rcornercontext(i,4)=0;
            rcornerdiff(i,4)=0;
            rcornerindex(i,5:11)=cornerindex(i,4:10);
            rcornercontext(i,5:11)=cornercontext(i,4:10);
            rcornerdiff(i,5:11)=cornerdiff(i,4:10);
        case 2
            rcornerindex(i,4:11)=cornerindex(i,4:11);
            rcornercontext(i,4:11)=cornercontext(i,4:11);
            rcornerdiff(i,4:11)=cornerdiff(i,4:11);
        case 3
            rcornerindex(i,4:11)=cornerindex(i,5:12);
            rcornercontext(i,4:11)=cornercontext(i,5:12);
            rcornerdiff(i,4:11)=cornerdiff(i,5:12);
        case 4
            rcornerindex(i,4:11)=cornerindex(i,6:13);
            rcornercontext(i,4:11)=cornercontext(i,6:13);
            rcornerdiff(i,4:11)=cornerdiff(i,6:13);
        case 5 
            rcornerindex(i,4:10)=cornerindex(i,7:13);
            rcornercontext(i,4:10)=cornercontext(i,7:13);
            rcornerdiff(i,4:10)=cornerdiff(i,7:13);
            rcornerindex(i,11)=0;
            rcornercontext(i,11)=0;
            rcornerdiff(i,11)=0; 
    end
    switch cornerindex(i,3)
         case 0 
            rcornerindex(i,12:13)=0;
            rcornercontext(i,12:13)=0;
            rcornerdiff(i,12:13)=0;
            rcornerindex(i,14:19)=cornerindex(i,14:19);
            rcornercontext(i,14:19)=cornercontext(i,14:19);
            rcornerdiff(i,14:19)=cornercontext(i,14:19);
        case 1
            rcornerindex(i,12)=0;
            rcornercontext(i,12)=0;
            rcornerdiff(i,12)=0;
            rcornerindex(i,13:19)=cornerindex(i,14:20);
            rcornercontext(i,13:19)=cornercontext(i,14:20);
            rcornerdiff(i,13:19)=cornerdiff(i,14:20);
        case 2
            rcornerindex(i,12:19)=cornerindex(i,14:21);
            rcornercontext(i,12:19)=cornercontext(i,14:21);
            rcornerdiff(i,12:19)=cornerdiff(i,14:21);
        case 3
            rcornerindex(i,12:19)=cornerindex(i,15:22);
            rcornercontext(i,12:19)=cornercontext(i,15:22);
            rcornerdiff(i,12:19)=cornerdiff(i,15:22);
        case 4
            rcornerindex(i,12:19)=cornerindex(i,16:23);
            rcornercontext(i,12:19)=cornercontext(i,16:23);
            rcornerdiff(i,12:19)=cornerdiff(i,16:23);
        case 5 
            rcornerindex(i,12:18)=cornerindex(i,17:23);
            rcornercontext(i,12:18)=cornercontext(i,17:23);
            rcornerdiff(i,12:18)=cornerdiff(i,17:23);
            rcornerindex(i,19)=0;
            rcornercontext(i,19)=0;
            rcornerdiff(i,19)=0; 
    end
        
end
[rguess,~]=size(guesscheck);
indices=zeros(1,numpul*2);
%for startsz--> 1=7, 2=8, 3=9, 4=10 
%--> therefore 6+guesscheck=col# in cornercontext
%for endsz--> 1=17, 2=18, 3=19, 4=20 --> 16+ 
k=1;%the number of cell events
%i=the number of current pulses *2 (where there are 9 pulses per cell event) 
% for i: (1) startsz (2) endsz (3) startsq1 (4) endsq1 (5) startsq2 (6)
% endsq2 (7) startsq3 (8) endsq3 (9) startrecov1 (10) endrecov1 (11)
% startrecov2 (12)endrecov2 (13) startrecov3 (14) endrecov3 (15)
% startrecov4 (16) endrecov4 (17) startrecov5 (18) endrecov5
for i=1:rguess
    if cornercontext(i,1)==1 %is sizing
        [startindex,endindex]=getindices(cornerindex(i,:),guesscheck(i,:));
        indices(k,1)=startindex;
        indices(k,2)=endindex; 
        %now check if the next one is a sq pulse 
        if i+1<=rguess
            if cornercontext(i+1,1)~= 2 %if the next pulse is not a squeeze pulse move onto the next cell event
            k=k+1;
            end
        end
            
    elseif cornercontext(i,1)==2 %is sq pulse
        [startsqindex,endsqindex]=getindices(cornerindex(i,:),guesscheck(i,:));
        %next section should be DEBUGGED
        if cornercontext(i-1,1)==1 && startsqindex-indices(k,2)>100 %this squeeze pulse is actually not right after the previous size pulse 
            k=k+1;
        elseif cornercontext(i-1,1)==2 %could be second or third sq pulse
            if cornercontext(i-2,1)==2 %third squeeze pulse 
                indices(k,7)=startsqindex;
                indices(k,8)=endsqindex;
            else
                indices(k,5)=startsqindex;
                indices(k,6)=endsqindex;
            end
        else
            indices(k,3)=startsqindex;
            indices(k,4)=endsqindex;
        end
        
        if i+1<=rguess 
            if numsq==3 &&  i-2>0
                if cornercontext(i-2,1)==2 && cornercontext(i-1,1)==2 && cornercontext(i+1,1)~= 3 
                %if this is the third sq pulse and the next pulse isn't recovery go to the next indices row (meaning next cell event) 
                    k=k+1;
                end
            elseif numsq==1
                if cornercontext(i+1,1)~= 3 
                    k=k+1;
                end
            end
        end
      
    elseif cornercontext(i,1)==3
        if indices(k,9)==0 
            [startri,endri]=getindices(cornerindex(i,:),guesscheck(i,:)); 
            indices(k,9)=startri;
            indices(k,10)=endri;
        elseif indices(k,11)==0
            [startri,endri]=getindices(cornerindex(i,:),guesscheck(i,:)); 
            indices(k,11)=startri;
            indices(k,12)=endri;
        elseif indices(k,13)==0
            [startri,endri]=getindices(cornerindex(i,:),guesscheck(i,:)); 
            indices(k,13)=startri;
            indices(k,14)=endri;
        elseif indices(k,15)==0
            [startri,endri]=getindices(cornerindex(i,:),guesscheck(i,:)); 
            indices(k,15)=startri;
            indices(k,16)=endri;
        elseif indices(k,17)==0
            [startri,endri]=getindices(cornerindex(i,:),guesscheck(i,:)); 
            indices(k,17)=startri;
            indices(k,18)=endri;
        end
       if i+1<=rguess
            if cornercontext(i+1,1)~= 3 %if the next pulse is not a recov pulse go to next indices row 
            k=k+1;
            end 
        end
    end
end

[numind,~]=size(indices);
x=1;
while x<=numind
    if isempty(find(indices(x,:),1))
        indices(x,:)=[];
        numind=numind-1;
        x=x-1;
    end
    x=x+1;
end
%pulses columns: 1) filename, 2) sizing start, 3) size end, 4) Iavg sz, 5)
%Iavgbaseline, 6)diameter, 7)uflow (sizing), 8) sqstart1, 9) sqend1, 10-13)
%same as 8,9 but for sq2,3, 14-16) vratio1-3, 17-19) strain1-3, 20) tau,
%21) Rsquared, 22) numpast sizing, 23) pressure, 24) geometry, 25) device
%num, 26) voltage, 27) dinput 
pulses=zeros(numind,28);  
[recovavg,recovtimept,fitrecov,pulses,baseymatrix]=calculaterecoveryver6(indices,ym,yasls,stderror,pulses,numrecov,numsq,numpul);
save(saveto,'recovavg','recovtimept','pulses','dinput','-append');

[pulses]=calculateforNPSver6(indices,ym,yasls,sampleRate,N,Length,Deffective,szlength,sqlength,hchannel,baseymatrix,pulses,recovavg,stderror,numpul,info,numsq);

pulsesforprint=cell(numind,28);
pulsesforprint(:,1)=cellstr(name);
[pulses]=recordanalysisinfo(pulses,dinput,info,noise);
pulsesforprint(:,2:end)=num2cell(pulses(:,2:end));
%change device#.0 or .1 to L vs. R 
devicenum=pulses(1,25);
devicenumr=round(devicenum);
rorl=round(devicenum-devicenumr,1);
if rorl==0
    devicenumcell=strcat(num2str(devicenumr),'L');
elseif rorl==0.1
    devicenumcell=strcat(num2str(devicenumr),'R');
end
pulsesforprint(:,25)=cellstr(devicenumcell);
    
    
    

save(saveto,'rcornercontext','rcornerindex','rcornerdiff','pulses','indices','pulsesforprint','-append')



