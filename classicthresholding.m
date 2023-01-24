function [threshpts]=classicthresholding(yasls3,upperlimit,threshpts,i,k,lengthpost,recovstop,numrecov)
%i is the index of the data
%k is the index of the threshpts
%upperlimit is the line that yasls3 must pass to create a threshpts
%recovstop is only important if you are using this algorithm to RECALCULATE
%the threshpts for recovery pulses (during recovery pulse analysis) 
while i<lengthpost-2 && recovstop~=2
    if yasls3(i)==upperlimit(i) % if a data point== the threshold
        threshpts(k,1)=i;
%         xline(a2,i,'k');
        %determine whether trending down or up? 
        if yasls3(i-1)>upperlimit(i) && yasls3(i+1)<upperlimit(i+1) %trending down 
            threshpts(k,2)=0;
        elseif yasls3(i-1)<upperlimit(i) && yasls3(i+1)>upperlimit(i+1) %trending up 
            threshpts(k,2)=1;
        else
            threshpts(k,2)=3;
        end
            
        k=k+1;
    elseif yasls3(i)>upperlimit(i) && yasls3(i+1)<upperlimit(i+1) && yasls3(i+2)<upperlimit(i+2) % if btwn two points going down
        threshpts(k,1)=i;
        threshpts(k,2)=0;
%         xline(a2,i,'k');
        k=k+1;
    elseif yasls3(i)<upperlimit(i) && yasls3(i+1)>upperlimit(i+1) && yasls3(i-1)<upperlimit(i-1) % if btwn two points going up 
        threshpts(k,1)=i+1;
        threshpts(k,2)=1;
%         xline(a2,i+1,'k');
        k=k+1;
    end
    if recovstop==1
        %this stops the function after 10 pulses if this is a recovery
        %round
        if k==(numrecov+1)*2
            recovstop=2;
        end
    end
    i=i+1;
end
