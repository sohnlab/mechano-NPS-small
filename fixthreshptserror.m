function [threshpts,rowthresh]=fixthreshptserror(threshpts,interval,i,rowthresh)

if threshpts(i,2)==1
%     xline(a2,threshpts(i),'b');
    threshpts(i,:)=[];
    rowthresh=rowthresh-1;
end
while i<rowthresh
    if threshpts(i+1,1)-threshpts(i,1)<interval
        if threshpts(i,2)==0 && threshpts(i+1,2)==1
%             xline(a2,threshpts(i),'b');
%             xline(a2,threshpts(i+1),'b');
            threshpts(i,:)=[];
            threshpts(i,:)=[];
            rowthresh=rowthresh-2;
        elseif threshpts(i,2)==1     %idk what this section is?? it might be bad start
%             xline(a2,threshpts(i),'b');
            threshpts(i,:)=[];
            rowthresh=rowthresh-1;
        end
    else
        i=i+2;
    end
    
end