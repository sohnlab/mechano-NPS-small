function [startindex,endindex]=getindices(cornerindex,guesscheck)
    if guesscheck(1,2) ==-1
        startindex=guesscheck(1,1);
        startindex=cornerindex(1,startindex+6);
    else 
        startindex=guesscheck(1,2);
        startindex=cornerindex(1,startindex+6);
    end
    if guesscheck(1,4)==-1
            endindex=guesscheck(1,3);
            endindex=cornerindex(1,endindex+16);
    else
            endindex=guesscheck(1,4);
            endindex=cornerindex(1,endindex+16);
    end
    