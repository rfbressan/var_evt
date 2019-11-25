% =========================================================================
% ===  DUREE
% ==== PURPOSE: Compute the duration series from the sequences of 
%          violations as in Christoffersen & Pelletier (2004)
% =========================================================================
%  Function : DUREE(I)
%
%  where - I is the vector (T,1) of violations series
% =========================================================================
%  Output : - res.duree
%  Tokpavi Sessi and C. Hurlin
%  July, 2007. 
%  LEO, University of Orleans
% =========================================================================

function [res]=Duree(I)

if sum(I)==0
 
    res.duree=NaN;             % Cas with 0 violations
    
    res.censure=NaN;
    
else
    
    nb=sum(I)-(I(1)==1)-(I(end)==1)+1;       % Number of Duration

    duree=ones(nb,1)*NaN;

    t=(1:1:length(I))';

    date_violation=t(I==1);
    
    if (sum(I)==1)==1

        if (I(1)==1)==1
            
            duree(1)=0;
            
        else

            duree(1)=date_violation(1)-1;
            
        end
        
    else
        
        if (I(1)==1)==1
        
        duree(1)=date_violation(2)-1;
        
        else
        
        duree(1)=date_violation(1)-1;
        
        end
        
    end
        
    for i=2:nb-1
        
        j=i+(I(1)==1);
    
        duree(i)=date_violation(j)-date_violation(j-1);
    
    end
    
    if (I(end)==1)==0

    duree(end)=length(I)-date_violation(end);
    
    else
        
        if (sum(I)==1)==1
            
        duree(end)=date_violation(end);
        
        else
            
            duree(end)=date_violation(end)-date_violation(end-1);
            
        end
                
    end

    res.duree=duree;                 

    censure=zeros(length(duree),1);

    censure(1)=(I(1)==0);

    censure(end)=(I(end)==0);

    res.censure=censure;

end

