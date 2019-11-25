%=========================================================================
% ==== PURPOSE: Test for Distributional Asumptions 
%               "Moments Based Tests for Discrete Distributions", Bontemps (2006)
%               Independence Test 
% ========================================================================
%
%  Function : TDA_Geometric(Y,beta,p)
%
%  where - Y is the vector of number of  failures before the first Hit supported on the set { 0, 1, 2, 3, ... }. 
%        - p the number of orthogonal conditions (see Bontemps, 2006)
%        - beta = estimated failure rate
% ========================================================================
%  Output : - res.tda:  The statistic of the test
%           - res.pvalue : The corresponding P-value 
% =================================================================== 
% Hurlin Christophe 
% August, 2007. 
% LEO, University of Orleans
% ==================================================================

function [res]=TDA_Geometric_IND(Y,beta,p)
warning off
%=========================================================
%=== The p equations of Meixner Polynomials 
%=== for geometric distribution under the null
%==== These equations defines the orthogonality conditions
%=========================================================

alpha=beta;                            % Estimated failure rate
    
mu=(1-alpha)/alpha;                     % Transformed parameter            

nobs=length(Y);                         % Sample size

m=zeros(nobs,p);                        % Matrix of Residuals

m1=0;                                   % Initial Condition M-1(y)

m0=1;                                   % Initial Condition M0(y)       

m(:,1)=(mu-Y)/sqrt(mu*(mu+1));          % M1(y)

if p>1                                 
    
    j=1;m(:,2)=(mu*(2*j+1)+(j-Y))/((j+1)*sqrt(mu*(mu+1))).*m(:,1)-(j/(j+1)).*m0;
    
    for j=2:p-1
    
        m(:,j+1)=(mu*(2*j+1)+(j-Y))/((j+1)*sqrt(mu*(mu+1))).*m(:,j)-(j/(j+1)).*m(:,j-1);
            
    end

end

stat=(((1/sqrt(nobs))*sum(m)).^2)';

res.tda=sum(stat);

res.pvalue=1-chi2cdf(res.tda,p);




