%=================================================================
% Christophe Hurlin
% April, 2011
% LEO, Universit� d'Orl�ans
%================================================================

function f1=RunMyCode_GMM(VaR,r,alpha,pmax) 

warning off

%=================
%=== Controls  ===
%=================

if (abs(length(VaR)-length(r))>0) | (alpha>0.10) | (alpha<0.001)

    if abs(length(VaR)-length(r))>0 
        
        disp(' ') ,disp('    Check your data: the VaR and r series do not have the same length')  
        
        f1=NaN;
        
    else
        disp(' ') ,disp('    Check your data: alpha must range between 0.001 and 0.10 (excluded) ')  
        
        f1=NaN;
        
    end
    
    
else
    
%==================
%=== Parameters ===
%==================

T=length(r);    % sample sizes

ns=9999;
      
%==============================================
% ====Calcul Stat LR and Stat BM under H0 ====
%==============================================

statgeoH0=NaN*ones(ns,pmax);                         % Vector of TDA (UC and CC) stat under H0        

statgeoH0_IND=statgeoH0;                             % Vector of TDA (independence) stat under H0          

Ui=unifrnd(0,1,ns,1);                                 % Uniform random number
                
s=0;                                                  % Simulation number
                
while s<ns                                            % While Loop on Simulation number

    s=s+1;                                            % Simulation Number
            
    I=binornd(1,alpha,T,1);                           % Hit function under H0
                 
    res=Duree(I);                                     % Durations
                 
    d=res.duree;                                      % Durations
                 
    cens=res.censure;                                 % Cenored Durations 
            
    Y=d-1;                                            % Number of Failures 

    beta_hat=sum(I)/T;
    
    for indic_p=1:pmax                                % Loop on Laguerre Polynomial
    
        p=indic_p;                                     % Number of Laguerre Polynomial
               
        res2=TDA_Geometric(Y,alpha,indic_p);           % TDA Statistic (Geometric Distribution)
                 
        statgeoH0(s,indic_p)=res2.tda;                 % TDA Statistic
            
        res3=TDA_Geometric_IND(Y,beta_hat,indic_p);           % TDA Statistic (Geometric Distribution)
                 
        statgeoH0_IND(s,indic_p)=res3.tda;                 % TDA Statistic
         
    end

    if sum(isfinite(statgeoH0(s,:)))<pmax;                     % Selection Rate
                    
        s=s-1;                                    % If the computation is not done s=s-1
                    
    end                                             % End of if condition Selection Rate
                                     
end                                                   % Loop on ns

        
%========================
%==== Duration series ===
%========================
       
I=(r<VaR);                                           % Hit variable
                                
res=Duree(I);                                           % Computation of Durations
                
d=res.duree;                                            % Durations
                
cens=res.censure;                                       % Binomial variable for cesnure
   
Y=d-1;                                                  % Number of Failures             

beta_hat=sum(I)/T;                                      % Estimation of beta wit historical data
       
%==================================================
%==== Calcul Stat TDA (Geometric Distribution )====
%==================================================
    
stat_tda=zeros(1,pmax);

pvalue_tda=stat_tda;

pvalue_tda_corrected=stat_tda;

for indic_p=1:pmax                                      % Loop on  Number of Laguerre Polynomials 
                
    U0=unifrnd(0,1);                                        % Unifrom randm number 

    p=indic_p;                                          % Number of Laguerre Polynomials
    
    res=TDA_Geometric(Y,alpha,p);                       % TDA Test (Geometric Distribution)
                 
    pvalue_tda(indic_p)=res.pvalue;                    % pvalue 
                
    GN_geo=1-(1/ns)*sum(statgeoH0(:,indic_p)<=res.tda)+(1/ns)*sum((statgeoH0(:,indic_p)==res.tda).*(Ui>=U0));
                
    pvalue_tda_corrected(indic_p)=(ns*GN_geo+1)/(1+ns); % Corrected TDA statistic
                
    stat_tda(indic_p)=res.tda;
                
end                                                      % End of Loop on  Number of Laguerre Polynomials 
     
%========================================
%==== Calcul Stat TDA (Independence )====
%========================================

stat_tda_IND=zeros(1,pmax);

pvalue_tda_IND=stat_tda;

pvalue_tda_corrected_IND=stat_tda;

for indic_p=1:pmax                                      % Loop on  Number of Laguerre Polynomials 
                
    U0=unifrnd(0,1);                                    % Unifrom randm number 

    p=indic_p;                                          % Number of Laguerre Polynomials
    
    res=TDA_Geometric_IND(Y,beta_hat,p);                       % TDA Test (Geometric Distribution)
                 
    pvalue_tda_IND(indic_p)=res.pvalue;                    % pvalue 
                
    GN_geo_IND=1-(1/ns)*sum(statgeoH0_IND(:,indic_p)<=res.tda)+(1/ns)*sum((statgeoH0_IND(:,indic_p)==res.tda).*(Ui>=U0));
                
    pvalue_tda_IND_corrected(indic_p)=(ns*GN_geo_IND+1)/(1+ns); % Corrected TDA statistic
                
    stat_tda_IND(indic_p)=res.tda;
                
end                                                      % End of Loop on  Number of Laguerre Polynomials 
     
%==================
%=== Graphiques ===
%==================

e=(1:1:length(r))';

f1=figure;

plot(e,r,e,VaR)

legend ('Returns','Value-at-Risk',0) 

%================
%== Affichage ===
%================

disp (' ===============================================================')
disp (' ==   Backtesting Value-at-Risk: A GMM Duration Based Test    ==')
disp (' == Candelon B., C. Colletaz, C. Hurlin and S. Tokpavi,(2011) ==' )
disp (' ===============================================================')

disp(' ')

disp(sprintf('  Coverate Rate of Value-at-Risk = %1.2f',alpha)), disp(' ')

disp(sprintf('  Number of violations = %2.0f',sum(I))), disp(' ')

disp(sprintf('  Empiricial frequency of violations = %1.4f',sum(I)/T)), disp(' ')

disp(sprintf('  Maximal number of polynomials = %1.0f',pmax)), disp(' ')

disp(sprintf('  Sample size = %4.0f',T)), disp(' ')

disp(sprintf('  Number of simulations (for Dufour''s correction)= %1.0f',ns)), disp(' ')

disp(' **************************************************************')
disp(' *** Conditional and Unconditional (p=1) Coverage (UC) Tests **')
disp(' **************************************************************')
disp(' ')

for i=1:pmax
    
    if i==1
        
        disp(sprintf('  p = %1.0f   J_UC statistic = %2.4f   J_UC pvalue (corrected) = %2.4f ' ,i, stat_tda(i),pvalue_tda_corrected(i) )),disp(' ')

    else

        disp(sprintf('  p = %1.0f   J_CC statistic = %2.4f   J_CC pvalue (corrected) = %2.4f ' ,i, stat_tda(i),pvalue_tda_corrected(i) )),disp(' ')

    end
    
end


disp(' ')
disp(' *******************************')
disp(' *** Independence (IND) Tests **')
disp(' *******************************')
disp(' ')

for i=1:pmax

    disp(sprintf('  p = %1.0f   J_IND statistic = %2.4f   J_IND pvalue (corrected) = %2.4f ' ,i, stat_tda_IND(i),pvalue_tda_IND_corrected(i) )),disp(' ')

end

end
