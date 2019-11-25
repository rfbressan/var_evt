
%== GMM : version 1.00 
%== 13 avril 2011

clc

clear all

close all

warning off

%=== data demo ===

data=xlsread('DataLaCAixa.xls');

VaR=data(:,2);

r=data(:,1);

alpha=0.01;

pmax=6;

% === Appel de fonction ==

f1=RunMyCode_GMM(VaR,r,alpha,pmax) ;


