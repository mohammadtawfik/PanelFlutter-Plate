%This function evaluates linear matrices
% [Kb]  [Mb]  [Kt] [Ka]
%It may be used for the finite element 
% program for plates of C1 continuity

%Functions will work on Octave, FreeMat
% and Matlab
%Created by Mohammad Tawfik
%mohammad.tawfik@gmail.com 
%In assotiation with research papers
% published on ResearchGate.Net
%Author: Mohammad Tawfik

%Title: Panel Flutter
%DOI: 10.13140/RG.2.1.1537.6807
%Updated text link:
%https://www.researchgate.net/publication/275712979_Panel_Flutter

%More code and functions related to
% this problem may be downloaded from:
% https://github.com/mohammadtawfik/PanelFlutter-Plate

function [KB,MB,KT,KA]=CalcLinear(D,Lx,Ly)

GCn=8; %Number of Gauss integration points
%Get the Gauss constants
GaussConstants=GetGC(GCn);
%Initialize the matrices
KB=zeros(16,16);
MB=zeros(16,16);
KT=zeros(16,16);
KA=zeros(16,16);

%Start the numerical integrration procedure
for Xi=1:GCn
  %Evaluating the physiscal value of X
  X = Lx * (GaussConstants(2, Xi) + 1) / 2;
  for Yi=1:GCn
    %Evaluating the physical value of Y
    Y= Ly * (GaussConstants(2,Yi) + 1) / 2;
    Hw =CalcHw(X,Y);
    Hwx=CalcHwx(X,Y);
    Cb =CalcCb(X,Y);
    Ct =CalcCt(X,Y);
    
    Kb= Cb'*D*Cb;
    Mb= Hw'*Hw;
    Kt= Ct'*Ct;
    Ka= Hw'*Hwx;
    %performing the weighted summation
    KB=KB+GaussConstants(1,Xi)*GaussConstants(1,Yi)*Kb;
    MB=MB+GaussConstants(1,Xi)*GaussConstants(1,Yi)*Mb;
    KT=KT+GaussConstants(1,Xi)*GaussConstants(1,Yi)*Kt;
    KA=KA+GaussConstants(1,Xi)*GaussConstants(1,Yi)*Ka;
    %End of Calculation loop body
  endfor
endfor 
%Multiplying the resulting matrices by Jacobian
KB = KB * Lx * Ly / 4;
MB = MB * Lx * Ly / 4;
KT = KT * Lx * Ly / 4;
KA = KA * Lx * Ly / 4;
endfunction
