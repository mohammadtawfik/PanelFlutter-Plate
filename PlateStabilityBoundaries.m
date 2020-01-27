%This program finds the stability boundaries for  
% a panel subject to Aerodynamic and thermal
% loading using C1 rectangular elements

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

%Clearing the memory
clear all
close all
clc

%Problem Data
%Aluminum material properties are used
Modulus  = 71e9;   %GPa    - Modulus of elasticity
Nu       = 0.3;    %         Poisson's ration
Rho      = 2700;   %Kg/m^3 - Density
Alpha    = 22.5e-6;%         Thermal expansion coefficient
%Geometric Data
LengthX  = 1;         % m
LengthY  = 1;         % m
Thickness= 0.001;     % m
%Finite Element Problem Data
Nx = 5; %number of elements in the x-direction
Ny = 5; %number of elements in the x-direction
%Support conditions are user-defined
% 0=CCCC
% 1=SSSS
SupportConditions=1;

%Evaluating the basic quantities
Lx=LengthX/Nx; %element dimensions
Ly=LengthY/Ny;
Qq=Modulus/(1-Nu*Nu)* ...
    [1,Nu,0;Nu,1,0;0,0,(1-Nu)/2];
Dd=Qq*Thickness*Thickness*Thickness/12;
Nt=Modulus*Thickness*Alpha/(1-Nu);
Wo=sqrt(Dd(1,1)/Rho/Thickness/LengthX/LengthX/LengthX/LengthX);
TbInv=CalcTbInv(Lx,Ly);

[Kb,Mb,Kt,Ka]=CalcLinear(Dd,Lx,Ly);
Kb=TbInv'*Kb*TbInv;
Ga=TbInv'*Mb*TbInv;
Mb=Ga*Rho*Thickness;
Ga =Dd(1,1)/Wo/(LengthX^4)*Ga; 
Kt=Nt*TbInv'*Kt*TbInv;
Ka=Dd(1,1)*TbInv'*Ka*TbInv/LengthX/LengthX/LengthX;
%***************************************
%Creating the Mesh (Noedes and Elements)
%***************************************
%Initializing the Nodes registry!
% Columns 1&2 contain the 
%  coordinates of the nodes (x,y)
% Columns 3 to 6 contain the numbers
%  of the degrees of freedom at the node
Nodes=zeros((Nx+1)*(Ny+1),6);
%filling the Nodes registry
for jj=1:Ny+1
  for ii=1:Nx+1
    Pointer=(jj-1)*(Nx+1) + ii; %Node number
    Nodes(Pointer,1)=(ii-1)*Lx;
    Nodes(Pointer,2)=(jj-1)*Ly;
    Nodes(Pointer,3)=Pointer*4-3; %w
    Nodes(Pointer,4)=Pointer*4-2; %wx
    Nodes(Pointer,5)=Pointer*4-1; %wy
    Nodes(Pointer,6)=Pointer*4;   %wxy
  endfor
endfor

%Initializing the Elements registry
% Columns 1 to 4 contain the global node
%  numbers associated with the element
Elements=zeros(Nx*Ny,4);
for jj=1:Ny
  for ii=1:Nx
    Pointer=(jj-1)*Nx + ii;
    Elements(Pointer,1)=(jj-1)*(Nx+1)+ii;
    Elements(Pointer,2)=(jj-1)*(Nx+1)+ii+1;
    Elements(Pointer,3)=(jj)*(Nx+1)+ii+1;
    Elements(Pointer,4)=(jj)*(Nx+1)+ii;
  endfor
endfor
%
%Creating the Boundary Conditions vector
% This vector will contain the numbers of
%  the fixed degrees of freedom;
%  ie. the columns and rows that will be
%  eleminated from the global matrix
BCs=[];
for jj=1:Ny+1
  for ii=1:Nx+1
    Pointer=(jj-1)*(Nx+1) + ii; %Node number
    %If the node is on the sides ...
    if or(or(jj==1,jj==Ny+1),or(ii==1,ii==Nx+1))
      if SupportConditions==0
        %The following line counts the BCs
        % for a CCCC plate
        BCs=[BCs,Nodes(Pointer,3:6)];
      elseif SupportConditions==1
        %The following lines counts for BCs
        % for SSSS plate
        %if the node is a corner node
        if and(or(jj==1,jj==Ny+1),or(ii==1,ii==Nx+1))
          BCs=[BCs,Nodes(Pointer,3:6)];
        %if the node is on either horizontal sides
        elseif or(jj==1,jj==Ny+1)
          BCs=[BCs,Nodes(Pointer,3:4)];
        %if the node is on either vertical sides
        elseif or(ii==1,ii==Nx+1)
          BCs=[BCs,Nodes(Pointer,3)];
          BCs=[BCs,Nodes(Pointer,5)];
        endif
      endif
    endif
  endfor
endfor
%The complementary coundary conditions
% is a vector that contains the numbers
% of the degrees of freedom that are FREE
BCsC=1:4*(Nx+1)*(Ny+1);
%BCs
BCsC(BCs)=[];

%*****************************************
%The global Matrices
%*****************************************
%Initialization
KG=zeros(4*(Nx+1)*(Ny+1),4*(Nx+1)*(Ny+1));
MG=zeros(4*(Nx+1)*(Ny+1),4*(Nx+1)*(Ny+1));
KT=zeros(4*(Nx+1)*(Ny+1),4*(Nx+1)*(Ny+1));
KA=zeros(4*(Nx+1)*(Ny+1),4*(Nx+1)*(Ny+1));
GA=zeros(4*(Nx+1)*(Ny+1),4*(Nx+1)*(Ny+1));

%Looping for the elements
for ii=1:Nx*Ny
  %Looping for the element nodes
  for jj=1:4
    DOFsj=Nodes(Elements(ii,jj),3:6);
    %Looping AGAIN for the elements nodes;
    for kk=1:4
      DOFsk=Nodes(Elements(ii,kk),3:6);
      KG(DOFsj,DOFsk)=KG(DOFsj,DOFsk)+ ...
        Kb(4*jj-3:4*jj,4*kk-3:4*kk);
      KT(DOFsj,DOFsk)=KT(DOFsj,DOFsk)+ ...
        Kt(4*jj-3:4*jj,4*kk-3:4*kk);
      MG(DOFsj,DOFsk)=MG(DOFsj,DOFsk)+ ...
        Mb(4*jj-3:4*jj,4*kk-3:4*kk);
      GA(DOFsj,DOFsk)=GA(DOFsj,DOFsk)+ ...
        Ga(4*jj-3:4*jj,4*kk-3:4*kk);
      KA(DOFsj,DOFsk)=KA(DOFsj,DOFsk)+ ...
        Ka(4*jj-3:4*jj,4*kk-3:4*kk);
    endfor
  endfor
endfor
%Applying the boundary conditions
KReduced=KG(BCsC,BCsC);
TReduced=KT(BCsC,BCsC);
MReduced=MG(BCsC,BCsC);
GReduced=GA(BCsC,BCsC);
AReduced=KA(BCsC,BCsC);


%Evaluating the buckling temperature
% at different air-flow speeds (dynamic pressure)
DLamda=11;
for ii=0:20
  Lamda=ii*DLamda;
  LamdaVect(ii+1)=Lamda;
  BuckTemp(ii+1)= ...
    min(eig(inv(TReduced)*(KReduced+Lamda*AReduced))) ...
    *12*(1+Nu)*Alpha/pi/pi/Thickness/Thickness;
endfor
%******************************
%Evaluating the flutter speed
% by searching for it at different temperatures
%******************************
%Evaluating the buckling temperature
ttt=min(eig(inv(TReduced)*KReduced));
DeltaT=ttt/10; %Temperature increment
JInitial=300;
%Looping for temperatures
for ii=0:9
  Ttt=ii*DeltaT*2; %Current Temperature
  TVect(ii+1)=Ttt*12*(1+Nu)*Alpha/pi/pi/Thickness/Thickness;
  %Looping for Dynamic pressure
  for jj=JInitial:1000
    Lamda=jj;
    %Total Stiffness Matrix
    KTotal=KReduced-Ttt*TReduced+Lamda*AReduced;
    %Eigenvalues
    Kappa=eig(inv(MReduced)*KTotal);
    %Checking if flutter ocured
    %Due to numerical inacurecies, the imaginary
    % values may be in the order of 10^-5
    % that is why we set the threshold at 0.001
    % instead of zero
    if max(imag(Kappa))>0.001
      JInitial=Lamda-100; %resetting the initial pressure
      %Storing the flutter pressure
      LamdaFlutter(ii+1)=Lamda;
      break
    endif
  endfor
endfor

figure(2)
plot(BuckTemp,LamdaVect,TVect,LamdaFlutter)
grid
xlabel("Normalized Temp")
ylabel("Nondimensional Dynamic Pressure")