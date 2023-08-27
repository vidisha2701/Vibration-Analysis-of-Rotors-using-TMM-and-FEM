% Roll Number: 224103438
% Birth Month: January

clc
clear all
syms D L G w E rho 

A=pi*D^2/4;              % m^2
I=pi*D^4/64;             % m^4

disp('FINITE ELEMENT METHOD FORMULATION')

% Elemental equation for Element 1:
m1=(rho*A*L/420)*[156 22*L 54 -12*L;22*L 4*L^2 13*L -3*L^2;54 13*L 156 -22*L;-12*L -3*L^2 -22*L 4*L^2];
m1=subs(m1,{D,L,E,rho},{0.05,0.5,2.1*10^11,7800});
k1=(E*I/L^3)*[12 6*L -12 6*L;6*L 4*L^2 -6*L 2*L^2;-12 -6*L 12 -6*L;6*L 2*L^2 -6*L 4*L^2];
k1=subs(k1,{D,L,E,rho},{0.05,0.5,2.1*10^11,7800});
k1=double(k1);
m1=double(m1);

% Elemental Equation for Element 2:
m2=(rho*A*L/420)*[156 22*L 54 -12*L;22*L 4*L^2 13*L -3*L^2;54 13*L 156 -22*L;-12*L -3*L^2 -22*L 4*L^2];
m2=subs(m2,{D,L,E,rho},{0.05,0.5,2.1*10^11,7800});
k2=(E*I/L^3)*[12 6*L -12 6*L;6*L 4*L^2 -6*L 2*L^2;-12 -6*L 12 -6*L;6*L 2*L^2 -6*L 4*L^2];
k2=subs(k2,{D,L,E,rho},{0.05,0.5,2.1*10^11,7800});
k2=double(k2);
m2=double(m2);

% Global Mass Matrix
Global_M=zeros(6,6);
Global_M(1:4,1:4)=m1;
Global_M(3:6,3:6)=Global_M(3:6,3:6)+m2

% Global Stiffness Matrix
Global_K=zeros(6,6);
Global_K(1:4,1:4)=k1;
Global_K(3:6,3:6)=Global_K(3:6,3:6)+k2

% Reduced Mass and Stiffness Matrix
Global_M=Global_M(3:4,3:4);
Global_K=Global_K(3:4,3:4);

% Applying Boundary Condition
Global_M=(-w^2*Global_M);
A2=(Global_M + Global_K);

% Solution: Natural Frequency
A2=det(A2);
wnf=solve(A2,w);
wnf=vpa(wnf([2,4]))






