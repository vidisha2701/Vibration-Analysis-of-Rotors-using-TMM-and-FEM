% Assignment 4:224103438
% Birth Month:January

clc
clear all
syms D1 L1 D2 L2 G w Ip z t

% Using Transfer Matrix Method

% Polar MOI of Shaft in (m^4)
J1=pi*D1^4/32;
J2=pi*D2^4/32;

% Stiffness of Shaft in (N/m)
kt1=G*J1/L1;
kt2=G*J2/L2;

% Point and Field Matrix
P=[1 0;-w^2*Ip 1];
F1=[1 1/kt1;0 1];
F2=[1 1/kt2;0 1];

% Applying Boundary Condition:
s=[z;t];
A=P*F1*F2*s;
A1=subs(A,{z},{0});
A1=subs(A,{D1,L1,D2,L2,G,Ip,z},{0.01,0.6,0.03,0.4,8*10^10,0.01,0});

% SOLUTION
wnf= solve(A1(2),w);
disp('TRANSFER MATRIX METHOD')
fprintf('\n')
disp('The Frequency Using Transfer Matrix Method :')
wnf=vpa(wnf(2))


% Using Finite Element Method:

disp('Using Finite Element Method')

% Element 1
M1=[0 0;0 0];
k1=(G*J2/L2)*[1 -1;-1 1];
k1=subs(k1,{G,D2,L2},{8*10^10,0.03,0.4});

% Element 2
k2=(G*J1/L1)*[1 -1;-1 1];
M2=[0 0;0 Ip];
M2=subs(M2,Ip,0.01);
k2=subs(k2,{G,D1,L1},{0.8*10^11,0.01,0.6});

% Global Mass matrix
M=zeros(3,3);
M(1:2,1:2)=M1;
M(2:3,2:3)=M(2:3,2:3)+M2

% Global Stiffness matrix
K=zeros(3,3);
K(1:2,1:2)=k1;
K(2:3,2:3)=K(2:3,2:3)+k2

% Boundary condition
M=vpa(-w^2*M);
A2=vpa(M+K);

% Reduced Mass and Stiffness mastrix
A2=A2(2:3,2:3);

%solving for frequency
A2=det(A2);
wnf2=solve(A2,w);
disp('The Frequency using Finite Element Method:')
wnf2=wnf2(2)










