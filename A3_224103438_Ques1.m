% Birth Month:January
clc
clear all
syms D Ip1 Ip2 L G w

Output = fopen('A3_224103438_Ques1.txt', 'w');
fprintf (Output, 'Two Disc Rotor System\n');
fprintf (Output, 'Using TMM\n');

% Two-Disc Rotor System
% Formula:
J=pi*D^4/32;
Kt=G*J/L;

% Point and Field Matrices
% Disc1:
P1=[1 0;-w^2*Ip1 1];

% Disc2:
P2=[1 0;-w^2*Ip2 1];

% Shaft:
F=[1 1/Kt;0 1];

% Overall Transfer Matrix
T=P2*F*P1;

% Input Data and Boundary Condition:
T=subs(T,{D,Ip1,Ip2,L,G},{0.01,0.0008,0.002,0.075,0.8*10^11});
f=T(2);
f=solve(f,w);

% Natural Frequency:
wnf=vpa(f)

% Finding Mode Shapes:

M=[Ip1,0;0 Ip2];
K=[Kt Kt;Kt Kt];
A=inv(M)*K;

% Data Input:
A=subs(A,{D,Ip1,Ip2,L,G},{0.01,0.0008,0.002,0.075,0.8*10^11});
[eig_vec eig_value]=eig(A);


% Mode Shapes:
mode_shape=eig_vec


% Writing Result in Output file:
%Printing the results to output file
P(:,:,1) = P1;
P(:,:,2) = P2;
F(:,:,:) = F;


fprintf (Output,'Point matrices :\n');
for i = 1:2
    fprintf (Output,'[P]%d :\n',i);
    for a=1:2
        for b = 1:2
            fprintf (Output,'%s\t\t',char(P(a,b,i)));
        end
        fprintf (Output,'\n');
    end
    fprintf (Output,'\n');
end
fprintf (Output,'Field matrix :\n');
for i=1:1
        fprintf (Output,'[F]%d :\n',i);
        for a=1:2
            for b = 1:2
                fprintf (Output,'%s\t\t',char(F(a,b,i)));
            end
            fprintf (Output,'\n');
        end
        fprintf (Output,'\n');
end

fprintf (Output,'Free-Free boundary condition\n\n');
fprintf(Output,'t(2,1)=0\n\n');
fprintf(Output,'Natural Frequecies:\n\n');
fprintf(Output,'%f\n',wnf);
fprintf(Output,'Mode Shapes:\n\n');
fprintf(Output,'%f\n',mode_shape);
fclose(Output);
