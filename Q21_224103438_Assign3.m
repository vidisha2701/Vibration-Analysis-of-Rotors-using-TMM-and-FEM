% Vidisha Singh:224103438

clc
clear all

% Three disc stepped shaft system
% Given System parameters
%All Polar mass moment of inertia values are in kg-m^2
Ip1 = 0.03;       
Ip2 = 0.02;
Ip3 = 0.04;

%All shaft diameter values are in m.
da = 0.04;         
db = 0.1;
dc = 0.1;
d = [da, db, dc];

%Polar second moment of area in m^4.
J = pi*(d.^4)/32; 

%All lengths of the shafts in m.
la = 0.6;             
lb = 0.4;
lc = 0.4;

%Rigidity modulus of shafts in N/m^2.
G = 80*10^9;         

Output = fopen('Q21_224103438_Assign3.txt', 'w');
fprintf (Output, 'Given system is a 3 disc stepped shaft system\n\n');
fprintf (Output, 'Using Direct Approach\n');

disp('Using Direct Method')
% Finding the length of equivalent uniform shaft system between disc 1 and 2
le = (J(2)/J(1))*la + lb;
fprintf (Output, ['The length of equivalent uniform shaft system of diameter 0.4 m. between disc 1 and disc 2 is ' num2str(le) ' m.\n']);

%Shaft stiffness calculation
Kte = G*J(2)/le;    
Ktc = G*J(3)/lc;

%Solving for eigen values and eigen vectors
%Inertia matrix
M = [Ip1, 0, 0; 0, Ip2, 0; 0, 0, Ip3];  

%Stiffness matrix
K = [Kte, -Kte, 0; -Kte, Kte+Ktc, -Ktc; 0, -Ktc, Ktc];  

D = M\K;
[p, q] = eig(D);

%Natural frequency of the system
Wnf = sqrt(q)    

%Mode shapes of the system
p(:,1) = p(:,1)/p(3,1);      
p(:,2) = p(:,2)/p(3,2);
p(:,3) = p(:,3)/p(3,3);
Modeshapes = p

fprintf (Output, 'Obtained natural frequencies and mode shapes are shown in pdf copy of solution.\n');
fprintf (Output, 'One of the natural frequency is zero, corresponding to rigid body mode, where all the discs will have same angular displacement and hence no nodes.\n');
fprintf (Output, 'Of the remaining two flexible modes, one has single node and 2nd one has two nodes.\n');
fprintf (Output, 'Locations of these nodes are calculated and shown pictorially in pdf solution copy.\n');
disp('Node locations')

% For Wnf2
lne1_2 = G*J(2)/((Wnf(2,2)^2)*Ip1)        %Location of node in equivalent shaft system from disc 1 to the right of it
% For Wnf3
lne1_3 = G*J(2)/((Wnf(1,1)^2)*Ip1)        %Location of node 1 in equivalent shaft system from disc 1 to the right of it
lne2_3 = G*J(3)/((Wnf(1,1)^2)*Ip3)

    %%%   Using Indirect Method    %%%
disp('Using Indirect method')
fprintf (Output, '\n\nUsing Indirect method\n');
syms x
a = Ip3 + ((Ip2+Ip3)*Ip3)/Ip1;
b = Ip2*le + Ip3*(le+lc) + Ip3*Ip2*lc/Ip1;
c = Ip2*le*lc;
sol = solve(a*(x^2)- (b*x) + c, x);
digits(6);
disp('Node locations are obtained as below:');
l3_2 = vpa(sol)         %Node locations in equivalent shaft system
l1_2 = l3_2*Ip3/Ip1
l3_2 = double(l3_2);
fprintf (Output,['The values of l3(2) are obtained as ' num2str(l3_2(1)) ' m. and ' num2str(l3_2(2)) ' m.\n']);
fprintf (Output, 'Based on above values the node locations are pictorially shown in pdf solution copy.\n');
% Mode shapes based on the concept of similar triangles
%For flexible mode having two nodes
disp('For mode shape having two nodes');
phiz3 = 1             %Taking unit displacement for disc 3
l1 = le;
l2 = lc;
phiz2 = (l2-l3_2(1))/l3_2(1)
phiz1 = (phiz2*l1_2(1))/(l1-l1_2(1))
%Natural frequencies
Wnf1 = 0
Wnf2 = sqrt(G*J(2)/(l1_2(2)*Ip1))
Wnf3 = sqrt(G*J(2)/(l1_2(1)*Ip1))
fprintf (Output, 'Mode shapes and natural frequencies obtained are same as that of Direct method.\n');

      %%% TMM %%%
nst = 3;        %No. of stations 
d = [dc dc];    %diameter of shafts in m.
L = [le lc];    %length of shafts in m
Ipolar = [Ip1, Ip2, Ip3];     %polar mass moment of inertia of disc in kg-m^2

syms w
T = eye(2);
for i =1:nst    % loop for finding field and point matrix
    if(i~=1)
        J = pi*d(i-1)^4/32;
        l =L(i-1);
        kt=G*J/l;
        F(:,:,i) = [1  1/kt; %field matrix
                    0  1];
       T = F(:,:,i)*T;
    end
    Ip=Ipolar(i);
    P(:,:,i)  = [1         0;  %point matrix
                -w^2*Ip 1];
    T = P(:,:,i)*T;
end
f = T(2,1);     %frequency equation for free-free BC
n = polynomialDegree(f)/2;
delta = 1000;     %interval width
x1 = 0;         %initial value of w
k = 1;          
e = 0.0001;
m = n;          %total number of critical speeds
wnf = zeros(1,m);
if(subs(f,w,x1)==0)
    wnf(1) = 0; %1st critical speed if it is zero
    k = 2;
end

for i = k:m     %loop for other critical speeds
    while(wnf(i)==0)    
        x2 = x1+delta;   %taking arbitrary interval
        fx1 = subs(f,w,x1);
        fx2 = subs(f,w,x2);
        
        if(fx1*fx2<0)   %checking whether root is available in interval
            while(abs(x1-x2)>e) %bisection method
                xm = (x1+x2)/2;
                fxm = subs(f,w,xm);
                
                if (fx1*fxm<0)
                    x2 = xm;
                    fx2 = fxm;
                elseif(fx1*fxm>0)
                    x1 = xm;
                    fx1 = fxm;
                else
                    x1 = xm;
                    x2 = xm;
                    fx1 = fxm;
                    fx2 = fxm;
                end
            end
            wnf(i) = xm;    %critical speed or root of function
            x1 =x2;         % going to next interval
        else
            x1 = x2;
            if x2>15000         %loop to in case of less number of roots
                break;
            end
        end
    end
end

disp('using TMM method');
disp('Natural frequencies are');
disp(wnf);

 S(:,1) = [1; 0];
 S = sym(S);
 
for i=1:nst         %loop for finding corresponding state vector
   if i~=1
        S(:,2*i-1) = F(:,:,i)*S(:,2*(i-1));
   end
        S(:,2*i) = P(:,:,i)*S(:,2*i-1);
end

phiz =[];
Tz = [];
nzw = [];   %nonzero natural frequencies
%Substituting obtained natural frequencies in symbolic form of corresponding state vectors obtained above, to give numerical values of state variables
%Collecting displacements in phiz, torques in Tz
for i = 1:length(wnf)          
    if wnf(i)~=0               
        s = subs(S,w,wnf(i));
        phiz = [phiz;s(1,:)];
        Tz = [Tz;s(2,:)];
        nzw = [nzw,wnf(i)];
    end
end

phiz(1,:) = phiz(1,:)/phiz(1,6);
phiz(2,:) = phiz(2,:)/phiz(2,6);
L = [1, 0.4];
[r,c]  = size(phiz);   %r = number of non zero critical speed
x = [0,0];
for i=1:nst-1
    x = [x, x(2*i)+L(i),x(2*i)+L(i)]; %station coordinates
end

h = figure(1);
set(gcf, 'Position', get(0,'Screensize'));
for i = 1:r
    if(i==1)
         plot(x,phiz(i,:),'k-', 'Displayname', ['\omega=',num2str(nzw(i))],'LineWidth', 2);
         hold on;
    elseif(i==2)
         plot(x,phiz(i,:),'k-.','Displayname',['\omega=',num2str(nzw(i))], 'LineWidth', 2);
         hold on;
    elseif(i==3)
         plot(x,phiz(i,:),'k--','Displayname',['\omega=',num2str(nzw(i))], 'LineWidth', 2);
         hold on;
    else
         plot(x,phiz(i,:),'k:','Displayname',['\omega=',num2str(nzw(i))], 'LineWidth', 2);
         hold on;
    end
end
grid on
title ('Relative angular displacement with shaft length','fontsize',14);
ylabel('Relative Phiz');
xlabel('Shaft length');
legend('show');

sympref('floatingPointOutput',true);
fprintf (Output, '\n\nUsing Transfer Matrix Method\n');
fprintf (Output,'Number of stations = %d\n',nst);
fprintf (Output,'Point matrices :\n');
for i = 1:nst
    fprintf (Output,'[P]%d :\n',i);
    for a=1:2
        for b = 1:2
            fprintf (Output,'%s\t\t',char(P(a,b,i)));
        end
        fprintf (Output,'\n');
    end
    fprintf (Output,'\n');
end

fprintf (Output,'Field matrices :\n');
for i = 2:nst
    fprintf (Output,'[F]%d :\n',i-1);
    for a=1:2
        for b = 1:2
            fprintf (Output,'%4f\t',F(a,b,i));
        end
        fprintf (Output,'\n');
    end
    fprintf (Output,'\n');
end
fprintf (Output,'Natural frequencies:\n');
fprintf (Output,'%.2f\n',wnf);
fclose(Output);