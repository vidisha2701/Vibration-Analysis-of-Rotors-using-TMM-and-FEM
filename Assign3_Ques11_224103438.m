% Birth Month:January

% A Motor connected to a Flywheel by a Two-stage Gear-Box
% Input Parameters
% Gear ratio
n1 = 3; 
n2 = 4;

% Modulus of Rigidity of Shafts in N/m^2
G = 8*10^10;   

Output = fopen('Ques_11_224103438_TMM.txt', 'w');
fprintf (Output,'A Motor connected to a Flywheel by a Two-stage Gear-Box\n');
fprintf (Output, 'Using TMM\n');

% No. of stations
nst = 6; 

% Polar Mass Moment of Inertia of the Motor and Flywheel in kg-m^2
IP = [0.01, 0, 0, 0, 0, 0.04]; 

% Diameter of the Shafts in m.
d = [0.015 0.012 0.01];    

% Length of the Shafts in m.
L = [0.2 0.4 0.2]; 

% Stiffness
KT = zeros(size(d));
for i=1:(size(d,2))
    KT(i) = G*pi()*d(i)^4/(32*L(i))
end

syms w
% Point Matrices and Field Matrices
P1 = [1 0;
      -w^2*IP(1) 1];
P2 = [1 0;
      -w^2*IP(2) 1];
P3 = [1 0;
      -w^2*IP(3) 1];
P4 = [1 0;
      -w^2*IP(4) 1];
P5 = [1 0;
      -w^2*IP(5) 1];
P6 = [1 0;
      -w^2*IP(6) 1];
F1 = [1 1/KT(1);
      0 1];
F2 = [1 1/KT(2);
      0 1];
F3 = [1 1/KT(3);
      0 1];

%  Gear Matrix

G1 = [1/n1 0;
     0   n1];
G2 = [1/n2 0;
     0   n2];

% Overall Transfer Matrix
T = P6*F3*P5*G2*P4*F2*P3*G1*P2*F1*P1;       

% Frequency Equation
f = T(2,1);       
n = polynomialDegree(f)/2;

delta = 50;     % Interval width
x1 = 0;         % Initial value of w
k = 1;          
e = 0.0001;     % Tolerance chosen
m = n;          % Total number of critical speeds
wnf = zeros(1,m);
if(subs(f,w,x1)==0)
    wnf(1) = 0;      %1st critical speed if it is zero
    k = 2;
end
for i = k:m     % loop for other critical speeds
    while(wnf(i)==0)    
        x2 = x1+delta;   % taking arbitrary interval
        fx1 = subs(f,w,x1);
        fx2 = subs(f,w,x2);
        
        if(fx1*fx2<0)   % checking whether root is available in interval
            while(abs(x1-x2)>e) % bisection method
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
            wnf(i) = xm    %critical speed or root of function
            x1 =x2;         % going to next interval
        else
            x1 = x2;
            if x2>1000         %loop to in case of less number of roots
                break;
            end
        end
    end
end

% Printing the results to output file
P(:,:,1) = P1;
P(:,:,2) = P2;
P(:,:,3) = P3;
P(:,:,4) = P4;
P(:,:,5) = P5;
P(:,:,6) = P6;
F(:,:,1) = F1;
F(:,:,2) = F2;
F(:,:,2) = F3;
GP(:,:,1) = G1;
GP(:,:,1) = G2;

sympref('floatingPointOutput',true);
fprintf (Output,'Number of stations = %d\n',nst);

fprintf (Output,'Point matrices :\n');
for i = 1:6
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
for i = 1:3
    if (i~=3)
        fprintf (Output,'[F]%d :\n',i);
        for a=1:2
            for b = 1:2
                fprintf (Output,'%4f\t', F(a,b,i));
            end
            fprintf (Output,'\n');
        end
        fprintf (Output,'\n');
    end
end

fprintf (Output,'Gear Ratio Transformation Matrix :\n\n');
fprintf (Output, '[G]:\n');
for a=1:2
    for b = 1:2
        fprintf (Output,'%4f\t', GP(a,b,1));
    end
    fprintf (Output,'\n');
end
fprintf (Output,'Free-Free boundary condition\n');
fprintf(Output,'t(2,1)=0\n\n');
fprintf(Output,'Natural Frequecies:\n');
fprintf(Output,'%f\n',wnf);
fclose(Output);