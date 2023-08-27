clc
clear all

%Given system parametrs
G = 80*10^9;        %Rigidity modulus of the shaft in N/m^2.
d = [0.01, 0.01];  %Diameter of the shafts in m.
rho = 7800;        %density of the shaft in kg/m^3.

nele = 2;           %No of elements into which the system is discretised
connect = [1, 1, 2; %1st column is element no
           2, 2, 3]; %2nd and 3rd columns are node numbers in sequence
coord = [1, 0;      %1st column is node number
         2, 0.5;    %2nd column is x coordinate in m.
         3, 1];
Ip = [0.02, 0, 0.08];   %Polar mass moment of inertia of the discs in kg-m^2.

M = zeros(nele+1, nele+1);    %Initializing Global mass and stiffness matrices
K = zeros(nele+1, nele+1);

%Without considering mass of shaft
for i=1:nele        %Loop for finding elemental & global mass and stiffness matrices 
    nd1 = connect(i,2);
    nd2 = connect(i,3);
    x1 = coord(nd1,2);
    x2 = coord(nd2,2);
    l = x2-x1;
    J = pi*d(i)^4/32;
    kele(:,:,i) = G*J/l*[1, -1;-1, 1];
    mele(:,:,i) = [Ip(nd1), 0; 0, Ip(nd2)];
    vec = [nd1, nd2];       %Global DOF vector for assembly
     for ii = 1:2
        for jj = 1:2
            K(vec(ii), vec(jj)) = K(vec(ii), vec(jj))+kele(ii,jj,i);
            M(vec(ii), vec(jj)) = M(vec(ii), vec(jj))+mele(ii,jj,i);
        end
    end
end

%Free-free boundary condition
Kred = K;      %Reduced stiffness and mass matrices
Mred = M;

D = Kred\Mred;    
[eig_vec,eig_val] = eig(D);

for i=1:nele+1
    wnf(i) = sqrt(1/eig_val(i,i));      %Natural frequency
    mode(:,i) = eig_vec(:,i)/max(abs(eig_vec(:,i)));   %Normalised eigen vector
end

for p=1:length(wnf)-1           %Arranging natural frequencies in ascending order
  for i =1:length(wnf)-1        %And accordingly arranging mode shapes as well 
    if wnf(i)>wnf(i+1)
        temp_wnf = wnf(i);
        wnf(i) = wnf(i+1);
        wnf(i+1) = temp_wnf;
        temp_mode(:,1) = mode(:,i);
        mode(:,i) = mode(:,i+1);
        mode(:,i+1) = temp_mode(:,1);
    end
  end
end

wnf = wnf(:, size(wnf, 2)-2);
mode = mode(:, size(mode, 2)-2);

%Post processing
disp('Solution is printed to a text file "Assign4_Q5_224103438.txt"');
Output = fopen('Assign4_Q5_224103438.txt', 'w');
fprintf (Output, 'Given system is two disc rotor shaft system\n');
fprintf (Output, 'Using FEM (Neglecting shaft mass)\n');
fprintf (Output, '-----------------------------------\n');
fprintf (Output, 'No.of elements = %d\n', nele);
fprintf (Output, 'Elemental mass and stiffness matrices are shown below:\n');
for i=1:nele
    fprintf (Output,'Element-%d\n',i);
    fprintf (Output,'----------\n');
    fprintf (Output,'Mass Matrix [M]%d\n',i);
    for ii=1:2
        for jj=1:2
            fprintf(Output,'%.4f\t', mele(ii,jj,i));
        end
        fprintf(Output,'\n');
    end
    fprintf(Output,'\n');
    fprintf(Output,'Stiffness matrix [K]%d\n',i);
    for ii=1:2
        for jj=1:2
            fprintf(Output,'%.4f\t', kele(ii,jj,i));
        end
        fprintf(Output,'\n');
    end
    fprintf(Output,'\n');
end

fprintf (Output, 'Global Mass matrix is shown below:\n');
for i=1:nele+1
    for j= 1:nele+1
        fprintf (Output,'%.4f \t', M(i,j));
    end
    fprintf (Output,'\n');
end
fprintf (Output,'\n');
fprintf (Output, 'Global Stiffness matrix is shown below:\n');
for i=1:nele+1
    for j= 1:nele+1
        fprintf (Output,'%1.2e \t', K(i,j));
    end
    fprintf (Output,'\n');
end
fprintf (Output,'\n');
fprintf (Output, 'Free-Free Boundary condition\n\n');

fprintf (Output,'Obtained Natural frequencies(in rad/s) are:\n');
fprintf (Output,'%.3f \n', wnf);
fprintf (Output,'\n');

fprintf (Output,'Normalised eigen vector matrix is:\n');
for i=1:nele+1
    for j= 1:length(wnf)
        fprintf(Output,'%.4f \t', mode(i,j));
    end
    fprintf(Output,'\n');
end
fprintf(Output,'\n');

for i = 1:length(wnf)       %This loop is for finding displacements at various intermediate positions of shaft to plot smooth mode shapes
    x_z = [];               %For a particular frequency, x values of various interval points on the shaft
    phi_z=[];               %For a particular frequency, displacements at all x_z over the entire length of shaft  
    for j = 1:nele
        nd1 = connect(j,2);  
        nd2 = connect(j,3);
        x1 = coord(nd1,2);
        x2 = coord(nd2,2);
        l = x2-x1;
        x = x1:l/10:x2;     %Dividing each element into smaller elements for smooth mode shape
        phi1 = mode(j,i);   %Nodal displacement values of jth element for ith natural frequency
        phi2 = mode(j+1,i);
        for k=1:length(x)
            z = x(k)-x1;    
            N1 = 1-(z/l);   %Shape fucntions
            N2 = z/l;
            phi(k)=[N1 N2]*[phi1; phi2];    %Displacement values at various interval points in jth element for ith natural frequency
        end
        x_z =[x_z,x];           %Collecting x coordinate values for all elements at ith frequency
        phi_z = [phi_z,phi];    %Collecting displacement values from all elements at ith frequency
        
    end       
    mode_1(:,i)=phi_z/max(abs(phi_z));  %Collecting Normalised mode shapes for all natural frequencies
end

h = figure(1);
set(gcf, 'Position', get(0,'Screensize'));

if(size(mode_1,2) < 4)
    nmode=size(mode_1,2);
else
    nmode=4;
end

for i = 1: nmode    %length(wnf)
    if(i == 1)
        plot(x_z, mode_1(:,i), '-k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 2)
         plot(x_z, mode_1(:,i), ':k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 3)
         plot(x_z, mode_1(:,i),'-.k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 4)
         plot(x_z, mode_1(:,i), '--k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    else
         plot(x_z, mode_1(:,i), 'DisplayName',['\omega=',num2str(wnf(i))]);
    end
    hold on;
end
grid on;
xlabel('Shaft length(m)','fontsize',12);
ylabel('Relative amplitude','fontsize',12);
title('Mode shapes','fontsize',16);
ylim([-1.2 1.2]);
legend('show');
text(0.55,1.15,'(Without considering shaft mass)','FontSize',8);

%Considering mass of shaft
for i=1:nele        %Loop for finding elemental & global mass and stiffness matrices 
    nd1 = connect(i,2);
    nd2 = connect(i,3);
    x1 = coord(nd1,2);
    x2 = coord(nd2,2);
    l = x2-x1;
    J = pi*d(i)^4/32;
    kele(:,:,i) = G*J/l*[1, -1;-1, 1];
    mele(:,:,i) = [Ip(nd1), 0; 0, Ip(nd2)] + (rho*J*l/6)*[2 1; 1 2];
    vec = [nd1, nd2];       %Global DOF vector for assembly
     for ii = 1:2
        for jj = 1:2
            K(vec(ii), vec(jj)) = K(vec(ii), vec(jj))+kele(ii,jj,i);
            M(vec(ii), vec(jj)) = M(vec(ii), vec(jj))+mele(ii,jj,i);
        end
    end
end

%Free-free boundary condition
Kred = K;      %Reduced stiffness and mass matrices
Mred = M;

D = Mred\Kred;    
[eig_vec,eig_val] = eig(D);

for i=1:nele+1
    wnf(i) = sqrt(eig_val(i,i));      %Natural frequency
    mode(:,i) = eig_vec(:,i)/max(abs(eig_vec(:,i)));   %Normalised eigen vector
end

for p=1:length(wnf)-1           %Arranging natural frequencies in ascending order
  for i =1:length(wnf)-1        %And accordingly arranging mode shapes as well 
    if wnf(i)>wnf(i+1)
        temp_wnf = wnf(i);
        wnf(i) = wnf(i+1);
        wnf(i+1) = temp_wnf;
        temp_mode(:,1) = mode(:,i);
        mode(:,i) = mode(:,i+1);
        mode(:,i+1) = temp_mode(:,1);
    end
  end
end

%Post processing
fprintf (Output, 'Using FEM (Considering shaft mass)\n');
fprintf (Output, '-----------------------------------\n');
fprintf (Output, 'No.of elements = %d\n', nele);
fprintf (Output, 'Elemental mass and stiffness matrices are shown below:\n');
for i=1:nele
    fprintf (Output,'Element-%d\n',i);
    fprintf (Output,'----------\n');
    fprintf (Output,'Mass Matrix [M]%d\n',i);
    for ii=1:2
        for jj=1:2
            fprintf(Output,'%.4f\t', mele(ii,jj,i));
        end
        fprintf(Output,'\n');
    end
    fprintf(Output,'\n');
    fprintf(Output,'Stiffness matrix [K]%d\n',i);
    for ii=1:2
        for jj=1:2
            fprintf(Output,'%.4f\t', kele(ii,jj,i));
        end
        fprintf(Output,'\n');
    end
    fprintf(Output,'\n');
end

fprintf (Output, 'Global Mass matrix is shown below:\n');
for i=1:nele+1
    for j= 1:nele+1
        fprintf (Output,'%.4f \t', M(i,j));
    end
    fprintf (Output,'\n');
end
fprintf (Output,'\n');
fprintf (Output, 'Global Stiffness matrix is shown below:\n');
for i=1:nele+1
    for j= 1:nele+1
        fprintf (Output,'%1.2e \t', K(i,j));
    end
    fprintf (Output,'\n');
end
fprintf (Output,'\n');
fprintf (Output, 'Free-Free Boundary condition\n\n');

fprintf (Output,'Obtained Natural frequencies(in rad/s) are:\n');
fprintf (Output,'%.3f \n', wnf);
fprintf (Output,'\n');

fprintf (Output,'Normalised eigen vector matrix is:\n');
for i=1:nele+1
    for j= 1:length(wnf)
        fprintf(Output,'%.4f \t', mode(i,j));
    end
    fprintf(Output,'\n');
end
fclose(Output);

for i = 1:length(wnf)       %This loop is for finding displacements at various intermediate positions of shaft to plot smooth mode shapes
    x_z = [];               %For a particular frequency, x values of various interval points on the shaft
    phi_z=[];               %For a particular frequency, displacements at all x_z over the entire length of shaft  
    for j = 1:nele
        nd1 = connect(j,2);  
        nd2 = connect(j,3);
        x1 = coord(nd1,2);
        x2 = coord(nd2,2);
        l = x2-x1;
        x = x1:l/10:x2;     %Dividing each element into smaller elements for smooth mode shape
        phi1 = mode(j,i);   %Nodal displacement values of jth element for ith natural frequency
        phi2 = mode(j+1,i);
        for k=1:length(x)
            z = x(k)-x1;    
            N1 = 1-(z/l);   %Shape fucntions
            N2 = z/l;
            phi(k)=[N1 N2]*[phi1; phi2];    %Displacement values at various interval points in jth element for ith natural frequency
        end
        x_z =[x_z,x];           %Collecting x coordinate values for all elements at ith frequency
        phi_z = [phi_z,phi];    %Collecting displacement values from all elements at ith frequency
        
    end       
    mode_1(:,i)=phi_z/max(abs(phi_z));  %Collecting Normalised mode shapes for all natural frequencies
end

h = figure(2);
set(gcf, 'Position', get(0,'Screensize'));

if(size(mode_1,2) < 4)
    nmode=size(mode_1,2);
else
    nmode=4;
end

for i = 1: nmode    %length(wnf)
    if(i == 1)
        plot(x_z, mode_1(:,i), '-k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 2)
         plot(x_z, mode_1(:,i), ':k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 3)
         plot(x_z, mode_1(:,i),'-.k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    elseif (i == 4)
         plot(x_z, mode_1(:,i), '--k', 'LineWidth', 2, 'DisplayName',['\omega=',num2str(wnf(i))]);
    else
         plot(x_z, mode_1(:,i), 'DisplayName',['\omega=',num2str(wnf(i))]);
    end
    hold on;
end
grid on;
xlabel('Shaft length(m)','fontsize',12);
ylabel('Relative amplitude','fontsize',12);
title('Mode shapes','fontsize',16);
ylim([-1.2 1.2]);
legend('Location','south');
text(0.35,1.1,'(Considering shaft mass)','FontSize',8);