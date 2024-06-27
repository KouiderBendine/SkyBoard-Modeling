%% Author : Kouider Bendine
%  Last Update: 05-06-2024
%  MATLAB Code for free style Composite Ski mechanical behaviour  
% The present code is relied on the modelling approach proposed by:
% {A. A. Nordt, ‘Computing the mechanical properties of alpine skis’, 
% Sports Eng., vol. 2, no. 2, p. 65, May 1999, doi:
% 10.1046/j.1460-2687.1999.00026.x.}
% The code is able to perform three type of analysis : Static, Modal and
% contact snow-Ski (in the case of contact the snow is assumed to have a pure elastic behaviour)
% Required Functions are : Geometry_Xcoord, Geometry_Ycoord and BendingStifness
%% Email : kouider.bendine@list.lu

clear all
clc
close all
tic

% Define Geometrical parameters

Lc = 1100e-3;   % Contact Length
Ls = 175e-3;    % Shovel length
Lh = 175e-3;    % Heel Length
Hs = 55e-3;     % Tip Height
Hf = 20e-3;     % Free Bottom Camber
Hh = 55e-3;     % Heel Height
bm = 240e-3;    % width at the center of the ski
Wh = 290e-3;    % width at the heel
Ws = 290e-3;    % width at the shovel

prompt = "Please type the required analysis? -  1 : Static  - 2 : Modal - 3 : for Contact  : ";
Analysis = input(prompt);
       if isempty(Analysis)
          Analysis = 1;
       end
disp('The x and z coordinate of the ski  is undergoing')
[PtXcoo,PtZcoo,PtZcoou,Xf_filtered,Lc]=Geometry_XZcoord(Lc,Ls,Lh,Hs,Hf,Hh);
disp('The Y coordinate of the ski  is undergoing')
[PtXcoob,PtYcoob]=Geometry_Ycoord(Xf_filtered,PtZcoo,Lc,bm,Wh,Ws);
disp('The calculation of the bending stifness is undergoing')
[Mass,EI,RC,IPT,YL,kL,BendingStiff]=BendingStifness(PtXcoo,PtYcoob,PtZcoo,PtZcoou);
disp('The construction of the Finite element model is on running')
num_elements = length(YL) ;
BcModal = floor(length(YL) / 2);
n = 2;                                      % degrees of freedom per node
num_dofs = n * (num_elements + 1);           % Total Number of degrees of freedom
element_nodes = zeros(num_elements, 2);
el = 2*max(YL) / num_elements;
%el = length(YL) / num_elements;
% Generate element nodes using a for loop
for i = 1:num_elements
    node1 = 2 * i - 1;
    node2 = 2 * i;
    element_nodes(i, :) = [node1, node2];
end
%disp('Element Nodes:');
%disp(element_nodes);

function [KE, Force] = element_stiffness_matrix(EI_elements, el, P)
    KE = (EI_elements/el^3) * [12  ,  6*el, -12, 6*el;
                               6*el,  4*el^2, -6*el, 2*el^2;
                              -12  , -6*el, 12, -6*el;
                               6*el,  2*el^2, -6*el, 4*el^2];
    Force = transpose([P*el/2, P*el^2/12, P*el/2, -P*el^2/12]);
end

% Assemble global stiffness matrix
global_stiffness = zeros(num_dofs, num_dofs);
global_Force = zeros(num_dofs, 1);
KE_iter = zeros(4, 4, num_elements + 1);

for i = 1:(num_elements-1)
    x1 = YL(i);
    x2 = YL(i + 1);
    z1 = kL(i);
    z2 = kL(i + 1);

    lengt = sqrt((x2 - x1)^2 + (z2 - z1)^2);
    lx = (x2 - x1) / 2;
    mx = (z2 - z1) / 2;
    cosa = lx / lengt;
    sena = mx / lengt;

    T = [cosa, sena, 0, 0;
         -sena, cosa, 0, 0;
         0, 0, cosa, sena;
         0, 0, -sena, cosa];

    [KE_iter(:,:,i), FE] = element_stiffness_matrix( BendingStiff(i), el, 0);
    % KE_iter(:,:,i) = T' * KE_iter(:,:,i) * T;
end

global_mass = eye(num_dofs);
for itr = 1:(num_elements - 1)
    global_mass(2*itr, 2*itr) = Mass(itr, 1);
    global_mass(2*itr+1, 2*itr+1) = IPT(itr, 1);
end

for itr = 1:2:(2 * (num_elements + 1) - 2)
    global_stiffness(itr:itr+3, itr:itr+3) = global_stiffness(itr:itr+3, itr:itr+3) + KE_iter(:,:,i);
    global_Force(itr:itr+3) = global_Force(itr:itr+3) + FE;
end
disp('Entring the Analysis')
switch Analysis
    case 1
    disp('Static Deflection ');
   

    indices = [2*11, 2*11+1, 2*32, 2*32+1];
    global_stiffness(indices, :) = [];
    global_stiffness(:, indices) = [];
    global_mass(indices, :) = [];
    global_mass(:, indices) = [];
    global_Force(indices) = [];

    %global_Force(2*16) = -200
    %global_Force(2*17) = -200;
    global_Force(2*18) = 60;
    %global_Force(2*19) = -200;
    global_Force(2*22) = -60;
    %global_Force(2*21) = -200;
    %global_Force(2*22) = -200
    %global_Force(2*23) = -200;
    [displacements, reactions] = solve_displacements(global_stiffness, global_Force);

    %Uz_Disp = [0; displacements(2:2:num_dofs-4); 0];
    %Ux_Disp = [displacements(1:2:num_dofs-4); displacements(num_dofs-4)];

    Uz_Disp = [displacements(2:2:end)];
    Ux_Disp = [displacements(1:2:end)];

    disp('The maximum Deflection is :')
    disp(max(abs(Ux_Disp)));

    figure('Color', 'w')
    set(gca, 'FontSize', 20)
    set(gca, 'FontName', 'Times New Roman')
    fontsize(20,"points")
    
    plot(YL(1:end-1), kL(1:end-1), 'k--*');
    hold on;
    plot(YL(1:end-1), kL(1:end-1)'-8*((Ux_Disp(1:end))),'-');
    ylabel('Displacements [m]');
    xlabel('X Coordinate [m]');
    grid on;
    %legend('Location', 'upper center');
    xlim([-1.1, 1.1]);
    ylim([-0.04, 0.22]);
    hold off;

    case 2
    disp('Modal Analysis');

    [global_mass, global_stiffness, frequencies, evecs] = modal_analysis(global_mass, global_stiffness, BcModal);

    figure('Color', 'w')
    set(gca, 'FontSize', 20)
    set(gca, 'FontName', 'Times New Roman')
    fontsize(20,"points")
    subplot(2, 2, 1);
    plot(YL(2:end-1), kL(2:end-1)' + 0.1 * evecs(1:2:end-4, 1), 'r--.');
    hold on;
    plot(YL(2:end-1), kL(2:end-1)' + 0.1 * evecs(1:2:end-4, 2), 'b--.');
    plot(YL(2:end-1), kL(2:end-1), 'k--*');
    xlabel('X Coordinate [m]');
    legend('Afterbody','Forebody','Undeformed');
    grid on;
    %title('1st Bending Mode');
    hold off;

    subplot(2, 2, 2);
    plot(YL(2:end-1), kL(2:end-1)' + 0.08 * evecs(1:2:end-4, 5), 'r--.');
    hold on;
    plot(YL(2:end-1), kL(2:end-1)' + 0.05 * evecs(1:2:end-4, 6), 'b--.');
    plot(YL(2:end-1), kL(2:end-1), 'k--*');
    xlabel('X Coordinate [m]');
    legend('Afterbody','Forebody','Undeformed');
    grid on;
    %title('2nd Bending Mode');
    hold off;

    subplot(2, 2, 3);
    plot(YL(2:end-1), kL(2:end-1)' + 0.1 * evecs(1:2:end-4, 7), 'r--.');
    hold on;
    plot(YL(2:end-1), kL(2:end-1)' - 0.1 * evecs(1:2:end-4, 8), 'b--.');
    plot(YL(2:end-1), kL(2:end-1), 'k--*');
    xlabel('X Coordinate [m]');
    legend('Afterbody','Forebody','Undeformed');
    grid on;
    %title('3rd Bending Mode');
    hold off;

    subplot(2, 2, 4);
    plot(YL(2:end-1), kL(2:end-1)' + 0.1 * evecs(1:2:end-4, 12), 'r--.');
    hold on;
    plot(YL(2:end-1), kL(2:end-1)' + 0.1 * evecs(1:2:end-4, 13), 'b--.');
    plot(YL(2:end-1), kL(2:end-1), 'k--*');
    xlabel('X Coordinate [m]');
    legend('Afterbody','Forebody','Undeformed');
    grid on;
    %title('4th Bending Mode');
    hold off;

    case 3  
    disp('Contact Ski-Snow');

clear length

% Initialize NewYL array with zeros
NewYL = zeros(2 * length(PtZcoo(1:end-1)), 1);

% Fill the array of zeros with the elements of the original vector
NewYL(1:2:end) = PtZcoo(1:end-1);

% Assign to c
c = PtZcoo(1:end-1);

% Reduce the size of global_stiffness
global_stiffness = global_stiffness(1:2:end-2, 1:2:end-2);

U = 1;
fc = 100;
Ks1 = 10e6;     % (thesis) brenane MODELING THE MECHANICAL CHARACTERISTICS (page 32 and 64)
toux = 15e3;    % (thesis) brenane MODELING THE MECHANICAL CHARACTERISTICS

% Function to calculate displacements


num_load_steps = 30;
Ks = zeros(size(global_stiffness));
Ks(9, 9) = Ks1;
Ks(33, 33) = Ks1;

Fz = zeros(size(global_stiffness, 1), num_load_steps);
Fz(17:19, 1) = -610;
Fz(22:24, 1) = -610;

d = zeros(size(global_stiffness));
P = zeros(size(global_stiffness, 1), num_load_steps);
Pload = zeros(size(global_stiffness, 1), num_load_steps);
dz = zeros(size(global_stiffness, 1), num_load_steps);
DeltZ = zeros(size(global_stiffness, 1), num_load_steps);

for t = 2:num_load_steps
    dz(:, t) = calculate_displacements(Fz(:, t-1), global_stiffness, Ks);
    dz(:, t) = dz(:, t) + dz(:, t-1);
    Fz(17:19, t) = Fz(17:19, t-1) - 50;
    Fz(22:24, t) = Fz(22:24, t-1) - 50;
    DeltZ(:, t) = PtZcoo(1:end)' + dz(1:end, t);
    
    for it = 1:numel(DeltZ(:, t))
        if DeltZ(it, t) < 0
            DeltZ(it, t) = DeltZ(it, t);
        else
            DeltZ(it, t) = 0;
        end
    end
    
    Pload(:, t) = (abs(DeltZ(:, t)))' * Ks./0.5;
    
    for i = 1:numel(dz(:, t))
        if abs(dz(i, t)) >= PtZcoo(i) || abs(dz(i, t)) > PtZcoo(i)
            Ks(i, i) = Ks1;
            else
            
        end
    end
end

%% Plotting

figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
for t = 2:num_load_steps
plot(YL(5:end-4), PtZcoo(5:end-4) + dz(5:end-4, t-1)','LineWidth', 2,'color','k');
hold on;
end
plot(YL(5:end-4), PtZcoo(5:end-4),'LineWidth', 2,'color', 'r');
ylabel('Displacements [m]');
xlabel('X Coordinate [m]');
ylim([-0.01, 0.04]);
xlim([-0.8, 0.8]);
legend('Deformed Shape', 'Undeformed Shape');
grid on;

figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
fontsize(20,"points")
for t = 2:num_load_steps
    plot(YL(1:end-1), Pload(1:end-1, t)/1e3,'LineWidth', 2,'color', 'k');
    hold on;
end
ylabel('Pressure [KPa]');
xlabel('X Coordinate [m]');
xlim([-0.6, 0.6]);
grid on;

figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
fontsize(20,"points")
for t = 2:num_load_steps
    plot3( YL(8:end-7), t * ones(size(YL(8:end-7))), PtZcoo(8:end-7) + dz(8:end-7, t-1)','LineWidth', 2, 'color','k');
    hold on;
end
xlabel('X Coordinate');
ylabel('Load Steps');
zlabel('Displacements [m]');
xlim([-0.6, 0.6]);
zlim([-0.0055, 0.025]);
grid on;

figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
fontsize(20,"points")
for t = 2:num_load_steps
    plot3( YL(8:end-7), t * ones(size(YL(8:end-7))), Pload(8:end-7, t)/1e3,'LineWidth', 2, 'color','k');
    hold on;
end
xlabel('X Coordinate');
ylabel('Load Steps');
zlabel('Pressure [KPa]');
xlim([-0.65, 0.65]);
grid on;

figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
fontsize(20,"points")
    filename=['SResuls_APDL.txt'];
    filename1=['DResuls_APDL.txt'];
    A=importdata(filename);
    B=importdata(filename1);
    AnsysDis=A(:,1);
    AnsysCtr=A(:,2);
    AnsysDisD=B(:,1);
    AnsysCtrD=B(:,2);
 plot(YL(1:2:end),AnsysDis, 'k*--')  
 hold on
 plot(YL(1:2:end),AnsysDisD, 'k*--')
 plot(YL(8:end-7),   dz(8:end-7, t-1)','LineWidth', 2,'color','k')
xlabel('X Coordinate [m]');
ylabel('Displacement [m]');

xlim([-0.40, 0.40]);
legend('Ansys Static','Ansys Dynamic','Present Model')
grid on;
figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
fontsize(20,"points")
plot(YL(1:2:end),AnsysCtr/150, 'k*--') 
hold on
xlim([-0.40, 0.40]);
plot(YL(1:2:end),AnsysCtrD/300, 'b*--') 
hold on
plot(YL(1:end-1), Pload(1:end-1, t),'LineWidth', 2,'color', 'k');
legend('Ansys Static','Ansys Dynamic','Present Model')
xlabel('X Coordinate [m]');
ylabel('Contact Pressure [Pa]');
grid on;
disp([max(abs( dz(8:end-7, t-1))),max(abs(AnsysDis(8:14)))])
disp([max(abs( Pload(1:end-1, t))),max(abs(AnsysCtr/200))])
ErrorDisplacement =(0.0207-0.0184)/0.0207
ErrorPressure =(max(abs(AnsysCtr/200))-2.8532)/max(abs(AnsysCtr/200))
end

%% Required Functions

function [displacements, reactions] = solve_displacements(global_stiffness, global_Force)
        displacements = global_stiffness \ global_Force;
        reactions = global_stiffness * displacements;
end

function [global_mass, global_stiffness, frequencies, evecs] = modal_analysis(global_mass, global_stiffness, BcModal)
        
        BcModal = floor(BcModal);
        indices = [2*BcModal, 2*BcModal-1];
        global_stiffness(indices, :) = [];
        global_stiffness(:, indices) = [];
        global_mass(indices, :)      = [];
        global_mass(:, indices)      = [];
        [evecs, evals] = eig(global_stiffness, global_mass);
        frequencies = sqrt(diag(evals));
        disp('The first 10 Frequencies are :')
        disp(frequencies(2:5:35));
end

function d = calculate_displacements(Fz, global_stiffness, Ks)
         d = (global_stiffness + Ks) \ Fz;
end
toc