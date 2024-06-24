% This Code is used to validate the propsed model with the model of Nordt
% et al. A. A. Nordt, ‘Computing the mechanical properties of alpine skis’, Sports Eng., vol. 2, no. 2, p. 65, May 1999, doi: 10.1046/j.1460-2687.1999.00026.x.
% Last Update 10-05-2024
clear all
clc
close all
clear length
% Define data
Validationx = [-0.738095043,-0.571428528,-0.36507948,-0.190475974,-0.015873074,0.190475974,0.373015866,0.579364914,0.761904805,1];
Validationy = [40.67824007,118.6440371,315.2543885,505.0849255,559.3221479,464.406944,305.0849255,115.2545178,30.50877703,6.779943749];
ValiNordtx  =  [-0.992063615,-0.833333485,-0.642857208,-0.515873225,-0.396825627,-0.285714415,-0.222222121,-0.079365368,0.015873074,0.111110909,0.182539589,0.214285736,0.230158507,0.301587186,0.35714249,0.373015866,0.36507948,0.396825022,0.436507554,0.484127078,0.571428528,0.674602749,0.761904805,0.904761559,1.007936385];
ValiNordty  =  [-3.389519284,27.11886982,64.40707335,152.5423334,247.4577959,379.6612032,447.4577959,538.9832218,522.0340736,464.406944,406.7798144,518.6441664,623.7289626,583.0509812,535.5933146,454.237481,379.6612032,301.6951476,227.1188698,179.6612032,115.2545178,57.62738823,33.89829632,23.72909191,0];

% Global parameters
leng = 1.886; % Length of the Board in meters
Lr = 1.645;
YACP = -0.8285;
YFCP = 0.8285;
Ytail = -0.867;
Ytip = 1.019;

thickness = 2e-3; % Thickness of the Board in meters

% Array for the SKi K2 ChK 204 Dimensions
YL = [Ytail,YACP,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,YFCP,1,Ytip]; % Keypoints along the length
tL = 1e-3 * [66,67,73,113,158,172,151,114,71,65,65,65]; % Total Thickness along X direction
b = 1e-4 * [832,867,834,721,659,649,690,782,925,974,220,0]; % Keypoints along the width
kL = 1e-3 * [6,0,.6,2.9,4.2,4.70,4.2,2.90,.6,0,24.5,30.1]; % Keypoints along the z direction

num_elements = length(YL); % Number of elements


HPly = 1e-4 * [0.205,0.205,0.255,5.550,5.550,0.255,0.205,0.205]; % Thickness of composite layers
TotalThickness = sum(HPly);

% Material properties
% ABS
E_abs = 0.055e6;
G_abs = 0.02e6;
nu_abs = 0.4;
h_abs = 8.89e-4;
h_left = 3.89e-2;
% Steel
h_steel = 0.305e-2;
E_steel = 200e9; % Young's modulus in Pa of steel
G_steel = 80e6; % Shear modulus in Pa
nu_steel = 0.3;
% CFRP
E_Grfp = 40.510e9; % Young's modulus in Pa
E2_Grfp = 7.845e9; % Young's modulus in Pa (perpendicular to fibers)
nu_Grfp = 0.283; % Poisson's ratio
G_Grfp = 3.220e9; % Shear modulus in Pa
t_Gfrp = thickness; % Thickness of each layer in meters
h_Grfp = 8.64e-4;
% Wood
E_wood = 14e9; % Young's modulus in Pa for wood
G_wood = 0.710e9;
nu_wood = 0.3; % Poisson's ratio for wood
h_wood = thickness; % Thickness of wood layer in meters

b2 = 8.89e-4;
b4 = 3e-3;
b1 = b - 2 * b4;
h1 = TotalThickness;
h2 = TotalThickness - 2.03e-3;
h4 = 1.5e-3;
r1 = tL / 2; % The distance between the global centroid and the portion center
r4 = h4 / 2;
r2 = r4 + h2 / 2;


E_Materail = [E_Grfp,E_Grfp,E_Grfp,E_wood,E_wood,E_Grfp,E_Grfp,E_Grfp];
G_Materail = [G_Grfp,G_Grfp,G_Grfp,G_wood,G_wood,G_Grfp,G_Grfp,G_Grfp];
nu_Materail = [nu_Grfp,nu_Grfp,nu_Grfp,nu_wood,nu_wood,nu_Grfp,nu_Grfp,nu_Grfp];

el = leng / num_elements;
P = 1;

r1 = tL / 2; % The distance between the global centroid and the portion center
r4 = h4 / 2;
r2 = r4 + h2 / 2;


angles = [90,90,0,0,0,90,90]; % Angle of each CFRP layer in degrees

% Function to calculate ABD matrices
function [A, B, D, inverse_Stifness_matrix] = calculate_ABD_matrices(E_mat, nu_mat, G_mat, h_mat, angles)
    % Initialize matrices A, B, and D
    A = zeros(3, 3);
    B = zeros(3, 3);
    D = zeros(3, 3);
    
    % Calculate A, B, and D matrix elements
    for i = 1:length(h_mat)-1
        if i == 4 || i == 5
            Q_mat = [E_mat(i)/(1 - nu_mat(i)^2), nu_mat(i)*E_mat(i)/(1 - nu_mat(i)^2), 0;
                     nu_mat(i)*E_mat(i)/(1 - nu_mat(i)^2), E_mat(i)/(1 - nu_mat(i)^2), 0;
                     0, 0, G_mat(i)];
        else
            Q_mat = [E_mat(i)/(1 - nu_mat(i)^2), nu_mat(i)*E_mat(i)/(1 - nu_mat(i)^2), 0;
                     nu_mat(i)*E_mat(i)/(1 - nu_mat(i)^2), E_mat(i)/(1 - nu_mat(i)^2), 0;
                     0, 0, G_mat(i)];
            angle = deg2rad(angles(i));
            cosA = cos(angle);
            sinA = sin(angle);
            cs = cosA * sinA;
            cc = cosA^2;
            ss = sinA^2;
            T = [cc, ss, cs;
                 ss, cc, -cs;
                 -2*cs, 2*cs, cc-ss];
            Q_mat = T' * Q_mat * T;
        end
        
        A = A + Q_mat * (h_mat(i+1) - h_mat(i));
        B = B + 0.5 * Q_mat * (h_mat(i+1)^2 - h_mat(i)^2);
        D = D + (1/3) * Q_mat * (h_mat(i+1)^3 - h_mat(i)^3);
    end
    
    % Calculate the compliance matrix
    A_B = [A, B; B, D];
    inverse_Stifness_matrix = inv(A_B);
end

% Function to calculate the bending stiffness of a composite laminate

   
function EI = calculate_bending_stiffness(inSM1, E_steel, E_abs, b1, b2, b4, r1, r2, r4, t4, t2)
    beta1 = inSM1(1:3, 4:6);
    delta1 = inSM1(4:6, 4:6);
    alfa1 = inSM1(1:3, 1:3);
    
    EI1 = b1 * (2 * beta1(3, 3) * r1 - delta1(3, 3) * r1^2 - alfa1(3, 3)) / (beta1(3, 3)^2 - delta1(3, 3) * alfa1(3, 3));
    EI2 = 2 * E_abs * b2 * t2 * (t2^2 / 12 + r2^2);
    EI3 = 2 * E_steel * b4 * t4 * (t4^2 / 12 + r4^2);

    G1 = (b1 * beta1(3, 3) - r1 * delta1(3, 3) * alfa1(3, 3)) / (beta1(3, 3)^2 - delta1(3, 3) * alfa1(3, 3));
    G4 = (r4 * b4 * E_abs * t4);
    G2 = (r4 * b4 * E_steel * t4);

    G3 = -(b1 * beta1(3, 3) - b1 * alfa1(3, 3)) / (beta1(3, 3)^2 - delta1(3, 3) * alfa1(3, 3));
    G5 = (r2 * b2 * E_abs * t2);
    G6 = (r4 * b4 * E_steel * t4);

    RC = (G1 + 2 * G2 + 2 * G4) / (G3 + 2 * G5 + 2 * G6);
    EI = EI1 + EI2 + EI3;
end

% Initialize variables
element_length = ones(1, num_elements) * (leng / num_elements);
BendingStiff = zeros(length(b), 1);

% Calculate bending stiffness
for j = 1:length(b)
    [A1, B1, D1, inverseM1] = calculate_ABD_matrices(E_Materail, nu_Materail, G_Materail, HPly, angles);
    EI = calculate_bending_stiffness(inverseM1, E_steel, E_abs, b1(j), b2, b4, r1(j), r2, r4, h4, h2);
    BendingStiff(j) = EI;
end

% Plot results
figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
plot(YL, BendingStiff, 'k-', 'DisplayName', 'Present Model');
hold on;
plot(Validationx, Validationy, 'k*', 'DisplayName', 'Nordt et.all (Data)');
plot(ValiNordtx, ValiNordty, 'k--', 'DisplayName', 'Nordt et.all (Model)');
xlim([-0.850, 1.0]);
legend;
xlabel('X-coordinate (m)');
ylabel('Bending Stiffness EI (Nm^2)');
grid on;
%exportgraphics('Figure 1','Validation.png','Resolution',300)
