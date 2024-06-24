function [Mass,EI,RC,IPT,YL,kL,BendingStiff]=BendingStifness(PtXcoo,PtYcoob,PtZcoo,PtZcoou)
%% This function used to calculate the bending stifness, mass, polar moment of inertia based on the proposed geometry
%% Last Update 10-05-2024
YL = PtXcoo;                   % Coordinate along the length
tL = abs((PtZcoou - PtZcoo));  % Total Thickness along X direction
phi = 0.40;
delta = 0.05;

% Initialize variables
b = PtYcoob;               % Coordinate along the width
kL = PtZcoo;               % Kyepoints along the z direction
num_elements = length(YL); % Number of elements
b(1) = b(2);
b(end) = b(end-1);

% Material properties
% ABS
E_abs = 1.7e9;
G_abs = 0.6e9;
nu_abs = 0.49;
Roh_Abs = 1050;
h_abs = 4e-3;
h_left = 3.89e-3;
% Steel
h_steel = 2e-3;
E_steel = 200e9;  % Young's modulus in Pa of steel
G_steel = 80e9;   % Shear modulus in Pa
nu_steel = 0.3;
Roh_steel = 7900;
% CFRP
E_Grfp = 30.00e9;  % Young's modulus in Pa
E2_Grfp = 30.00e9; % Young's modulus in Pa (perpendicular to fibers)
nu_Grfp = 0.283;   % Poisson's ratio
G_Grfp = 6.7e9;    % Shear modulus in Pa
Roh_Grfp = 1660;
h_Grfp = 1e-3;
% Wood
E_wood = 14e9;     % Young's modulus in Pa for wood
G_wood = 5.5e9;
nu_wood = 0.3;     % Poisson's ratio for wood
Roh_wood = 500;
h_wood = 4e-3;

% Calculate HPlyh, HPlyc, HPlys arrays
HPlyh = [-(h_wood/3+h_Grfp/3), -(h_wood/3+2*h_Grfp/3), -(h_wood/3+h_Grfp/3), -h_wood/3, 0, h_wood/2, (h_wood/3+h_Grfp/3), (h_wood/3+2*h_Grfp/3), (h_wood/3+h_Grfp/2)];
HPlyc = [-(h_wood/2+h_Grfp/2), -(h_wood/2+2*h_Grfp/3), -(h_wood/2+h_Grfp/3), -h_wood/2, 0, h_wood/2, (h_wood/2+h_Grfp/3), (h_wood/2+2*h_Grfp/3), (h_wood/2+h_Grfp/2)];
HPlys = [-(h_wood/3+h_Grfp/3), -(h_wood/3+2*h_Grfp/3), -(h_wood/3+h_Grfp/3), -h_wood/3, 0, h_wood/2, (h_wood/3+h_Grfp/3), (h_wood/3+2*h_Grfp/3), (h_wood/3+h_Grfp/2)];

% Create HPly matrix
HPly = zeros(length(HPlyc), length(b));
HPly(:, 1:10)   = [HPlyh;HPlyh;HPlyh;HPlyh;HPlyh;HPlyh;HPlyh;HPlyh;HPlyh;HPlyh]';  % Fill first 4 columns with vector H
HPly(:, 11:20)  = [HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc]';  % Fill columns 4 to 6 with vector g
HPly(:, 21:30)  = [HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc]';
HPly(:, 31:end) = [HPlys;HPlys;HPlys;HPlys;HPlys;HPlys;HPlys;HPlys;HPlys;HPlys]';

% TotalThickness = sum(HPly);
TotalThickness = sum(sum(HPly));
thickness = TotalThickness; % 5e-3;

b4 = 4e-3;
b2 = thickness - b4;

b1 = b - 2*b4;
h1 = TotalThickness; % The Laminate total thickness
h2 = b4;             % The Abs total thickness
h4 = b4;             % The Steel total thickness
r1 = tL;             % The distance btw the global centroid and the portion center
r4 = h4;             % The distance btw the global centroid and the Steel portion center
r2 = r4 + h2/2;      % The distance btw the global centroid and the Abs portion center

% PolplaWood
E_Pwood = 8.5e9;     % Young's modulus in Pa for Binder
nu_Pwood = 0.3;      % Poisson's ratio for Binder
G_Pwood = 0.610e9;

E_Materail = [E_Grfp,E_Grfp,E_Grfp,E_wood,E_wood,E_Grfp,E_Grfp,E_Grfp];
G_Materail = [G_Grfp,G_Grfp,G_Grfp,G_wood,G_wood,G_Grfp,G_Grfp,G_Grfp];
nu_Materail = [nu_Grfp,nu_Grfp,nu_Grfp,nu_wood,nu_wood,nu_Grfp,nu_Grfp,nu_Grfp];

leng = max(abs(PtXcoo))*2;
el = leng / num_elements;

P = 1;
k = 30; % The curvature of the neutral longitudinal (d2w/dx2)

angles = [0,0,0,0,0,0,0,0]; % Angle of each CFRP layer in degrees



zjGrfp = abs((h_wood + h_Grfp/2) - 0.00406858); %% note that 0.00406858 is the center of mass
zjabs = abs(h_abs - 0.00406858);
zjwood = abs(h_wood - 0.00406858);
zjSteel = abs(h_steel - 0.00406858);
xjSteel = abs(b/2 - b4);
xjabs = abs(b/2 - b4);

clear length
% Calculate Mass and Polar Moment of Inertia
[Mass, IPT] = mass_fun(Roh_wood, Roh_steel, Roh_Abs, Roh_Grfp, b, h_steel, h_wood, h_Grfp, h_abs, el, zjGrfp, zjabs, zjwood, zjSteel, xjSteel, xjabs);

% Initialize Bending Stiffness and Neutral Axis
BendingStiff = zeros(length(b), 1);
Rcenter = zeros(length(b), 1);
for j = 1:length(b)
    [A1, B1, D1, inverseM1] = calculate_ABD_matrices(E_Materail, nu_Materail, G_Materail, HPly(:, j), angles);
    [EI, RC] = calculate_bending_stiffness(inverseM1, E_steel, E_abs, b1(j), b2, b4, r1(j), r2, r4, h4, h2);
    
    BendingStiff(j) = EI;
    Rcenter(j) = RC;
end

% Plot Bending Stiffness
figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
plot(YL, BendingStiff, 'k*');
xlabel('X-coordinate (m)');
ylabel('Bending Stiffnes EI (Nm^2)');
xlim([-0.650, 0.650]);
ylim([0.0, 800]);
grid on;

% Calculate total mass
total_mass = sum(Mass);
%disp(total_mass)


function [A, B, D, inverse_Stifness_matrix] = calculate_ABD_matrices(E_mat, nu_mat, G_mat, h_mat, angles)
    A = zeros(3, 3);
    B = zeros(3, 3);
    D = zeros(3, 3);

    for i = 1:length(angles)
        if i == 3 || i == 4
            Q_mat = [E_mat(i)/(1 - nu_mat(i)^2), nu_mat(i)*E_mat(i)/(1 - nu_mat(i)^2), 0;
                     nu_mat(i)*E_mat(i)/(1 - nu_mat(i)^2), E_mat(i)/(1 - nu_mat(i)^2), 0;
                     0, 0, G_mat(i)];
        else
            Q_mat = [E_mat(i)/(1 - nu_mat(i)^2), nu_mat(i)*E_mat(i)/(1 - nu_mat(i)^2), 0;
                     nu_mat(i)*E_mat(i)/(1 - nu_mat(i)^2), E_mat(i)/(1 - nu_mat(i)^2), 0;
                     0, 0, G_mat(i)];
            angle = deg2rad(angles(i));
            cos_val = cos(angle);
            sin_val = sin(angle);
            cs = cos_val * sin_val;
            cc = cos_val^2;
            ss = sin_val^2;
            T = [cc, ss, cs;
                 ss, cc, -cs;
                 -2*cs, 2*cs, cc-ss];
            Q_mat = T' * Q_mat * T;
        end

        A(1, 1) = A(1, 1) + Q_mat(1, 1) * (h_mat(i+1) - h_mat(i));
        A(2, 2) = A(2, 2) + Q_mat(2, 2) * (h_mat(i+1) - h_mat(i));
        A(1, 2) = A(1, 2) + Q_mat(1, 2) * (h_mat(i+1) - h_mat(i));
        A(2, 1) = A(2, 1) + Q_mat(2, 1) * (h_mat(i+1) - h_mat(i));
        A(3, 3) = A(3, 3) + Q_mat(3, 3) * (h_mat(i+1) - h_mat(i));

        D(1, 1) = D(1,1) + Q_mat(1, 1) * (h_mat(i+1)^3 - h_mat(i)^3) / 3;
        D(2, 2) = D(2, 2) + Q_mat(2, 2) * (h_mat(i+1)^3 - h_mat(i)^3) / 3;
        D(1, 2) = D(1, 2) + Q_mat(1, 2) * (h_mat(i+1)^3 - h_mat(i)^3) / 3;
        D(2, 1) = D(2, 1) + Q_mat(2, 1) * (h_mat(i+1)^3 - h_mat(i)^3) / 3;
        D(3, 3) = D(3, 3) + Q_mat(3, 3) * (h_mat(i+1)^3 - h_mat(i)^3) / 3;

        B(1, 1) = B(1, 1) + 0.5 * Q_mat(1, 1) * (h_mat(i+1)^2 - h_mat(i)^2);
        B(2, 2) = B(2, 2) + 0.5 * Q_mat(2, 2) * (h_mat(i+1)^2 - h_mat(i)^2);
        B(1, 2) = B(1, 2) + 0.5 * Q_mat(1, 2) * (h_mat(i+1)^2 - h_mat(i)^2);
        B(2, 1) = B(2, 1) + 0.5 * Q_mat(2, 1) * (h_mat(i+1)^2 - h_mat(i)^2);
        B(3, 3) = B(3, 3) + 0.5 * Q_mat(3, 3) * (h_mat(i+1)^2 - h_mat(i)^2);
    end

    A_B = [A, B; B, D];
    inverse_Stifness_matrix = inv(A_B);

end 

function [Mass, IPT] = mass_fun(Roh_wood, Roh_steel, Roh_abs, Roh_GRP, b, h_steel, h_wood, h_Grfp, h_abs, el, zjGrfp, zjabs, zjwood, zjSteel, xjSteel, xjabs)
    Mass = zeros(length(b), 1);
    IPT = zeros(length(b), 1);
    
    for i = 1:length(b)
        Mass(i) = (2 * (h_steel * b(i) * Roh_steel) + h_wood * b(i) * Roh_wood + h_Grfp * b(i) * Roh_GRP + 2 * (h_abs * b(i) * Roh_abs)) * el;

        TPhorGrfp = Roh_GRP * ((1 / 12) * b(i)^3 * h_Grfp + (1 / 3) * b(i) * (3 * zjGrfp^2 * h_Grfp + 1 / 4 * h_Grfp^3));  % Polar mass moment of inertia Grfp
        TPhorwood = Roh_wood * ((1 / 12) * b(i)^3 * h_wood + (1 / 3) * b(i) * (3 * zjwood^2 * h_wood + 1 / 4 * h_wood^3));  % Polar mass moment of inertia wood
        IPhor = TPhorwood + TPhorGrfp;

        TPhorSteel = Roh_steel * (xjSteel(i)^2 * b(i) * h_steel + (zjSteel^2 * h_steel * b(i)) + (1 / 12 * h_steel^3 * b(i)) + 1 / 12 * h_steel * b(i)^3);  % Polar mass moment of inertia steel
        TPhorabs = Roh_abs * (xjabs(i)^2 * b(i) * h_abs + (zjabs^2 * h_abs * b(i)) + (1 / 12 * h_abs^3 * b(i)) + 1 / 12 * h_abs * b(i)^3);  % Polar mass moment of inertia abs
        IPver = 2 * TPhorSteel + 2 * TPhorabs;
        IPT(i) = IPhor + IPver;
    end
end

function [EI, RC] = calculate_bending_stiffness(inSM1, E_steel, E_abs, b1, b2, b4, r1, r2, r4, t4, t2)
    beta1 = inSM1(1:3, 4:6);
    delta1 = inSM1(4:6, 4:6);
    alfa1 = inSM1(1:3, 1:3);

    EI1 = b1 * (2 * beta1(3, 3) * r1 - delta1(3, 3) * r1^2 - alfa1(3, 3)) / (beta1(3, 3)^2 - delta1(3, 3) * alfa1(3, 3));
    EI2 = 2 * E_abs * b2 * t2 * (t2^2 / 12 + r2^2);
    EI3 = 2 * E_steel * b4 * t4 * (t4^2 / 12 + r4^2);

    G1 = (b1 * beta1(3, 3) - r1 * delta1(3, 3) * alfa1(3, 3)) / (beta1(3, 3)^2 - delta1(3, 3) * alfa1(3, 3));
    G4 = (r4 * b4 * E_abs * t4);
    G2 = (r4 * b4 * E_steel * t4);

    G3 = (-b1 * alfa1(3, 3)) / (beta1(3, 3)^2 - delta1(3,3) * alfa1(3, 3));
    G5 = (b2 * E_abs * t2);
    G6 = (b4 * E_steel * t4);
    %disp([G1,G2,G3,G4,G5,G6])

    RC = (G1 + 2 * G2 + 2 * G4) ./ (G3 + 2 * G5 + 2 * G6);
    EI = EI1 + EI2 + EI3;
end
end