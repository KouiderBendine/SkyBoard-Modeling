% Define arrays
%clear all
%clc
close all
clear length
Brennanx = [-0.460975511, -0.44268283, -0.413414707, -0.395122026, -0.365853624, -0.347560943, -0.29634138, -0.230487895, ...
-0.186585293, -0.135366009, -0.084146445, -0.01097572, 0.058536524, 0.113414568, 0.168292611, 0.252439057, 0.2926829, ...
0.310975581, 0.354877904, 0.398780507, 0.428048629];

Brennany = [132.6086306, 160.8695686, 185.8696036, 204.3478326, 220.6521816, 236.9565306, 250.0000181, 242.3913136, 222.8261031, ...
205.4348141, 190.2174051, 190.2174051, 196.7391696, 215.2173986, 235.8695906, 230.4348076, 207.6087356, 191.3043866, 165.2174116, ...
140.2174181, 121.7390646];

BrennanEx = [-0.451359565, -0.422356709, -0.40060422, -0.371601364, -0.342598231, -0.306344453, -0.273715997, -0.233836896, ...
-0.197583118, -0.172205584, -0.135951806, -0.077945817, -0.034441117, 0.009063306, 0.059818651, 0.099697751, 0.132326207, ...
0.183081552, 0.237462219, 0.270090675, 0.306344453, 0.328096664, 0.378852009, 0.396978898, 0.418731109];

BrennanEy = [132.2492426, 159.3495952, 182.1138053, 207.0460867, 231.9783267, 244.9864645, 253.6585426, 241.7343991, 223.3062486, ...
211.3821051, 201.6260327, 192.9539133, 186.4498651, 184.2818352, 192.9539133, 208.130081, 222.222213, 234.1463566, 233.062321, ...
214.6341292, 197.289973, 178.8617812, 145.2574631, 133.3333195, 115.9890807];

% Global parameters
leng = 1.886;
Lr = 1.645;
YACP = -0.5588;
YFCP = 0.5588;
Ytail = -0.7226;
Ytip = 0.7365;
width = 8.89e-4;
thickness = 2e-3;
num_elements = 10;

% Array for the SKi K2 ChK 204 Dimenssions
YL = [Ytail, YACP, -0.5, -0.25, 0, 0.25, 0.5, YFCP, Ytip];
tLcore = 1e-3 * [1.9, 2.0, 2.8, 6.3, 5.8, 6.5, 3.5, 2.2, 1.9];
b = 1e-3 * [76.2, 271.9, 287.9, 242.9, 231.1, 236.8, 287.6, 272.6, 76.2];
kL = 1e-3 * [55.0, 0, 0.4, 5.1, 7.0, 5.0, 0.3, 0.0, 55.0];
tL = 1e-3 * [4.8, 4.8, 5.7, 9.2, 8.7, 9.4, 6.4, 5.1, 4.8];
HPly = 1e-4 * [0.205, 0.255, 5.550, 5.550, 0.255, 0.205];
TotalThickness = sum(HPly);

% Material properties
E_abs=1.7e9;
G_abs=0.6e9;
nu_abs=0.4;
h_abs=0.81;
h_left=3.89e-2;
Roh_abs=1050;
%steel
h_steel=2.03e-3;
E_steel=200e9   ;       % Young's modulus in Pa of steel
G_steel=80e6   ;       % Shear modulus in Pa
nu_steel=0.3;
Roh_steel=7900;
% CFRP
E_Grfp = 10e9  ;       % Young's modulus in Pa
E2_Grfp = 6.7e9 ;       % Young's modulus in Pa (perpendicular to fibers)
nu_Grfp = 0.26 ;      % Poisson's ratio
G_Grfp = 3.220e9;         % Shear modulus in Pa
t_Gfrp = thickness;     % Thickness of each layer in meters
h_Grfp = 8.64e-4;
Roh_Grfp=2160;
% Wood
E_wood = 14e9;            % Young's modulus in Pa for wood
G_wood = 5.5e9; 
nu_wood = 0.26;           % Poisson's ratio for wood
h_wood = thickness;      % Thickness of wood layer in meters
Roh_wood=500;

% Calculate HPlyh, HPlyc, HPlys arrays
HPlyh = [-(h_wood/3 + 2*h_Grfp/3), -(h_wood/3 + h_Grfp/3), -h_wood/3, 0, h_wood/2, (h_wood/3 + h_Grfp/3), (h_wood/3 + 2*h_Grfp/3)];
HPlyc = [-(h_wood/2 + 2*h_Grfp/3), -(h_wood/2 + h_Grfp/3), -h_wood/2, 0, h_wood/2, (h_wood/2 + h_Grfp/3), (h_wood/2 + 2*h_Grfp/3)];
HPlys = [-(h_wood/3 + 2*h_Grfp/3), -(h_wood/3 + h_Grfp/3), -h_wood/3, 0, h_wood/2, (h_wood/3 + h_Grfp/3), (h_wood/3 + 2*h_Grfp/3)];


% Create HPly matrix
HPly = zeros(length(HPlyc), length(b));
HPly(:, 1:2) = [HPlyh;HPlyh]';  % Fill first 4 columns with vector H
HPly(:, 3:8) = [HPlyc;HPlyc;HPlyc;HPlyc;HPlyc;HPlyc]';  % Fill columns 4 to 6 with vector g
HPly(:, 9:10) = [HPlyc;HPlyc]';




% Updating thickness
TotalThickness = sum(HPly);
thickness = TotalThickness;

b4 = 3e-3;
b2 = thickness - b4;
b1 = b - 2 * b4;
h1 = TotalThickness;
h2 = b4;
h4 = b4;
r1 = tL;
r4 = h4;
r2 = r4 + h2 / 2;

E_Pwood = 8.5e9;
nu_Pwood = 0.3;
G_Pwood = 0.610e9;

E_Materail = [E_Grfp, E_Grfp, E_wood, E_wood, E_Grfp, E_Grfp];
G_Materail = [G_Grfp, G_Grfp, G_wood, G_wood, G_Grfp, G_Grfp];
nu_Materail = [nu_Grfp, nu_Grfp, nu_wood, nu_wood, nu_Grfp, nu_Grfp];

leng = max(abs(YL)) * 2;
el = leng / num_elements;

angles = [0, 0, 0, 0, 0, 0];



h = 2e-3;
zjGrfp = abs((h_wood + h_Grfp / 2) - 0.00406858);
zjabs = abs(h_abs - 0.00406858);
zjwood = abs(h_wood - 0.00406858);
zjSteel = abs(h_steel - 0.00406858);

xjSteel = abs(b / 2 - b4);
xjabs = abs(b / 2 - b4);



[Mass, IPT] = mass_fun(Roh_wood, Roh_steel, Roh_abs, Roh_Grfp, b, h_steel, h_wood,h_Grfp, h_abs, el, zjGrfp, zjabs, zjwood, zjSteel, xjSteel, xjabs);

element_length = ones(num_elements, 1) * (leng / num_elements);

BendingStiff = zeros(length(b), 1);
Rcenter = zeros(length(b), 1);
j = 1;

for itr = 1:length(b)
    [A1, B1, D1, inverseM1] = calculate_ABD_matrices(E_Materail, nu_Materail, G_Materail, HPly(:, j), angles);
    [EI, RC] = calculate_bending_stiffness(inverseM1, E_steel, E_abs, b1(j), b2, b4, r1(j), r2, r4, h4, h2);
    
    BendingStiff(j) = EI(itr);
    Rcenter(j) = RC(itr);
    j = j + 1;
end

figure(2);
plot(YL, BendingStiff, 'k*', 'DisplayName', 'Present Model');
hold on;
plot(Brennanx, Brennany, 'k*', 'DisplayName', 'Brennan et Al. (Data)');
plot(BrennanEx, BrennanEy, 'k--', 'DisplayName', 'Brennan et Al. (Model)');
legend;
xlabel('X-coordinate (m)');
ylabel('Bending Stiffness EI (Nm^2)');
xlim([-0.58, 0.58]);
ylim([-100.0, 800]);
grid on;
hold off;

disp(BendingStiff)




% Function to calculate ABD matrices
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
        end
        
        z_bot = h_mat(i) / 2 * (-1 + 2 * (i - 1) / length(h_mat));
        z_top = h_mat(i) / 2 * (-1 + 2 * i / length(h_mat));
        
        A = A + Q_mat * (z_top - z_bot);
        B = B + 0.5 * Q_mat * (z_top^2 - z_bot^2);
        D = D + (1/3) * Q_mat * (z_top^3 - z_bot^3);
    end
    
    inverse_Stifness_matrix = inv([A, B; B, D]);
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

function [Mass, IPT] = mass_fun(Roh_wood, Roh_steel, Roh_abs, Roh_Grfp, b, h_steel, h_wood, h_Grfp, h_abs, el, zjGrfp, zjabs, zjwood, zjSteel, xjSteel, xjabs)
    count = 1;
    Mass = zeros(length(b), 1);
    IPT = zeros(length(b), 1);
    
    for i = 1:length(b)
        Mass(count) = (2 * (h_steel * b(count) * Roh_steel) + h_wood * b(count) * Roh_wood + h_Grfp * b(count) * Roh_Grfp + 2 * (h_abs * b(count) * Roh_abs)) * el;  % Mass
        
        TPhorGrfp = Roh_Grfp * ((1/12) * b(count)^3 * h_Grfp + (1/3) * b(count) * (3 * zjGrfp^2 * h_Grfp + 1/4 * h_Grfp^3));  % Polar mass moment of inertia Grfp
        TPhorwood = Roh_wood * ((1/12) * b(count)^3 * h_wood + (1/3) * b(count) * (3 * zjwood^2 * h_wood + 1/4 * h_wood^3));  % Polar mass moment of inertia wood
        IPhor = TPhorwood + TPhorGrfp;
        
        TPhorSteel = Roh_steel * (xjSteel(count)^2 * b(count) * h_steel + (zjSteel^2 * h_steel * b(count)) + (1/12 * h_steel^3 * b(count)) + 1/12 * h_steel * b(count)^3);  % Polar mass moment of inertia steel
        TPhorabs = Roh_abs * (xjabs(count)^2 * b(count) * h_abs + (zjabs^2 * h_abs * b(count)) + (1/12 * h_abs^3 * b(count)) + 1/12 * h_abs * b(count)^3);  % Polar mass moment of inertia abs
        IPver = 2 * TPhorSteel + 2 * TPhorabs;
        IPT(count) = IPhor + IPver;
        
        count = count + 1;
    end
end