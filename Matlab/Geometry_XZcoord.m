function [PtXcoo,PtZcoo,PtZcoou,Xf_filtered,Lc]=Geometry_XZcoord(Lc,Ls,Lh,Hs,Hf,Hh)
%% Function used to generat the x and z coordinates


% Heel Curve
[Xh, Zh] = heel_curve(Lh, Hh, Lc);

% Shovel Curve
[Xs, Zs] = shovel_curve(Ls, Hs, Lc);

% Camber Curve
rh = (Lh^2 + Hh^2) / (2 * Hh);
rs = (Ls^2 + Hs^2) / (2 * Hs);
[Xf, Zf, Xf0, Zf0] = camber_curve(Lc, Hf, rh, rs, Xh, Xs);

% Filter camber curve
Xf_filtered = Xf(Zf >= 0);
Zf_filtered = Zf(Zf >= 0);
Xf_filtered = flip(Xf_filtered);
Zf_filtered = flip(Zf_filtered);

% Plot results
figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
plot(Xh, Zh, '*', 'DisplayName', 'Heel Curve');
hold on;
plot(Xs, Zs, '*', 'DisplayName', 'Shovel Curve');
plot(Xf_filtered, Zf_filtered, 'DisplayName', 'Camber Curve');
xlabel('X Coordinate');
ylabel('Z Coordinate');
title('Shovel Design');
xlim([-1500e-3, 1500e-3]);
ylim([-300e-3, 300e-3]);
grid on;
legend;

% Concatenate results
PtXcoo = [Xh, Xf_filtered, Xs];
PtZcoo = [Zh, Zf_filtered, Zs];
th = 4e-3;
ts = 4e-3;
tbody = 8e-3;
PtXcoou = [Xh, Xf_filtered, Xs];
PtZcoou = [Zh + th, Zf_filtered + tbody, Zs + ts];

% Plot results
figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
plot(PtXcoo, PtZcoo, 'k*');
hold on;
plot(PtXcoou, PtZcoou, 'k');
xlabel('X Coordinate');
ylabel('Z Coordinate');
ylim([-50e-3, 50e-3]);
xlim([-1.50, 1.50]);
grid on;

% Save results
writematrix(PtXcoo', 'testX.txt', 'Delimiter', ',');
writematrix(PtZcoo', 'testZ.txt', 'Delimiter', ',');
writematrix(PtXcoou', 'testXu.txt', 'Delimiter', ',');
writematrix(PtZcoou', 'testZu.txt', 'Delimiter', ',');

% Functions
function [Xh, Zh] = heel_curve(Lh, Hh, Lc)
    rh = (Lh^2 + Hh^2) / (2 * Hh);
    alfah = linspace(pi + acos(Lh / rh), 3 * pi / 2, 10);
    Xh = -Lc / 2 + rh * cos(alfah);
    Zh = rh + rh * sin(alfah);
end

function [Xs, Zs] = shovel_curve(Ls, Hs, Lc)
    rs = (Ls^2 + Hs^2) / (2 * Hs);
    alfas = linspace(3 * pi / 2, 2 * pi - acos(Ls / rs), 10);
    Xs = Lc / 2 + rs * cos(alfas);
    Zs = rs + rs * sin(alfas);
end

function [Xf, Zf, Xf0, Zf0] = camber_curve(Lc, Hf, rh, rs, Xh, Xs)
    Xf0 = ((Lc^2 - 2 * rs * Hf + 2 * rh * Hf) / (2 * Lc)) - Lc / 2;
    rf = ((Xf0 + Lc / 2)^2 - 2 * rh * Hf + Hf^2) / (2 * Hf);
    Zf0 = Hf - rf;
    Xhc = Xh(end);
    Xsc = Xs(1);
    gammah = asin(1 - Hf / rf);
    gammas = asin(1 - Hf / rf);
    gamma = linspace(gammas, -gammah + pi, 20);
    Xf = Xf0 + rf * cos(gamma);
    Zf = Zf0 + rf * sin(gamma);
end

end