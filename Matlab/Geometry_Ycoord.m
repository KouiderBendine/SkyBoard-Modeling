function [PtXcoob,PtYcoob]=Geometry_Ycoord(Xf_filtered,PtZcoo,Lc,bm,Wh,Ws)
%% Function used to generat the y  coordinate

% Calculate coordinates
[PtXcoob, PtYcoob] = calculate_coordinates(bm, Lc, Wh, Ws, Xf_filtered);

% Plot results
figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
hold on;
for i = 0
    [PtXcoob, PtYcoob] = calculate_coordinates(bm + i * 1e-3, Lc, Wh, Ws, Xf_filtered);
    plot(PtXcoob, PtYcoob, '*');
    plot(PtXcoob, -PtYcoob, '*');
end
xlabel('X Coordinate');
ylabel('Y Coordinate');
xlim([-1000e-3, 1000e-3]);
ylim([-300e-3, 300e-3]);
grid on;
hold off;

% Plot 3D curve
figure('Color', 'w')
set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times New Roman')
hold on;
plot3(PtXcoob, PtYcoob, PtZcoo, 'k*');
plot3(PtXcoob, -PtYcoob, PtZcoo, 'k*');
%plot3(PtXcoob, PtYcoob, PtZcoou, 'k*');
%plot3(PtXcoob, -PtYcoob, PtZcoou, 'k*');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');
view(50, 25);
xlim([-1000e-3, 1000e-3]);
ylim([-800e-3, 800e-3]);
zlim([-50e-3, 100e-3]);
grid on;
hold off;

% Save results
writematrix(PtXcoob', 'testXTop.txt', 'Delimiter', ',');
writematrix(PtYcoob', 'testYTop.txt', 'Delimiter', ',');

% Functions
function [PtXcoob, PtYcoob] = calculate_coordinates(bm, Lc, Wh, Ws, Xf_filtered)
    Whs = Ws;
    
    Rscu = ((Lc)^2 + bm^2 - Whs^2) / (4 * (Whs - bm));
    Xscu = 0;
    Yscu = Rscu + bm / 2;
    betahu = asin(Rscu / (Rscu + Wh / 2));
    betasu = asin(Rscu / (Rscu + Ws / 2));
    betau = linspace(pi + 10.5 * betasu / 8.5, 2 * pi - 10.5 * betasu / 8.5, length(Xf_filtered));
    Xscu = Xscu + Rscu * cos(betau);
    Yscu = Yscu + Rscu * sin(betau);

    Rsc = ((Lc)^2 + bm^2 - Whs^2) / (4 * (Whs - bm));
    Xsc = 0;
    Ysc = Rsc + bm / 2;
    betah = asin(Rsc / (Rsc + Wh / 2));
    betas = asin(Rsc / (Rsc + Ws / 2));
    beta = linspace(2 * pi - betas, pi + betah, length(Xf_filtered));
    Xsc = -Xsc - Rsc * cos(beta);
    Ysc = -Ysc - Rsc * sin(beta);

    betaheel = linspace(betah, pi, 11);
    Xht = -Lc / 2 + Wh / 2 * cos(betaheel);
    Yht = Wh / 2 * sin(betaheel);
    Xht(1) = [];
    Yht(1) = [];
    betash = linspace(betas, pi, 11);
    Xst = Lc / 2 - Ws / 2 * cos(betash);
    Yst = Ws / 2 * sin(betash);
    Xst(1) = [];
    Yst(1) = [];
    Xht = flip(Xht);
    Yht = flip(Yht);

    PtXcoob = [Xht, Xscu, Xst];
    PtYcoob = [Yht, Yscu, Yst];
end
end