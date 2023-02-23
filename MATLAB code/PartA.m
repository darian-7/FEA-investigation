clear
clc

% Darian Irani
% 26 Nov 2022
% Part A
% Convergence plot of No. of elements vs Von Mises Stress
% Quadratic hexahedral mesh

%% No. of elements vs Von Mises Stress (MPa)
Elements_Stress = readtable('elements_stress','Sheet','Sheet1');

x = Elements_Stress.Elements;
y = Elements_Stress.Stress;

scatter(x,y,'Filled');
yline(214.1)    % Analytical sol
axis ([0 35000 100 230])
xlabel('Number of elements');
ylabel('Von Mises stress (MPa)');
xticks([0 5000 10000 15000 20000 25000 32000]);
xticklabels({'0', '5000', '10000', '15000', '20000', '25000', '32000'})
grid on
figure

%% No. of elements vs Tip deflection (mm)

Elements_Deflection = readtable('elements_def','Sheet','Sheet1');

x = Elements_Deflection.Elements;
y = Elements_Deflection.Deflection;

scatter(x,y,"Filled");
axis ([0 35000 0.1478 0.1498])
xlabel('Number of elements');
ylabel('Tip deflection (mm)');
xticks([0 5000 10000 15000 20000 25000 32000]);
xticklabels({'0', '5000', '10000', '15000', '20000', '25000', '32000'})
grid on

