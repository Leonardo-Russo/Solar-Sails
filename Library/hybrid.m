%% Solar Sails - Leonardo Russo 2015563

close all
clear all
clc

addpath('Library/')

%% Options

savechoice = 0;

% Settings for fmincon()
OptionsFMIN = optimoptions('fmincon');
OptionsFMIN.Display = 'iter-detailed';
OptionsFMIN.StepTolerance = 1e-12;
OptionsFMIN.ConstraintTolerance = 1e-18;
OptionsFMIN.Algorithm = 'sqp';
OptionsFMIN.Diagnostics = 'on';

% Settings for particleswarm()
dim = 4;    % this will be the dimension plotted during the PSO
OptionsPSO = optimoptions('particleswarm', ...
    'PlotFcn', {@pswplotbestf, @(optimValues, state) dimPlot(optimValues, state, dim)}, ...
    'HybridFcn', {@fmincon, OptionsFMIN}, ...   % hybrid function setting
    'Display', 'iter', ...
    'InertiaRange', [0.1 2.1]);

% Settings for ode113()
Tol0 = 1e-11;
Tol1 = 1e-13;
OptionsODE = odeset('RelTol', Tol0, 'AbsTol', Tol1);

%% Initialization

% Define Constant Quantities
Dsol = 86400;           % s - Solar Day
rES = 1.496e8;          % km - Earth-Sun Distance
rVS= 1.08209e8;         % km - Venus-Sun Distance
c = 299792.458;         % km/s - Light speed
muS = 132712440018;     % km^3/s^2 - Sun Gravitational Parameter

sigma = 20;             % g/m^2 - Sailing Load
We = 1361;              % W/m^2 - Solar constant
Beta = 2e-3*We*rES^2/(c*sigma);

% Define Constants Vector
C = [Dsol rES rVS c muS sigma We Beta savechoice]';

% Define Initial Conditions
r0 = rES;
theta0 = 0;
vr0 = 0;
vt0 = sqrt(muS/r0);
X0 = [r0 theta0 vr0 vt0]';

% Define Target Conditions
r_tgt = rVS;
vr_tgt = 0;
vt_tgt = sqrt(muS/r_tgt);
X_tgt = [r_tgt vr_tgt vt_tgt]';

%% PSO Optimization

mindays = 200;
% ptol = 1;

lb = [-Inf, -Inf, -Inf, mindays*Dsol];
ub = [+Inf, +Inf, +Inf, +Inf];

W0 = particleswarm(@(W0) PenalizedObjectiveFunction(W0, X0, X_tgt, C), 4, lb, ub, OptionsPSO);

P0 = W0(1:3)';
tf = W0(end);

save('Optimized Workspace.mat')


%% Results Visualization

[tspan, X] = ode113(@(t, X0) DynamicalModel(t, X0, C), [0 tf], [X0; P0], OptionsODE);

M = length(tspan);
r = X(:, 1);
theta = X(:, 2);
vr = X(:, 3);
vt = X(:, 4);
p_r = X(:, 5);
p_vr = X(:, 6);
p_vt = X(:, 7);

delta_r = r(end) - r_tgt;
delta_vr = vr(end) - vr_tgt;
delta_vt = vt(end) - vt_tgt;

alpha = zeros(M, 1);

for j = 1 : M
    alpha(j) = atan2(-3*p_vr(j) - sqrt(9*p_vr(j)^2+8*p_vt(j)^2), 4*p_vt(j));
end


var_names = ["tf", "p_r0", "p_vr0", "p_vt0", "p_rf", "p_vrf", "p_vtf"]';
units = ["days", "", "", "", "", "", ""]';
summary = [tf/Dsol, p_r(1), p_vr(1), p_vt(1), p_r(end), p_vr(end), p_vt(end)]';

fprintf('\n\n\t\t<strong>Summary of Optimization</strong>\n\n')
disp(table(var_names, units, summary, 'VariableNames',["Quantity", "Units", "Value"]))


var_names = ["delta_r", "delta_vr", "delta_vt"]';
units = ["km", "km/s", "km/s"]';
summary = [delta_r, delta_vr, delta_vt]';

fprintf('\n\n\t\t<strong>Summary of Objective Function</strong>\n\n')
disp(table(var_names, units, summary, 'VariableNames',["Quantity", "Units", "Value"]))


figure('name', 'Plot of the 2D Optimized Trajectory')
DrawTraj2D(X, [savechoice, rES, rVS])

figure('name', 'Plot of the Control Angle Alpha')
plot(tspan, rad2deg(alpha), 'LineStyle','-', 'LineWidth', 1.4)
xlabel('$t \ [s]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\alpha \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == 1
    saveas(gcf, strcat('Output\alpha.jpg'))
end

