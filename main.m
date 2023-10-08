%% Solar Sails - Leonardo Russo 2015563

close all
clear all
clc

addpath('Library/')
addpath('Planets/')

%% Options

savechoice = 0;

% Settings for particleswarm()
OptionsPSO = optimoptions('particleswarm', ...
    'PlotFcn', {@pswplotbestf, @(optimValues, state) dimPlot(optimValues, state, 4)}, ...
    'SwarmSize', 1000, ...
    'Display', 'iter', ...
    'InertiaRange', [0.1 1.1], ...
    'UseParallel', true, ...
    'MaxIterations', 1000);

% Settings for fmincon()
OptionsFMIN = optimoptions('fmincon', ...
    'PlotFcn', 'optimplotfval', ...
    'Display', 'iter-detailed', ...
    'StepTolerance', 1e-12, ...
    'ConstraintTolerance', 1e-6, ...
    'Algorithm', 'sqp', ...
    'Diagnostics', 'off', ...
    'MaxFunctionEvaluations', 10000);

% Settings for ode113()
Tol0 = 1e-11;
Tol1 = 1e-13;
OptionsODE = odeset('RelTol', Tol0, 'AbsTol', Tol1);


%% Single Arc Model

% Define Constant Quantities
Dsol = 86400;           % s - Solar Day
rES = 1.496e8;          % km - Earth-Sun Distance
rVS= 1.08209e8;         % km - Venus-Sun Distance
c = 299792.458;         % km/s - Light speed
muS = 132712440018;     % km^3/s^2 - Sun Gravitational Parameter

rE = 6378.136;          % km - Earth Radius
rV = 6051.8;            % km - Venus Radius

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
maxdays = 300;

lb = [100, 100, 1e5, mindays*Dsol];
ub = [400, 400, 1e10, maxdays*Dsol];

W0 = particleswarm(@(W0) ObjectiveFunction(W0, X0, X_tgt, C), 4, lb, ub, OptionsPSO);

P0 = W0(1:3)';
tf = W0(end);

save('WorkspacePSO.mat')


%% FMINCON Optimization

A = [];
B = [];
Aeq = [];
Beq = [];

mindays = 200;
maxdays = 300;

lb = [-Inf, -Inf, -Inf, mindays*Dsol];
ub = [Inf, Inf, Inf, maxdays*Dsol];

max_iterations = 5;
obj_tgt = 1;

for i = 1 : max_iterations

    [W0, Jf, exitflag] = fmincon(@(W0) ObjectiveFunction(W0, X0, X_tgt, C), W0, A, B, Aeq, Beq, lb, ub, @(W0) nonlinconstr(W0, X0, X_tgt, C), OptionsFMIN);

    if Jf < obj_tgt
        break
    end

end

P0 = W0(1:3)';
tf = W0(end);

save('WorkspaceFMIN.mat')


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
    alpha(j) = atan((-3*p_vr(j) - sqrt(9*p_vr(j)^2+8*p_vt(j)^2))/(4*p_vt(j)));
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


r1S = [X(1, 1)*cos(X(1, 2)), X(1, 1)*sin(X(1, 2)), 0]';
r2S = [X(end, 1)*cos(X(end, 2)), X(end, 1)*sin(X(end, 2)), 0]';
figure('name', 'Plot of the 3D Optimized Trajectory')
DrawPlanetsScaled3D(r1S, rE, 'earth', r2S, rV, 'venus', 1)
DrawTraj2D(X, [savechoice, rES, rVS])
if savechoice == 1
    saveas(gcf, strcat('Output\traj.jpg'))
end

figure('name', 'Plot of the Control Angle Alpha')
plot(tspan/Dsol, rad2deg(alpha), 'LineStyle','-', 'LineWidth', 1.4)
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\alpha \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
if savechoice == 1
    saveas(gcf, strcat('Output\alpha.jpg'))
end


% input('Press Enter to continue...\n');


%% Three Arcs Model

% clc
% close all
% clear all
% load('Archive/WorkspaceFMIN.mat')

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

% Settings for fmincon()
OptionsFMIN = optimoptions('fmincon', ...
    'PlotFcn', {@optimplotfval}, ...
    'Display', 'iter-detailed', ...
    'StepTolerance', 1e-24, ...
    'ConstraintTolerance', 1e-12, ...
    'Algorithm', 'sqp', ...
    'Diagnostics', 'off', ...
    'MaxFunctionEvaluations', 10000);


% Settings for particleswarm()
OptionsPSO = optimoptions('particleswarm', ...
    'PlotFcn', {@pswplotbestf, @alphaTimesPlot}, ...
    'HybridFcn', {@fmincon, OptionsFMIN}, ...   % hybrid function setting
    'SwarmSize', 1000, ...
    'Display', 'iter', ...
    'InertiaRange', [0.1 1.6], ...
    'UseParallel', true, ...
    'MaxIterations', 1000);


mindays = 60;
maxdays = 140;

% lb = [-2*pi, -2*pi, -2*pi, mindays*Dsol, mindays*Dsol, mindays*Dsol];
% ub = [2*pi, 2*pi, 2*pi, maxdays*Dsol, maxdays*Dsol, maxdays*Dsol];

lb = [-deg2rad(60), -deg2rad(60), -deg2rad(60), mindays*Dsol, mindays*Dsol, mindays*Dsol];
ub = [-deg2rad(20), -deg2rad(20), -deg2rad(20), maxdays*Dsol, maxdays*Dsol, maxdays*Dsol];

W0 = particleswarm(@(W0) ObjectiveFunction3Arcs(W0, X0, X_tgt, C), 6, lb, ub, OptionsPSO);

save('WorkspaceHybrid.mat')


%% Results Visualization

alpha1 = W0(1);
alpha2 = W0(2);
alpha3 = W0(3);
t1 = W0(4);
t2 = W0(5);
t3 = W0(6);

% First Arc
C1 = [C; alpha1];
[tspan1, X1] = ode113(@(t, X0) DynamicalModel3Arcs(t, X0, C1), [0 t1], X0, OptionsODE);
X1f = X1(end, :);

% Second Arc
C2 = [C; alpha2];
[tspan2, X2] = ode113(@(t, X0) DynamicalModel3Arcs(t, X0, C2), [t1 t1+t2], X1f, OptionsODE);
X2f = X2(end, :);

% Second Arc
C3 = [C; alpha3];
[tspan3, X3] = ode113(@(t, X0) DynamicalModel3Arcs(t, X0, C3), [t1+t2 t1+t2+t3], X2f, OptionsODE);
X3f = X3(end, :);

tspan = [tspan1; tspan2; tspan3];
X = [X1; X2; X3];

M = length(tspan);
r = X(:, 1);
theta = X(:, 2);
vr = X(:, 3);
vt = X(:, 4);

delta_r = r(end) - r_tgt;
delta_vr = vr(end) - vr_tgt;
delta_vt = vt(end) - vt_tgt;


var_names = ["alpha1", "alpha2", "alpha3", "t1", "t2", "t3", "tf"]';
units = ["deg", "deg", "deg", "days", "days", "days", "days"]';
summary = [rad2deg(alpha1), rad2deg(alpha2), rad2deg(alpha3), t1/Dsol, t2/Dsol, t3/Dsol, (t1+t2+t3)/Dsol]';

fprintf('\n\n\t\t<strong>Summary of Optimization</strong>\n\n')
disp(table(var_names, units, summary, 'VariableNames',["Quantity", "Units", "Value"]))


var_names = ["delta_r", "delta_vr", "delta_vt"]';
units = ["km", "km/s", "km/s"]';
summary = [delta_r, delta_vr, delta_vt]';

fprintf('\n\n\t\t<strong>Summary of Objective Function</strong>\n\n')
disp(table(var_names, units, summary, 'VariableNames',["Quantity", "Units", "Value"]))


r1S = [X(1, 1)*cos(X(1, 2)), X(1, 1)*sin(X(1, 2)), 0]';
r2S = [X(end, 1)*cos(X(end, 2)), X(end, 1)*sin(X(end, 2)), 0]';
figure('name', 'Plot of the 3D Optimized Trajectory')
DrawPlanetsScaled3D(r1S, rE, 'earth', r2S, rV, 'venus', 1)
DrawTraj2D(X, [savechoice, rES, rVS])
if savechoice == 1
    saveas(gcf, strcat('Output\traj3Arcs.jpg'))
end


figure('name', 'Plot of the Control Angle Alpha')
hold on
grid on
plot(tspan1/Dsol, rad2deg(alpha1)*ones(length(tspan1), 1), 'LineStyle','-', 'LineWidth', 1.4, 'Color', '#22bf83')
plot(tspan2/Dsol, rad2deg(alpha2)*ones(length(tspan2), 1), 'LineStyle','-', 'LineWidth', 1.4, 'Color', '#227bbf')
plot(tspan3/Dsol, rad2deg(alpha3)*ones(length(tspan3), 1), 'LineStyle','-', 'LineWidth', 1.4, 'Color', '#ff7403')
plot([tspan1(1), tspan1(end)]./Dsol, rad2deg(alpha1*ones(2, 1)), 'LineStyle','none', 'Marker','.', 'MarkerSize',10, 'Color', '#22bf83')
plot([tspan2(1), tspan2(end)]./Dsol, rad2deg(alpha2*ones(2, 1)), 'LineStyle','none', 'Marker','.', 'MarkerSize',10, 'Color', '#227bbf')
plot([tspan3(1), tspan3(end)]./Dsol, rad2deg(alpha3*ones(2, 1)), 'LineStyle','none', 'Marker','.', 'MarkerSize',10, 'Color', '#ff7403')
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\alpha_i \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\alpha_1$', '$\alpha_2$', '$\alpha_3$', 'Interpreter','latex', 'FontSize', 12, 'location', 'best')
if savechoice == 1
    saveas(gcf, strcat('Output\alpha3arcs.jpg'))
end

