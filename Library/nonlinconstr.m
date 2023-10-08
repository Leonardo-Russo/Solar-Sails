function [ineq, equal] = nonlinconstr(W0, X0, X_tgt, C)
% Description: this function takes as input the unknown quantities inside X
% and returns the vector of the NON-LINEAR inequality and equality
% constraints.

% Matlab Help: accepts X and returns the vectors C and Ceq, representing the nonlinear 
% inequalities and equalities respectively. fmincon minimizes FUN such 
% that C(X) <= 0 and Ceq(X) = 0.

tol_ineq = 1e-3;

% Import Data from Input
P0 = W0(1:3)';
tf = W0(4);

Y0 = [X0; P0];

% Perform the Integration to Obtain Final Values
Tol0 = 1e-9;
Tol1 = 1e-11;
OptionsODE = odeset('RelTol', Tol0, 'AbsTol',Tol1);

[~, Y] = ode113(@(t, Y0) DynamicalModel(t, Y0, C), [0, tf], Y0, OptionsODE);

% Store Final Results
Yf = Y(end, :);
Xf = Yf(1:4)';
Pf = Yf(5:7)';

% Compute Hf
Hf = get_Hf(Pf, Xf, C);


% Assign Output Values
ineq = -(Hf + tol_ineq);
equal = [Xf(1); Xf(3:4)] - X_tgt;

% equal = [equal; Hf - 1];
% ineq = [];


end