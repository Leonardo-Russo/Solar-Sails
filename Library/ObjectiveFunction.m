function J = ObjectiveFunction(W0, X0, X_tgt, C)
% Description: this function defines the cost function.

% Retrieve Data from Input
P0 = W0(1:3)';
tf = W0(4);

r_tgt = X_tgt(1);
vr_tgt = X_tgt(2);
vt_tgt = X_tgt(3);

Dsol = C(1);

% Propagation to Final State
Tol0 = 1e-11;
Tol1 = 1e-13;
optionsODE = odeset('RelTol',Tol0,'AbsTol',Tol1);

Y0 = [X0; P0];
[~, Y] = ode113(@(t, Y0) DynamicalModel(t, Y0, C), [0 tf], Y0, optionsODE);

r = Y(:, 1);
theta = Y(:, 2);
vr = Y(:, 3);
vt = Y(:, 4);

% Compute the Cost Function
c1 = 1;
c2 = 1e6;
c3 = 1e6;

delta_r = r(end) - r_tgt;
delta_vr = vr(end) - vr_tgt;
delta_vt = vt(end) - vt_tgt;

J = tf/Dsol + c1*abs(delta_r) + c2*abs(delta_vr) + c3*abs(delta_vt);


end

