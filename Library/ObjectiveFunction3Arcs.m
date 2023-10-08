function J = ObjectiveFunction3Arcs(W0, X0, X_tgt, C)
% Description: this function defines the cost function in the case of three
% arcs being model being used.

% Retrieve Data from Input
alpha1 = W0(1);
alpha2 = W0(2);
alpha3 = W0(3);
t1 = W0(4);
t2 = W0(5);
t3 = W0(6);

r_tgt = X_tgt(1);
vr_tgt = X_tgt(2);
vt_tgt = X_tgt(3);

Dsol = C(1);

tf_max = 300;
tf_min = 270;

% Propagation to Final State
Tol0 = 1e-11;
Tol1 = 1e-13;
OptionsODE = odeset('RelTol',Tol0,'AbsTol',Tol1);


% First Arc
C1 = [C; alpha1];
[~, X1] = ode113(@(t, X0) DynamicalModel3Arcs(t, X0, C1), [0 t1], X0, OptionsODE);
X1f = X1(end, :);

% Second Arc
C2 = [C; alpha2];
[~, X2] = ode113(@(t, X0) DynamicalModel3Arcs(t, X0, C2), [t1 t1+t2], X1f, OptionsODE);
X2f = X2(end, :);

% Second Arc
C3 = [C; alpha3];
[~, X3] = ode113(@(t, X0) DynamicalModel3Arcs(t, X0, C3), [t1+t2 t1+t2+t3], X2f, OptionsODE);
X3f = X3(end, :);

r3f = X3f(1);
vr3f = X3f(3);
vt3f = X3f(4);

% Compute the Cost Function
c1 = 1;
c2 = 1e6;
c3 = 1e6;

delta_r = r3f - r_tgt;
delta_vr = vr3f - vr_tgt;
delta_vt = vt3f - vt_tgt;


J = (t1 + t2 + t3)/Dsol + c1*abs(delta_r)+c2*abs(delta_vr)+c3*abs(delta_vt);

if (t1+t2+t3)/Dsol > tf_max || (t1+t2+t3)/Dsol < tf_min
    J = J + 1e20;
end

end