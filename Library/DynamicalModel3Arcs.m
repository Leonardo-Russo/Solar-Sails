function dY = DynamicalModel3Arcs(~, Y, C)
% Description: this function represents the dynamical model for both the
% state and the costate.
% X = [r, theta, vr, vt];
% C = [Dsol rE rV c muS sigma We RE RV Beta ... alpha];

N = length(Y);

% Retrieve Data from Input
r = Y(1);
theta = Y(2);
vr = Y(3);
vt = Y(4);

Dsol = C(1);
rES = C(2);
rVS = C(3);
c = C(4);
muS = C(5);
sigma = C(6);
We = C(7);
Beta = C(8);
savechoice = C(9);
alpha = C(end);

% Initialize Derivatives Vector
dY = zeros(N, 1);

% Assign Derivative Values
dY(1) = vr;
dY(2) = vt/r;
dY(3) = vt^2/r - muS/r^2 + (Beta/r^2)*(cos(alpha))^3;
dY(4) = -(vr*vt)/r + (Beta/r^2)*sin(alpha)*(cos(alpha))^2;


end