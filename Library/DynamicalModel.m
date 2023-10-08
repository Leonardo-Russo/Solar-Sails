function dY = DynamicalModel(~, Y, C)
% Description: this function represents the dynamical model for both the
% state and the costate.
% X = [r, theta, vr, vt];
% u = alpha;
% P = [p_r p_vr p_vt];
% C = [Dsol rE rV c muS sigma We RE RV Beta];

N = length(Y);

% Retrieve Data from Input
r = Y(1);
theta = Y(2);
vr = Y(3);
vt = Y(4);
p_r = Y(5);
p_vr = Y(6);
p_vt = Y(7);

Dsol = C(1);
rES = C(2);
rVS = C(3);
c = C(4);
muS = C(5);
sigma = C(6);
We = C(7);
Beta = C(8);
savechoice = C(9);

% Initialize Derivatives Vector
dY = zeros(N, 1);

% Compute Local Variables
alpha = atan((-3*p_vr - sqrt(9*p_vr^2+8*p_vt^2))/(4*p_vt));
% alpha = atan2(-3*p_vr - sqrt(9*p_vr^2+8*p_vt^2), 4*p_vt);

% Assign Derivative Values
dY(1) = vr;
dY(2) = vt/r;
dY(3) = vt^2/r - muS/r^2 + (Beta/r^2)*(cos(alpha))^3;
dY(4) = -(vr*vt)/r + (Beta/r^2)*sin(alpha)*(cos(alpha))^2;

dY(5) = p_vr*(vt^2/r^2-2*muS/r^3+(2*Beta/r^3)*(cos(alpha))^3) + p_vt*(-vr*vt/r^2+(2*Beta/r^3)*(cos(alpha)^2)*sin(alpha));
dY(6) = -p_r + p_vt*vt/r;
dY(7) = -2*p_vr*vt/r + p_vt*vr/r;

end