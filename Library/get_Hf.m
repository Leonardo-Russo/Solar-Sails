function Hf = get_Hf(Pf, Xf, C)
% Description: this functions takes as input the raw output of the
% integration and evaluates Hf.

% Import Data from Input
p_rf = Pf(1);
p_tf = 0;
p_vrf = Pf(2);
p_vtf = Pf(3);
p_kf = 1;

rf = Xf(1);
thetaf = Xf(2);
vrf = Xf(3);
vtf = Xf(4);

muS = C(5);
Beta = C(8);

% Compute Local Variables
alpha = atan((-3*p_vrf - sqrt(9*p_vrf^2+8*p_vtf^2))/(4*p_vtf));


% Evaluate Hf
Hf = p_rf*vrf + p_tf*vtf/rf + p_vrf*(vtf^2/rf - muS/rf^2 + Beta/rf^2*cos(alpha)^3) + p_vtf*(-vrf*vtf/rf + Beta/rf^2*cos(alpha)^2*sin(alpha)) - p_kf;


end