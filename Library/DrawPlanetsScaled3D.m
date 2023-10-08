function DrawPlanetsScaled3D(r1S, r1, p1, r2S, r2, p2, sunflag)
% Description: this function creates a 3D Plot of the Trajectory from Earth
% to Mars.
% r1S = distance between Planet 1 and Sun
% r1 = radius of Planet 1
% p1 = name of Planet 1
% r2S = distance between Planet 2 and Sun
% r2 = radius of Planet 2
% p2 = name of Planet 2

switch p1
    case 'earth'
        im1 = 'earth.jpg';
    case 'mars'
        im1 = 'mars.jpg';
    case 'venus'
        im1 = 'venus.jpg';
    otherwise
        error('Invalid Planet 1')
end

switch p2
    case 'earth'
        im2 = 'earth.jpg';
    case 'mars'
        im2 = 'mars.jpg';
    case 'venus'
        im2 = 'venus.jpg';
    otherwise
        error('Invalid Planet 1')
end

rS = 696340;    % km
scaleP = 1e3;    % factor by which Planet radii are multiplied
scaleS = 1e1;    % factor by which Sun radius is multiplied

r1 = r1 * scaleP;
r2 = r2 * scaleP;
rS = rS * scaleS;

[x1, y1, z1] = sphere;
[x2, y2, z2] = sphere;
[xS, yS, zS] = sphere;

% Add the Position Vectors
x1 = r1*x1 + r1S(1);
y1 = r1*y1 + r1S(2);
z1 = r1*z1 + r1S(3);

x2 = r2*x2 + r2S(1);
y2 = r2*y2 + r2S(2);
z2 = r2*z2 + r2S(3);

nnu = 1000;
nu = linspace(0, 2*pi, nnu);

plot3(0, 0, 0, 'g*')
hold on

I1 = imread(im1);
P1 = surface(x1, y1, z1, flipud(I1), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'EdgeAlpha', 1, 'FaceAlpha', 1);

I2 = imread(im2);
P2 = surface(x2, y2, z2, flipud(I2), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'EdgeAlpha', 1, 'FaceAlpha', 1);

if sunflag
    IS = imread('sun.jpg');
    Sun = surface(rS*xS, rS*yS, rS*zS, flipud(IS), 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CDataMapping', 'direct', 'EdgeAlpha', 1, 'FaceAlpha', 1);
end

plot3(norm(r1S)*cos(nu), norm(r1S)*sin(nu), zeros(nnu, 1), 'linestyle', '-', 'color', '#c9c9c9', 'LineWidth', 0.5)
plot3(norm(r2S)*cos(nu), norm(r2S)*sin(nu), zeros(nnu, 1), 'linestyle', '-', 'color', '#c9c9c9', 'LineWidth', 0.5)

grid on
axis equal
xlabel('$x$ [km]' ,'Interpreter', 'latex', 'FontSize', 12)
ylabel('$y$ [km]' ,'Interpreter', 'latex', 'FontSize', 12)
zlabel('$z$ [km]' ,'Interpreter', 'latex', 'FontSize', 12)
view([100, 20])


end