function DrawTraj3D(rTraj)
% Description: this function creates a 3D Plot of the Trajectory from Earth
% to Mars.

x = rTraj(:, 1);
y = rTraj(:, 2);
z = rTraj(:, 3);

plot3(x, y, z,'Color','#ff7403', 'Linestyle', '-', 'linewidth', 0.1)
hold on
plot3(x(1), y(1), z(1), 'Color', '#ff7403', 'LineStyle','none', 'marker', '.', 'markersize', 15)
plot3(x(end), y(end), z(end), 'Color', '#ff7403', 'LineStyle','none', 'marker', '.', 'markersize', 15)

end