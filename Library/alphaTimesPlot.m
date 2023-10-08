function stop = alphaTimesPlot(optimValues, state)
% Description: this function is a custom plot to show the evolution
%  of a given dimension during the optimization process.

stop = false;

gifFilename = 'alphatimes.gif';

persistent alf1 alf2 alf3

% alpha1 = optimValues.swarm(:, 1);
% alpha2 = optimValues.swarm(:, 2);
% alpha3 = optimValues.swarm(:, 3);
% t1 = optimValues.swarm(:, 4);
% t2 = optimValues.swarm(:, 5);
% t3 = optimValues.swarm(:, 6);

alpha1 = optimValues.bestx(1);
alpha2 = optimValues.bestx(2);
alpha3 = optimValues.bestx(3);
t1 = optimValues.bestx(4);
t2 = optimValues.bestx(5);
t3 = optimValues.bestx(6);

Dsol = 86400;
    
% Plot
if optimValues.iteration > 0
    delete(alf1)
    delete(alf2)
    delete(alf3)
end
grid on 
hold on

alf1 = plot([0, t1]./Dsol, rad2deg(alpha1*ones(2, 1)), 'LineStyle','-', 'LineWidth', 1.5, 'Marker','.', 'MarkerSize',10, 'Color', '#22bf83');
alf2 = plot([t1, t1+t2]./Dsol, rad2deg(alpha2*ones(2, 1)), 'LineStyle','-', 'LineWidth', 1.5,  'Marker','.', 'MarkerSize',10, 'Color', '#227bbf');
alf3 = plot([t1+t2, t1+t2+t3]./Dsol, rad2deg(alpha3*ones(2, 1)), 'LineStyle','-', 'LineWidth', 1.5,  'Marker','.', 'MarkerSize',10, 'Color', '#ff7403');
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 12)
ylabel('$\alpha_i \ [deg]$', 'Interpreter','latex', 'FontSize', 12)
legend('$\alpha_1$', '$\alpha_2$', '$\alpha_3$', 'Interpreter','latex', 'FontSize', 12, 'location', 'best')
axis([0, 300, -361, 361])
title('Optimization of Alphas');
drawnow;

% Capture the plot as a frame
frame = getframe(gcf);
[imind, cm] = rgb2ind(frame.cdata, 256);

% Write to the GIF file
if optimValues.iteration == 0
    imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
else
    imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
end

    
end
