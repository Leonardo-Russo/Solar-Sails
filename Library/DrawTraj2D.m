function DrawTraj2D(X, P)
% Description: this function plots the 2D orbit of the S/C during the
% manoeuvre.

gif = 0;

r = X(:, 1);
theta = X(:, 2);
savechoice = P(1);
rES = P(2);
rVS = P(3);
M = size(X, 1);


% Create Visuals
m = 5000;
theta_span = linspace(0, 2*pi, m);
plot3(rES*cos(theta_span), rES*sin(theta_span), zeros(m, 1), 'Color','#3061ff')
hold on
plot3(rVS*cos(theta_span), rVS*sin(theta_span), zeros(m, 1), 'Color','#3cc0e8')
xlabel('$x \ [km]$', 'Interpreter','latex')
ylabel('$y \ [km]$', 'Interpreter','latex')
zlabel('$z \ [km]$', 'Interpreter','latex')

% Static Figure
if ~gif
    
    plot3(r.*cos(theta), r.*sin(theta), zeros(M, 1), 'Color', '#ff9421', 'LineStyle','-', 'LineWidth', 2)
    axis equal

end


% % Dynamic Figure
% if gif
% 
%     for i = 1 : M
% 
%         plot(r(i)*cos(theta(i)), r(i)*sin(theta(i)), 'Marker','.', 'MarkerSize',5, 'Color', '#ff9421')
%         xlabel('$x \ [km]$', 'Interpreter','latex')
%         ylabel('$y \ [km]$', 'Interpreter','latex')
%         if savechoice == 1
%             exportgraphics(gcf, 'Output/manoeuvre.gif', 'Append', true)
%         else
%             pause(0.01)
%         end
% 
%     end
% 
%     plot(r.*cos(theta), r.*sin(theta), 'Color', '#ff9421', 'LineStyle','-', 'LineWidth', 2)
%     legend('Terra', 'Venere')
% 
%     if savechoice == 1
%         exportgraphics(gcf, 'Output/manoeuvre.gif', 'Append', true)
%     end
% 
% end


end