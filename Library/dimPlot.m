function stop = dimPlot(optimValues, state, dim)
% Description: this function is a custom plot to show the evolution
%  of a given dimension during the optimization process.

    stop = false;
    
    persistent iterationHistory; % To keep track of values across iterations
    gifFilename = 'time_opt.gif';
    
    if isempty(iterationHistory)
        iterationHistory = [];
    end
    
    if strcmp(state,'init')
        iterationHistory = [];
    elseif strcmp(state,'iter')
        % Extract particle positions for the specified dimension
        dimensionValues = optimValues.swarm(:, dim);
        
        % Append to history
        iterationHistory = [iterationHistory; mean(dimensionValues)];

        grid on
        
        % Plot
        if dim == 4         % since 4th dimension is time in seconds
            plot(iterationHistory/86400, 'b-', 'marker', '.', 'markersize', 10);
        else
            plot(iterationHistory, 'b-', 'marker', '.', 'markersize', 10);
        end

        xlabel('Iteration');
        ylabel(['Values of Dimension ', num2str(dim)]);
        title(['Evolution of Dimension ', num2str(dim)]);
        drawnow;

        % Capture the plot as a frame
        frame = getframe(gcf);
        [imind, cm] = rgb2ind(frame.cdata, 256);

        % Write to the GIF file
        if isempty(iterationHistory) || length(iterationHistory) == 1
            imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
        else
            imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
        end
    end

end
