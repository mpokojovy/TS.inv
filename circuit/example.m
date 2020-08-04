%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      (C) Michael Pokojovy (2020)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function example
    % Set random seed
    rng(1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          Constants and parameters                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Time horizon
    T = 5;
    
    % grid size
    n = 1000;

    h = T/(n + 1);
    
    t_grid = (0:n)*h;
    t_grid = t_grid';
    
    % noise intensity
    sigma = 0.1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       Observation/measurement                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % system matrices
    A = [0 -1;
         1 -2];     
    C = [0 1];
    
    % true (unknown) initial condition
    x0 = [1; 3];
    
    % true (unknown) state
    x = zeros(2, n+1);
    x = par_to_state(x0);
    
    % unnoisy observation
    y = C*x(:, 2:end);
    
    % noisy observation
    y_noisy = y + sigma*randn(size(y));
    
    sigma_hat = 1.48/sqrt(2)*median(abs(y_noisy(2:end) - y_noisy(1:(end-1))));
    
    % taut-string denoising
    y_hat_TS = TSregression(y_noisy', sigma_hat*1.149/sqrt(n))';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        Reconstructing the state                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [x0_hat_TS, resnorm, residual, exitflag] = lsqnonlin(@(x0) fun(x0, y_hat_TS), [0; 0]); 
    x_hat_TS = par_to_state(x0_hat_TS);
    
    display(['x_0 = (', num2str(x0_hat_TS(1)), ', ', num2str(x0_hat_TS(2)), ')']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          Plotting the figures                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(1);
    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 28; ySize = 14;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);
    
    % Plotting the observation
    subplot_tight(1, 2, 1, [0.1 0.075]);
    hold on;
    
    p1 = plot(t_grid(2:end), y_noisy, '--', 'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
    plot(t_grid(2:end), y_noisy, 'o', 'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
    p2 = plot(t_grid(2:end), y, 'k-', 'LineWidth', 2);
    p3 = plot(t_grid(2:end), y_hat_TS, '-', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    
    xlabel('Time $t$', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('Observation/measurement $y(t)$', 'interpreter', 'latex', 'FontSize', 18);
    
    legend([p2 p1 p3], {'Actual (unnoisy) observations', 'Noisy observations', 'Taut-string smoothed noisy observations'}, ...
                       'Location', 'Best', 'interpreter', 'latex', 'FontSize', 18);
    
    % Plotting the reconstructed state
    subplot_tight(1, 2, 2, [0.1 0.075]);
    hold on;
    
    p1 = plot(x(1, 1), x(2, 1), 'k*', 'LineWidth', 2);
    plot(x(1, :), x(2, :), 'k', 'LineWidth', 2);
    p2 = plot(x_hat_TS(1, 1), x_hat_TS(2, 1), '*', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    plot(x_hat_TS(1, :), x_hat_TS(2, :), '-', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    
    xlabel('State vector component $x_{1}(t)$', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('State vector component $x_{2}(t)$', 'interpreter', 'latex', 'FontSize', 18);
    
    legend([p1 p2], {'Actual initial state', 'Recovered initial state'}, ...
                     'Location', 'Best', 'interpreter', 'latex', 'FontSize', 18);
    
    cx = 0.50*(min([x(1, :) x_hat_TS(1, :)]) + max([x(1, :) x_hat_TS(1, :)]));
    cy = 0.50*(min([x(2, :) x_hat_TS(2, :)]) + max([x(2, :) x_hat_TS(2, :)]));
    rx = 0.60*(max([x(1, :) x_hat_TS(1, :)]) - min([x(1, :) x_hat_TS(1, :)]));
    ry = 0.60*(max([x(2, :) x_hat_TS(2, :)]) - min([x(2, :) x_hat_TS(2, :)]));
        
    axis([cx - rx cx + rx cy - ry cy + ry]);
                 
    function x_hat = par_to_state(x0)
        x_hat = zeros(size(x));
        
        for i = 1:(length(t_grid))
            t = t_grid(i);
            x_hat(:, i) = expm(A*t)*x0;
        end
    end

    function y_hat = par_to_obs(x0)
        y_hat = zeros(size(y));
        
        for i = 2:(length(t_grid))
            t = t_grid(i);
            y_hat(:, i - 1) = C*expm(A*t)*x0;
        end
    end

    function F = fun(x0, y_obs)
        y_hat = par_to_obs(x0)';
        
        F = [sigma_hat/sqrt(n)*x0; (h*cumsum(y_hat - y_obs'))*sqrt(h)];
    end
end