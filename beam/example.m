%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      (C) Michael Pokojovy (2020)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function example;
    % Set random seed
    rng(123);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          Constants and parameters                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % grid size
    n = 100;
    
    % elastic modulus (in GPa), moment inertia in torsion (in m^4) and
    % density (in kg/m) of a European HE 100A I-beam
    E  = 200E9;
    I  = 5.24E-8;
    rho  = 16.7;
    
    % flexural rigidity
    FR = E*I;
    
    % beam length (in m)
    L = 3.6;

    h = L/(n + 1);
    
    x_grid = (1:n)*h;
    x_grid = x_grid';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                              Volumetic loads                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_r = 2E3*9.8/4.0581/L;
    
    % actual loads - a total load of 2 tons
    f   = (x_grid/L + exp(-(x_grid - 2*L/3).^2/(2*(L/8)^2)));
    
    % noisy loads 
    sigma = 0.2;
    eps   = randn(n, 1)*sigma;
    
    f_n = (f + eps)*f_r;
    
    % noisy corrupted loads
    alpha      = 0.05;
	ind        = randi(n, floor(n*alpha), 1);
    eps_p      = eps;    
    eps_p(ind) = eps(ind) + 1.5;
    
    f_p = (f + eps_p)*f_r;
    f   = f*f_r;
    
    sigma_n = 1.48/sqrt(2)*median(f_n(2:end) - f_n(1:(end - 1)))*10;
    sigma_p = 1.48/sqrt(2)*median(f_p(2:end) - f_p(1:(end - 1)))*10;
    
    % taut-string denoising
    f_n_TS = TSregression(f_n, sigma_n*1.149/sqrt(n));
    f_p_TS = TSregression(f_p, sigma_p*1.149/sqrt(n));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        Reconstructing the state                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D  = sparse(1:n,1:n,-2*ones(1,n),n,n);
    E  = sparse(2:n,1:n-1,ones(1,n-1),n,n);
    DL = E+D+E';
    % Dirichlet Laplacian
    DL = DL/h^2; 
    
    w      = actual_solution(f_n);
    
    w_p    = par_estimator(f_p);    
    w_n_TS = par_estimator(f_n_TS);
    w_p_TS = par_estimator(f_p_TS);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                          Plotting the figures                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%% Figure 1
    figure(1);
    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 28; ySize = 14;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);
    
    % Plotting the loads
    
    subplot_tight(1, 2, 1, [0.1 0.075]);
    hold on;
    p1 = plot(x_grid, f_n,    'o',  'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
         plot(x_grid, f_n,    '--', 'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
    p2 = plot(x_grid, f,      'k',  'LineWidth', 2);
    p3 = plot(x_grid, f_n_TS, '-',  'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    
    c = 0.5*(min(f_n) + max(f_n));
    r = 1.0*(max(f_n) - min(f_n));
        
    axis([0 L c - r c + r]);
    
    xlabel('Position $x$ on midline [m]', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('Distributed loads [N/m]', 'interpreter', 'latex', 'FontSize', 18);
    
    legend([p2 p1 p3], {'Actual loads', 'Noisy measurements', 'Taut-string smoothing'}, ...
                        'Location', 'Best', 'interpreter', 'latex', 'FontSize', 18);
    
    % Plotting the deflection 
    
    subplot_tight(1, 2, 2, [0.1 0.075]);
    hold on;

    plot(x_grid, w,      'k', 'LineWidth', 2);
    plot(x_grid, w_n_TS, '-', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);;
    
    c = 0.5*(min(w) + max(w));
    r = 1.0*(max(w) - min(w));
        
    axis([0 L c - r c + r]);
    
    xlabel('Position $x$ on midline [m]', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('Vertical displacement of midline [m]', 'interpreter', 'latex', 'FontSize', 18);
    
    legend({'Actual deflection', 'Estimated deflection'}, ...
            'Location', 'Best', 'interpreter', 'latex', 'FontSize', 18);    
    
    %%%% Figure 2
    figure(2);
    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 28; ySize = 14;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);
    
    % Plotting the loads
    
    subplot_tight(1, 2, 1, [0.1 0.075]);
    hold on;
    p1 = plot(x_grid(setdiff(1:n, ind)), f_p(setdiff(1:n, ind)), 'o', 'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
    p2 = plot(x_grid(ind),               f_p(ind),               '*', 'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
         plot(x_grid,                    f_p,                    '-', 'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
    p3 = plot(x_grid,                    f,                      'k', 'LineWidth', 2);
    p4 = plot(x_grid,                    f_p_TS,                 '-', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    
    c = 0.5*(min(f_p) + max(f_p));
    r = 1.0*(max(f_p) - min(f_p));
        
    axis([0 L c - r c + r]);
    
    xlabel('Position $x$ on midline [m]', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('Distributed loads [N/m]', 'interpreter', 'latex', 'FontSize', 18);
    
    legend([p3, p1, p2, p4], {'Actual loads', '95\% noisy measurements', '5\% corrupted measurements', 'Taut-string smoothing'}, ...
                              'Location', 'Best', 'interpreter', 'latex', 'FontSize', 18);
    
    % Plotting the deflection
    
    subplot_tight(1, 2, 2, [0.1 0.075]);
    hold on;

    plot(x_grid, w,      'k',  'LineWidth', 2);
    plot(x_grid, w_p,    '--', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    plot(x_grid, w_p_TS, '-',  'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    
    c = 0.5*(min(w) + max(w));
    r = 1.0*(max(w) - min(w));
        
    axis([0 L c - r c + r]);
    
    xlabel('Position $x$ on midline [m]', 'interpreter', 'latex', 'FontSize', 18);
    ylabel('Vertical displacement of midline [m]', 'interpreter', 'latex', 'FontSize', 18);
    
    legend({'Actual deflection', 'Estimated deflection w/o taut string', 'Estimated deflection w/ taut string'}, ...
            'Location', 'Best', 'interpreter', 'latex', 'FontSize', 18);     
        
    function par = par_estimator(f)
        w0 = ones(size(x_grid));
    
        [par, resnorm, residual, exitflag] = lsqnonlin(@energy_functional, w0); 
        
        function res = energy_functional(w)
            res = [h*DL*w*sigma*f_r/FR/sqrt(n); h*(DL*w - DL\(f/FR))];
        end
    end

    function sol = actual_solution(f)
        sol = DL\(DL\(f/FR));
    end
end