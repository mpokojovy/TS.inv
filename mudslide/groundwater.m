%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      (C) Michael Pokojovy (2020)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function groundwater
    global gwl_data rf_data time_data;
    
    N_historic = ceil(length(time_data)*0.5);
    
    gwl_historic_data  = gwl_data(1:N_historic);
    rf_historic_data   = rf_data(1:N_historic);
    time_historic_data = time_data(1:N_historic);
    
    gwl_future_data  = gwl_data((N_historic+1):end);
    rf_future_data   = rf_data((N_historic+1):end);
    time_future_data = time_data((N_historic+1):end);

    tau = 3;
    dt  = 1;
    
    n = N_historic - tau;
    
    m = length(time_data) - N_historic;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Groundwater level data                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h = sqrt(1/length(gwl_historic_data))*1.149;
    sigma = 1.48/sqrt(2)*median(abs(gwl_historic_data(2:end) - gwl_historic_data(1:(end-1))))/5;
    gwl = TSregression(gwl_historic_data, sigma*h);
    
    G0_0    = gwl(1);
    gamma_0 = 0.00;
    mu_0    = mean(gwl);
    K_0     = zeros(tau+1, 1); %1E-2*ones(tau+1, 1);
    
    input0  = [G0_0; gamma_0; mu_0; K_0];
    
    lb = [0; 0; 0; zeros(tau+1, 1)];
    rb = [];
    
    options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000);
    par = fmincon(@objective, input0, [], [], [], [], lb, rb, @constraints, options);
    
    display(['Estimated parameters: ']);
    display(['G_0 = ', num2str(par(1))]);
    display(['gamma = ',  num2str(par(2))]);
    display(['mu = ',  num2str(par(3))]);
    for k = 1:(tau+1)
        display(['K(', num2str(k - tau - 1), ') = ', num2str(par(tau + k))]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2);
    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 28; ySize = 14;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position', [0 0 xSize*50 ySize*50]);
    
    plot_groundwater(par);
    
    function res = objective(input)
        G0_r    = gwl(tau);
        gamma_r = 0;
        mu_r    = mean(gwl);
        K_r     = zeros(tau+1, 1);
        
        G0    = input(1);
        gamma = input(2);
        mu    = input(3);
        K     = input(4:end);
        
        G = recover_G(G0, gamma, mu, K);
        
        res = sum([(G0 - G0_r)*sigma/sqrt(n); ...
                   (gamma - gamma_r)*sigma/sqrt(n); ...
                   (mu - mu_r)*sigma/sqrt(n); ...
                   (K(1) - K_r(1))*sigma/sqrt(n); ...
                   derivative(K - K_r)*sqrt(dt)*sigma/sqrt(n); ...
                   antiderivative(G - gwl((tau+1):end))*sqrt(dt)].^2);
    end

    function [c, ceq] = constraints(input)
        K  = input(4:end);
        DK = (K(2:end) - K(1:(end - 1)))/dt;
        
        c = [];
        c   = -DK;
        ceq = [];
    end

    function G = recover_G(G0, gamma, mu, K)
        G = zeros(n, 1);
        
        C = 1/(1 + gamma*dt);
                
        % implicit Euler scheme
        for i = 1:n
            mem_int = dt*sum(K.*rf_data(i:(i + tau)));
                        
            if (i == 1)
                G(i) = C*G0 + C*gamma*mu*dt + C*dt*mem_int;
            else
                G(i) = C*G(i - 1) + C*gamma*mu*dt + C*dt*mem_int;
            end
        end
    end

    function G = predict_G(G0, gamma, mu, K)
        G = zeros(m, 1);
        
        C = 1/(1 + gamma*dt);
        
        % implicit Euler scheme
        for i = 1:m
            mem_int = dt*sum(K.*rf_data((N_historic+i-1-tau):(N_historic+i-1)));
                        
            if (i == 1)
                G(i) = C*G0 + C*gamma*mu*dt + C*dt*mem_int;
            else
                G(i) = C*G(i - 1) + C*gamma*mu*dt + C*dt*mem_int;
            end
        end
    end

    function D = derivative(f)
        dt = time_data(2) - time_data(1);
        D = (f(2:end) - f(1:(end - 1)))/dt;
    end

    function I = antiderivative(f)
        I = cumsum(f)*dt;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_groundwater(input)
        G0    = input(1);
        gamma = input(2);
        mu    = input(3);
        K     = input(4:end);
        
        G = recover_G(G0, gamma, mu, K);
        
        c = 0.5*(min(gwl_data) + max(gwl_data));
        r = 0.6*(max(gwl_data) - min(gwl_data));
        
        axis([min(time_data) max(time_data) c - r c + r]);
        
        hold on;
        xlabel('Time $t$ in hrs', 'interpreter', 'latex', 'FontSize', 18);
        ylabel('Groundwater level $G(t)$ in m', 'interpreter', 'latex', 'FontSize', 18);
        
        % plot recovered groundwater curve
        p1 = plot(time_historic_data,  gwl_historic_data, 'o', 'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
        plot(time_historic_data,  gwl_historic_data, '--', 'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
        p2 = plot(time_historic_data, [gwl_historic_data(1:tau); G], '-', 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);      
        plot(time_historic_data(tau), gwl_historic_data(tau), '*', 'MarkerSize', 15, 'LineWidth', 1, 'Color', [0.8500 0.3250 0.0980]);
        
        plot(time_future_data, gwl_future_data, 'o',  'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
        plot(time_future_data, gwl_future_data, '--', 'LineWidth', 1, 'Color', [0.0000 0.4470 0.7410]);
        
        G = predict_G(gwl_historic_data(end), gamma, mu, K);
        
        % plot predicted groundwater curve
        p3 = plot([time_historic_data(end); time_future_data], [gwl_historic_data(end); G], 'g--', 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880]);
        plot(time_historic_data(end), gwl_historic_data(end), 'g*', 'MarkerSize', 15, 'LineWidth', 1, 'Color', [0.4660 0.6740 0.1880]); 
        
        legend([p1 p2 p3], {'Taut-string smoothed groundwater data', 'Historic groundwater level reconstruction based on estimated parameters', 'Future groundwater prediction based on estimated parameters'}, ...
                'Location', 'Best', 'interpreter', 'latex', 'FontSize', 13);
    end
end