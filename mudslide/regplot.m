%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      (C) Michael Pokojovy (2020)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yr = regplot(y, lambda)
    set(gcf, 'PaperUnits', 'centimeters');
    xSize = 28; ySize = 14;
    xLeft = (21 - xSize)/2; yTop = (30 - ySize)/2;
    set(gcf,'PaperPosition', [xLeft yTop xSize ySize]);
    set(gcf,'Position',[0 0 xSize*50 ySize*50]);

    N = size(y, 1);
    
    x = linspace(1, N, N)'/N;
    x0 = [0; x];
    
    cy = [0; 1/N*cumsum(y)];
    
    cyl = cy - lambda; 
    cyl(1) = 0; 
    cyl(end) = cy(end);
    
    cyu = cy + lambda; 
    cyu(1) = 0; 
    cyu(end) = cy(end);
    
    [index, cyr] = TautString(x0, cyl, cyu);
    
    yr = zeros(N, 1);
    
    for i = 1:length(index)-1
        I = index(i):min(index(i+1), N);
        yr(I) = (cyr(i+1) - cyr(i))/(x0(index(i+1)) - x0(index(i)));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot_tight(1, 2, 1, [0.08 0.05]);
    hold on;
    xlabel('$t$', 'interpreter', 'latex', 'FontSize', 18);
    plot(x, y, 'bo');
	plot(x, yr, 'r');
    
    legend({'Raw data', 'Estimated signal $\hat{\mu}_{n}$'}, ...
            'Location', 'SouthWest', 'interpreter', 'latex', 'FontSize', 18);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot_tight(1, 2, 2, [0.08 0.05]);
    hold on;
    xlabel('$t$', 'interpreter', 'latex', 'FontSize', 18);
    
    X = []; Y = [];
	for i = 1:length(index)-1
        X = [X x0(index(i)) x0(index(i+1))];
        Y = [Y cyr(i) cyr(i+1)];
    end
    plot(X, Y, 'r');
    
    %plot([0; x], cy, 'b.');
    plot([0; x], cy - lambda, 'k');
    plot([0; x], cy + lambda, 'k');
    plot([0 1], [cyl(1) cyl(end)], 'k*');
    
    legend({'Taut string $s_{n}^{\ast}$', 'Functional tube $\mathcal{K}^{\theta_{n}}_{\varphi_{n}, \psi_{n}}$'}, ...
           'Location', 'SouthWest', 'interpreter', 'latex', 'FontSize', 18);
end