%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      (C) Michael Pokojovy (2020)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function yr = TSregression(y, lambda)
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
end