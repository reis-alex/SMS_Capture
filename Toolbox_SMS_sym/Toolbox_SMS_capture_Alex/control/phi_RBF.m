function phi = phi_RBF(w,centers,widths)
    %   w: The input state vector (Nx1) where N is the number of states
    %   centers     : The centers of the RBF functions (KxN) where K is the number of centers
    %   widths      : The widths of the RBF functions (1xK)
    
    %[K, N] = size(centers);
    % distanceSquared = sum((w - centers).^2, 2);
    % rbfOutput = exp(-distanceSquared ./ (2 * widths.^2));
    N = length(w); % Number of states (dimension of the state space)
    K = size(centers, 1); % Number of RBF centers
    expandedStateVector = repmat(w', K, 1);
    distanceSquared = sum((expandedStateVector - centers).^2, 2);
    rbfOutput = exp(-distanceSquared ./ (2 * widths.^2));
    

phi = rbfOutput;
