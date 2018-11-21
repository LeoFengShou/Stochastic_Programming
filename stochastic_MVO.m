function x_optimal = stochastic_MVO(mu, Q, targetRet)
    num_asset = size(mu, 1);
    rho = corrcov(Q);
    nPaths = 200;
    L = chol(rho, 'lower');
    T = 26;
    N = 1;
    dt = T/N;
    reward_per_dollor_surplus = -2; 
    % Because it is in a minimization problem, so minus means reward
    
    punishment_per_dollor_surplus = 1;
    % Because it is in a minimization problem, so positive means punishment
    
    risk_weight_coefficient = 400;
    
    S = zeros(num_asset, N+1, nPaths);  % Matrix of simulated price paths
    S(:, 1, :) = 100;
    % Generate paths
    for i = 1:nPaths
        for j = 1:N
            xi = L * randn(num_asset, 1); 
            for k = 1:num_asset
                S(k, j+1, i) = S(k, j, i) * exp( ( mu(k, 1) - 0.5 * Q (k, k) ) * dt ...
                                + sqrt(Q(k, k)) * sqrt(dt) * xi (k) );
            end
        end
    end 
    returns_sample = zeros(num_asset, nPaths); % returns_sample 20 * 2000
    for i = 1:nPaths
       for j = 1:num_asset
            returns_sample(j, i) = S(j, end, i) / S(j, 1, i) - 1;
       end
    end
    % Have 2000 scenarios and each scenatio has a surplus and shortfall
    % Also has num_asset = 20 weight variables
    % In total, have 2000 * 2 + 20 = 4020 variables.
    % First 2000 variables are surplus variables.
    % Second 2000 are shortfall variables.
    % The last 20 are weight variables.
    f = zeros(2 * nPaths + num_asset, 1);
    f(1: nPaths) = 1 / nPaths * reward_per_dollor_surplus;
    f(nPaths + 1: 2 * nPaths) = 1 / nPaths * punishment_per_dollor_surplus;
    Aeq = zeros(nPaths + 1, 2 * nPaths + num_asset);
    Aeq(1:nPaths, 1:nPaths) = -1 * eye(nPaths);
    Aeq(1:nPaths, nPaths + 1: 2 *nPaths) = 1 * eye(nPaths);
    for i = 1:nPaths
        Aeq(i, 2 * nPaths + 1: 2 * nPaths + num_asset) = (returns_sample(:, i))';
    end
    Aeq(nPaths + 1, 2 * nPaths + 1: 2 * nPaths + num_asset) = 1;
    beq = ones(nPaths + 1, 1);
    beq(1: nPaths, 1) = targetRet;
    Q_to_q = [zeros(2 * nPaths, 2 * nPaths), zeros(2 * nPaths, num_asset); ...
              zeros(num_asset, 2 * nPaths), risk_weight_coefficient * Q];
    lb = zeros(2 * nPaths + num_asset, 1);
    ub = ones(2 * nPaths + num_asset, 1);
    % ub(1: 2 * nPaths, 1) = Inf;
%    lb(2 * nPaths + 1: 2 * nPaths + num_asset) = -Inf;
    x = quadprog(Q_to_q, f, [], [], Aeq, beq, lb, ub);
    x_optimal = x(2 * nPaths + 1: 2 * nPaths + num_asset, 1);
end