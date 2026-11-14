%% tune_wiener.m - Reglage de L (ordre de Wiener) pour chaque sigma
% On utilise des graines de VALIDATION (990 a 998), differentes de la graine
% d'apprentissage (999) et des graines d'evaluation (1 a 100).

clear; clc;
params = get_params();

seeds_val = 990:998;
L_grid    = [5 10 20 50 100];
sigmas    = params.sigmas;
rmse      = zeros(numel(L_grid), numel(sigmas));

for j = 1:numel(sigmas)
    s = sigmas(j);
    p = params;  p.sigma_mes = s;
    for i = 1:numel(L_grid)
        p.L = L_grid(i);
        r = zeros(1, numel(seeds_val));
        for n = 1:numel(seeds_val)
            [x_true, z] = simulate_trajectory(s, seeds_val(n));
            x_est = filtre_wiener(z, p);
            m = compute_metrics(x_true, x_est, s);
            r(n) = m.rmse;
        end
        rmse(i, j) = mean(r);
    end
end

T = array2table(rmse, 'VariableNames', ...
    strcat('sigma_', string(sigmas)), 'RowNames', strcat('L=', string(L_grid)));
disp(T);

[~, i_best] = min(rmse, [], 1);
fprintf('Meilleur L pour chaque sigma (%s) : %s\n', num2str(sigmas), num2str(L_grid(i_best)));
fprintf('A reporter dans get_params.m : params.L_sigma = [%s];\n', num2str(L_grid(i_best)));