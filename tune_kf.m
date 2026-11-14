%% tune_kf.m  -  Reglage du parametre q du Kalman sur la graine de reglage 999
% Principe : on essaie plusieurs valeurs de q et on garde celle qui minimise
% le RMSE. On utilise UNIQUEMENT la graine 999, jamais les graines 1 a 100.

clear; clc;
params = get_params();

q_grid = [1e-4 1e-3 1e-2 1e-1 1 10 100];     % valeurs de q a essayer
sigmas = params.sigmas;
rmse   = zeros(numel(q_grid), numel(sigmas));

for j = 1:numel(sigmas)
    s = sigmas(j);
    [x_true, z] = simulate_trajectory(s, params.seed_tuning);
    p = params;  p.sigma_mes = s;
    for i = 1:numel(q_grid)
        p.q = q_grid(i);
        x_est = filtre_kf(z, p);
        m = compute_metrics(x_true, x_est, s);
        rmse(i, j) = m.rmse;
    end
end

% Affichage : une ligne par q, une colonne par sigma
T = array2table(rmse, 'VariableNames', ...
    strcat('sigma_', string(sigmas)), 'RowNames', strcat('q=', string(q_grid)));
disp(T);

[~, i_best] = min(mean(rmse, 2));
fprintf('q qui minimise le RMSE moyen sur tous les sigma : %g\n', q_grid(i_best));
fprintf('A reporter dans get_params.m : params.q = %g;\n', q_grid(i_best));
