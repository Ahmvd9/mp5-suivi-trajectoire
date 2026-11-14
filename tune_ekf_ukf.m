%% tune_ekf_ukf.m - regle qa et qo de l'EKF et de l'UKF de zero, graine 999
clear; clc;
params = get_params();
qa_grid = [1e-4 1e-3 1e-2 1e-1 1];
qo_grid = [1e-8 1e-6 1e-4];
rel = zeros(numel(qa_grid), numel(qo_grid), 2);       % 1 = EKF, 2 = UKF ; RMSE/sigma moyen
nS = numel(params.sigmas);

for j = 1:nS
    s = params.sigmas(j);
    [x_true, z] = simulate_trajectory(s, params.seed_tuning);
    p = params;  p.sigma_mes = s;
    for a = 1:numel(qa_grid)
        for o = 1:numel(qo_grid)
            p.qa_ekf = qa_grid(a);  p.qo_ekf = qo_grid(o);
            p.qa_ukf = qa_grid(a);  p.qo_ukf = qo_grid(o);
            m1 = compute_metrics(x_true, filtre_ekf(z, p), s);
            m2 = compute_metrics(x_true, filtre_ukf(z, p), s);
            rel(a,o,1) = rel(a,o,1) + m1.rmse/s/nS;
            rel(a,o,2) = rel(a,o,2) + m2.rmse/s/nS;
        end
    end
end

noms = {'EKF', 'UKF'};
for f = 1:2
    fprintf('--- %s (RMSE/sigma moyen) ---\n', noms{f});
    disp(array2table(rel(:,:,f), 'VariableNames', strcat('qo_', string(qo_grid)), ...
        'RowNames', strcat('qa=', string(qa_grid))));
    [~, idx] = min(reshape(rel(:,:,f), [], 1));
    [ia, io] = ind2sub([numel(qa_grid) numel(qo_grid)], idx);
    fprintf('Meilleur %s : qa = %g ; qo = %g\n\n', noms{f}, qa_grid(ia), qo_grid(io));
end