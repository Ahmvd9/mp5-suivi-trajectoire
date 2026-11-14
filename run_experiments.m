%% run_experiments.m - Evaluation finale sur les graines 1 a nMC
% 7 filtres : 4 "de zero" (Kalman, Wiener, EKF, UKF) + 2 "informes" + brut.
% Chaque filtre recoit EXACTEMENT les memes donnees (meme sigma, meme graine).

clear; clc; close all;
params = get_params();

noms    = {'Brut', 'Kalman', 'Wiener', 'EKF', 'UKF', 'EKF info', 'UKF info'};
filtres = {@(z,p) z, @filtre_kf, @filtre_wiener, @filtre_ekf, @filtre_ukf, ...
    @filtre_ekf_info, @filtre_ukf_info};

nF = numel(filtres);  nS = numel(params.sigmas);  nMC = params.nMC;
rmse  = zeros(nF, nS, nMC);   tconv = zeros(nF, nS, nMC);
eperm = zeros(nF, nS, nMC);   tcalc = zeros(nF, nS, nMC);

% --- Echauffement : evite que le 1er appel fausse les temps de calcul ---
p0 = params;  p0.sigma_mes = params.sigmas(1);  p0.sigma = params.sigmas(1);
[~, z0] = simulate_trajectory(params.sigmas(1), 1);
for f = 1:nF, filtres{f}(z0, p0); end

% --- Boucle Monte Carlo ---
for j = 1:nS
    s = params.sigmas(j);
    p = params;  p.sigma_mes = s;  p.sigma = s;      % les deux noms de sigma
    for n = 1:nMC
        [x_true, z] = simulate_trajectory(s, n);
        for f = 1:nF
            tic;
            x_est = filtres{f}(z, p);
            tcalc(f, j, n) = toc;
            m = compute_metrics(x_true, x_est, s);
            rmse(f, j, n)  = m.rmse;
            tconv(f, j, n) = m.t_conv;
            eperm(f, j, n) = m.err_perm;
        end
    end
    fprintf('sigma = %g termine\n', s);
end

save('resultats_mp5.mat', 'rmse', 'tconv', 'eperm', 'tcalc', 'noms', 'params');

% --- Tableaux de synthese ---
vn = strcat('sigma_', string(params.sigmas));
T1 = array2table(mean(rmse, 3),  'VariableNames', vn, 'RowNames', noms);
T2 = array2table(mean(eperm, 3), 'VariableNames', vn, 'RowNames', noms);
T3 = array2table(mean(tconv, 3, 'omitnan'), 'VariableNames', vn, 'RowNames', noms);
T4 = array2table(max(rmse, [], 3), 'VariableNames', vn, 'RowNames', noms);
T5 = table(1000*mean(mean(tcalc, 3), 2), 'VariableNames', {'ms'}, 'RowNames', noms);
disp('RMSE moyen (m) :');                       disp(T1);
disp('Erreur en regime permanent moyenne (m) :'); disp(T2);
disp('Temps de convergence moyen (s), NaN ignores :'); disp(T3);
disp('RMSE pire cas, max sur les graines (m) :'); disp(T4);
disp('Temps de calcul moyen par appel (ms) :');  disp(T5);

writetable(T1, 'tab_rmse.csv',   'WriteRowNames', true);
writetable(T2, 'tab_eperm.csv',  'WriteRowNames', true);
writetable(T3, 'tab_tconv.csv',  'WriteRowNames', true);
writetable(T4, 'tab_pirecas.csv','WriteRowNames', true);
writetable(T5, 'tab_temps.csv',  'WriteRowNames', true);

% --- Styles : les filtres "info" en pointilles ---
sty = repmat({'-o'}, 1, nF);
for f = 1:nF, if contains(noms{f}, 'info'), sty{f} = '--s'; end, end

% --- Figure 1 : RMSE en fonction du bruit ---
h1 = figure; hold on; grid on;
for f = 1:nF
    errorbar(params.sigmas, mean(rmse(f,:,:), 3), std(rmse(f,:,:), 0, 3), sty{f}, 'LineWidth', 1.3);
end
set(gca, 'XScale', 'log', 'YScale', 'log');  xticks(params.sigmas);
xlabel('Bruit de mesure \sigma (m)');  ylabel('RMSE (m)');
legend(noms, 'Location', 'northwest');
title(sprintf('RMSE en fonction du bruit (%d realisations)', nMC));

% --- Figure 2 : erreur en regime permanent ---
h2 = figure; hold on; grid on;
for f = 1:nF
    errorbar(params.sigmas, mean(eperm(f,:,:), 3), std(eperm(f,:,:), 0, 3), sty{f}, 'LineWidth', 1.3);
end
set(gca, 'XScale', 'log', 'YScale', 'log');  xticks(params.sigmas);
xlabel('Bruit de mesure \sigma (m)');  ylabel('Erreur en regime permanent (m)');
legend(noms, 'Location', 'northwest');
title('Erreur en regime permanent');

% --- Figure 3 : temps de calcul ---
h3 = figure;
bar(1000 * mean(mean(tcalc, 3), 2));
set(gca, 'XTickLabel', noms);  ylabel('Temps de calcul par appel (ms)');  grid on;
title('Cout de calcul (Wiener inclut son apprentissage)');

exportgraphics(h1, 'fig_rmse.png',  'Resolution', 200);
exportgraphics(h2, 'fig_eperm.png', 'Resolution', 200);
exportgraphics(h3, 'fig_temps.png', 'Resolution', 200);