%% run_kf_demo.m  -  Premier test du filtre de Kalman (graine de reglage 999)
% Lancer dans le dossier qui contient get_params.m, simulate_trajectory.m,
% compute_metrics.m et filtre_kf.m.

clear; clc; close all;
params = get_params();

sigma = 5;                                   % niveau de bruit teste (m)
params.sigma_mes = sigma;                    % le filtre doit connaitre le bruit du capteur

% Graine de REGLAGE (999) : on ne touche pas aux graines 1 a 100 ici
[x_true, z] = simulate_trajectory(sigma, params.seed_tuning);

tic;  x_est = filtre_kf(z, params);  t_calc = toc;

m_kf   = compute_metrics(x_true, x_est, sigma);
m_brut = compute_metrics(x_true, z,     sigma);

fprintf('--- sigma = %g m, graine %d ---\n', sigma, params.seed_tuning);
fprintf('Mesures brutes : RMSE = %.2f m | t_conv = %.1f s | err_perm = %.2f m\n', ...
        m_brut.rmse, m_brut.t_conv, m_brut.err_perm);
fprintf('Kalman (KF)    : RMSE = %.2f m | t_conv = %.1f s | err_perm = %.2f m | calcul = %.4f s\n', ...
        m_kf.rmse, m_kf.t_conv, m_kf.err_perm, t_calc);

% --- Figure 1 : trajectoire ---
figure;
plot(x_true(1,:), x_true(2,:), 'k-', 'LineWidth', 1.5); hold on;
plot(z(1,:), z(2,:), '.', 'Color', [0.7 0.7 0.7]);
plot(x_est(1,:), x_est(2,:), 'r-', 'LineWidth', 1.5);
axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)');
legend('Verite', 'Mesures bruitees', 'Estimation KF');
title(sprintf('Suivi par filtre de Kalman (sigma = %g m)', sigma));

% --- Figure 2 : erreur en fonction du temps ---
t = (0:size(z,2)-1) * params.dt;
err_brut = sqrt(sum((x_true(1:2,:) - z).^2, 1));
err_kf   = sqrt(sum((x_true(1:2,:) - x_est).^2, 1));
figure;
plot(t, err_brut, '-', 'Color', [0.7 0.7 0.7]); hold on;
plot(t, err_kf, 'r-', 'LineWidth', 1.2);
yline(2*sigma, 'k--');
grid on; xlabel('Temps (s)'); ylabel('Erreur de position (m)');
legend('Mesures brutes', 'KF', 'Seuil 2*sigma');
title('Erreur en fonction du temps');
