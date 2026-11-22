function m = compute_metrics(x_true, x_est, sigma)
%COMPUTE_METRICS  Metriques de performance d'un filtre de suivi.
%
%   m = compute_metrics(x_true, x_est, sigma)
%
%   Entrees
%     x_true : 4 x N (ou 2 x N) verite terrain ; seules les lignes 1:2
%              (positions x, y) sont utilisees
%     x_est  : 2 x N positions estimees par le filtre
%     sigma  : ecart-type du bruit de mesure (m), sert au seuil de convergence
%
%   Sortie : structure m avec
%     m.rmse     : RMSE de position sur toute la duree (m)
%     m.t_conv   : temps de convergence (s) = premier instant a partir duquel
%                  l'erreur LISSEE (moyenne glissante sur conv_window secondes)
%                  reste sous conv_factor*sigma jusqu'a la fin
%                  (NaN si ce n'est jamais le cas)
%     m.err_perm : erreur moyenne en regime permanent (m), sur les derniers
%                  perm_frac de la duree
%
%   Le temps de calcul se mesure a part, avec tic/toc autour du filtre.
%   Seuils et fenetres sont fixes dans get_params.m
%   (conv_factor, conv_window, perm_frac).

params = get_params();                        % charge la structure "params"

p_true = x_true(1:2, :);
N = size(p_true, 2);
assert(isequal(size(x_est), [2 N]), ...
    'x_est doit etre de taille 2 x N (N = %d).', N);

err = sqrt(sum((p_true - x_est).^2, 1));      % erreur de position a chaque instant (1 x N)

% RMSE global
m.rmse = sqrt(mean(err.^2));

% Temps de convergence, calcule sur l'erreur lissee (moyenne glissante) :
% un seul point isole au-dessus du seuil ne bloque pas la convergence.
w = max(1, round(params.conv_window / params.dt));  % taille de fenetre (echantillons)
err_s  = movmean(err, w);                     % erreur lissee (1 x N)
seuil  = params.conv_factor * sigma;
k_bad  = find(err_s > seuil, 1, 'last');      % dernier echantillon au-dessus du seuil
if isempty(k_bad)
    m.t_conv = 0;                             % toujours sous le seuil
elseif k_bad == N
    m.t_conv = NaN;                           % ne converge jamais
else
    m.t_conv = k_bad * params.dt;             % premier echantillon sous le seuil : k_bad+1 -> temps k_bad*dt
end

% Erreur en regime permanent : moyenne sur la fin de la simulation
n0 = N - round(params.perm_frac * N) + 1;
m.err_perm = mean(err(n0:N));
end
