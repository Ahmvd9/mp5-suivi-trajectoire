function x_est = filtre_wiener(z, params)
%FILTRE_WIENER  Filtre de Wiener FIR causal, applique a x et y separement.
%
%   x_est = filtre_wiener(z, params)
%
%   Entrees
%     z      : 2 x N mesures bruitees
%     params : structure issue de get_params(), avec en plus :
%                params.sigma_mes : ecart-type du bruit (m)          [OBLIGATOIRE]
%                params.L         : ordre du filtre FIR              [optionnel]
%                params.L_sigma   : un L par valeur de params.sigmas [optionnel]
%              (sans L ni L_sigma, on prend L = 20)
%
%   Sortie : x_est, 2 x N positions estimees
%
%   Principe : les poids h minimisent l'erreur quadratique moyenne.
%   Equations de Wiener-Hopf : toeplitz(r_zz) * h = r_zs.
%   Les correlations sont estimees sur une simulation d'APPRENTISSAGE
%   (graine params.seed_tuning), ou la verite est connue.

if ~isfield(params, 'sigma_mes')
    error('filtre_wiener : renseigner params.sigma_mes avant l''appel.');
end
sig = params.sigma_mes;

% --- Choix de l'ordre L ---
L = 20;
if isfield(params, 'L')
    L = params.L;
elseif isfield(params, 'L_sigma')
    idx = find(params.sigmas == sig, 1);
    if ~isempty(idx), L = params.L_sigma(idx); end
end

% --- 1) Simulation d'apprentissage : on connait la verite ---
[x_tr, z_tr] = simulate_trajectory(sig, params.seed_tuning);
s_tr = x_tr(1:2, :);                         % signal vrai (positions)
N_tr = size(z_tr, 2);

% --- 2) Correlations estimees (moyenne sur le temps, x et y regroupes) ---
r_zz = zeros(L, 1);                          % autocorrelation de la mesure
r_zs = zeros(L, 1);                          % correlation mesure / signal vrai
for tau = 0:L-1
    k = (tau+1):N_tr;
    A = z_tr(:, k)     .* z_tr(:, k - tau);
    B = z_tr(:, k-tau) .* s_tr(:, k);
    r_zz(tau+1) = mean(A(:));
    r_zs(tau+1) = mean(B(:));
end

% --- 3) Equations de Wiener-Hopf ---
h = toeplitz(r_zz) \ r_zs;                   % poids du filtre (L x 1)
h = h(:).';                                  % en ligne (1 x L)

% --- 4) Application au signal : filtrage de chaque axe ---
N = size(z, 2);
x_est = zeros(2, N);
for c = 1:2
    z_pad = [z(c,1) * ones(1, L-1), z(c, :)];   % debut complete par la 1ere mesure
    y = filter(h, 1, z_pad);
    x_est(c, :) = y(L:end);
end
end