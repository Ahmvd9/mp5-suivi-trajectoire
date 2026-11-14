function x_est = filtre_kf(z, params)
%FILTRE_KF  Filtre de Kalman lineaire a vitesse constante (2D).
%
%   x_est = filtre_kf(z, params)
%
%   Entrees
%     z      : 2 x N mesures bruitees de position (ligne 1 = x, ligne 2 = y)
%     params : structure issue de get_params(), avec en plus :
%                params.sigma_mes : ecart-type du bruit de mesure (m)  [OBLIGATOIRE]
%                params.q         : bruit de processus (m^2/s^3)       [optionnel, defaut 0.1]
%                params.P0_vit    : variance initiale des vitesses     [optionnel, defaut 100]
%
%   Sortie
%     x_est  : 2 x N positions estimees (ligne 1 = x, ligne 2 = y)
%
%   Modele : etat [x; vx; y; vy], la cible garde une vitesse constante, aux
%   perturbations pres (bruit d'acceleration, regle par q). Sur un cercle ce
%   modele est FAUX (la vitesse change de direction) : c'est le defaut du KF
%   lineaire que l'EKF et l'UKF corrigent.
%
%   Prediction : x = F*x          P = F*P*F' + Q
%   Correction : K = P*H'/(H*P*H' + R)
%                x = x + K*(z - H*x)      P = (I - K*H)*P

    % --- Parametres ---
    if ~isfield(params, 'sigma_mes')
        error('filtre_kf : renseigner params.sigma_mes avant l''appel.');
    end
    sig = params.sigma_mes;
    if isfield(params, 'q'),      q = params.q;           else, q = 0.1;  end
    if isfield(params, 'P0_vit'), P0v = params.P0_vit;    else, P0v = 100; end
    dt = params.dt;
    N  = size(z, 2);

    % --- Matrices du modele ---
    F1 = [1 dt; 0 1];                        % x_{k+1} = x_k + dt*vx_k
    F  = blkdiag(F1, F1);                    % etat [x vx y vy]
    H  = [1 0 0 0; 0 0 1 0];                 % on mesure x et y seulement
    Q1 = q * [dt^3/3, dt^2/2; dt^2/2, dt];   % bruit d'acceleration (blanc continu), par axe
    Q  = blkdiag(Q1, Q1);
    Rm = sig^2 * eye(2);                     % bruit de mesure

    % --- Initialisation : position = premiere mesure, vitesse inconnue ---
    x = [z(1,1); 0; z(2,1); 0];
    P = diag([sig^2, P0v, sig^2, P0v]);

    x_est = zeros(2, N);
    x_est(:, 1) = z(:, 1);

    % --- Boucle : prediction puis correction ---
    for k = 2:N
        x = F * x;                           % prediction de l'etat
        P = F * P * F' + Q;                  % prediction de l'incertitude

        S = H * P * H' + Rm;                 % incertitude sur la mesure predite
        K = (P * H') / S;                    % gain de Kalman
        x = x + K * (z(:, k) - H * x);       % correction par l'innovation
        P = (eye(4) - K * H) * P;

        x_est(:, k) = H * x;                 % position estimee
    end
end
