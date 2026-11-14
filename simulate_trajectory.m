function [x_true, z] = simulate_trajectory(sigma, seed)
%SIMULATE_TRAJECTORY  Cible en mouvement circulaire + mesures bruitees.
%
%   [x_true, z] = simulate_trajectory(sigma, seed)
%
%   Entrees
%     sigma : ecart-type du bruit de mesure (m)
%     seed  : graine aleatoire (entier) -> resultats reproductibles
%   Sorties
%     x_true : 4 x N  verite terrain [x; y; vx; vy]
%              x  = R*cos(omega*t),       y  = R*sin(omega*t)
%              vx = -R*omega*sin(omega*t), vy = R*omega*cos(omega*t)
%     z      : 2 x N  positions mesurees = positions vraies + sigma*randn
%              (ligne 1 = x, ligne 2 = y)
%
%   Les parametres (R, omega, dt, T) viennent de params.m.

    params = get_params();                   % charge la structure "params"
    rng(seed);                               % reproductibilite

    N = round(params.T / params.dt);         % nombre d'echantillons
    t = (0:N-1) * params.dt;                 % vecteur temps (1 x N)
    R = params.R;  w = params.omega;

    x_true = [ R*cos(w*t);                   % x
               R*sin(w*t);                   % y
              -R*w*sin(w*t);                 % vx
               R*w*cos(w*t) ];               % vy

    z = x_true(1:2, :) + sigma * randn(2, N);  % mesures bruitees
end
