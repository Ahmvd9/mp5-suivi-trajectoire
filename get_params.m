function params = get_params()
%GET_PARAMS  Parametres communs du MP5 (suivi de trajectoire bruitee).
%   Regle : ne PAS modifier seul. Toute modification se decide a trois,
%   car elle oblige a relancer toutes les experiences.

% --- Trajectoire (verite terrain) ---
params.R      = 100;            % rayon du cercle (m)
params.omega  = 0.1;            % vitesse angulaire (rad/s) -> vitesse lineaire R*omega = 10 m/s
params.dt     = 0.1;            % pas de temps du capteur (s)
params.T      = 60;             % duree de la simulation (s) -> N = T/dt = 600 echantillons

% --- Bruit de mesure (la difficulte) ---
params.sigmas = [1 2 5 10 20];  % ecarts-types du bruit gaussien (m), niveaux testes

% --- Protocole experimental ---
params.nMC         = 100;       % realisations Monte Carlo par sigma (graines 1 a nMC)
params.seed_tuning = 999;       % graine RESERVEE au reglage des filtres (q, P0, ...)
% Ne jamais l'utiliser pour l'evaluation finale.

% --- Metriques (utilisees par compute_metrics) ---
params.conv_factor = 2;         % seuil de convergence = conv_factor * sigma (m)
params.conv_window = 2;         % fenetre de lissage de l'erreur pour la convergence (s)
params.perm_frac   = 0.2;       % regime permanent = derniers 20 % de la duree

params.q = 1;               % bruit de processus du Kalman (regle sur la graine 999)
params.L_sigma = [10 20 20 50 50];   % ordre de Wiener par sigma (regle sur graines 990-998)
% --- Parametres propres aux filtres : a ajouter apres accord du groupe ---
