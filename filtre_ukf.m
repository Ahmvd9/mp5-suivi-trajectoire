function x_est = filtre_ukf(z, params)
%FILTRE_UKF  Filtre de Kalman "unscented" (UKF) DE ZERO.
%
%   Meme etat, meme P0 et meme Q que l'EKF de zero (comparaison equitable).
%   Idee de l'UKF : au lieu de LINEARISER le modele (Jacobienne), on choisit
%   2n+1 points ("points sigma") qui representent la gaussienne N(X,P), on les
%   fait passer TELS QUELS dans la fonction non lineaire, puis on recalcule
%   moyenne et covariance avec des poids. Aucune derivee n'est necessaire.

    % --- 1) Lecture de sigma (meme regle que l'EKF) ---
    if isfield(params, 'sigma_mes'),   sigma = params.sigma_mes;
    elseif isfield(params, 'sigma'),   sigma = params.sigma;
    else, error('Renseigner params.sigma_mes avant l''appel.');
    end
    dt = params.dt;  N = size(z, 2);
    n  = 5;                               % dimension de l'etat
    assert(size(z,1) == 2, 'z doit etre de taille 2 x N.');

    % --- 2) Bruit de processus (memes valeurs par defaut que l'EKF) ---
    qa = 1e-4;  qo = 1e-8;
    if isfield(params, 'qa_ukf'), qa = params.qa_ukf; end
    if isfield(params, 'qo_ukf'), qo = params.qo_ukf; end

    Rm = sigma^2 * eye(2);                % bruit de mesure
    Q1 = qa * [dt^3/3, dt^2/2; dt^2/2, dt];
    Q  = zeros(5);
    Q([1 3],[1 3]) = Q1;                  % axe x
    Q([2 4],[2 4]) = Q1;                  % axe y
    Q(5,5) = qo * dt;                     % omega

    % --- 3) Parametres de la transformation unscented ---
    alpha  = 1;                           % etalement des points sigma autour de la moyenne
    beta   = 2;                           % 2 = optimal pour une loi gaussienne
    kappa  = 0;                           % parametre secondaire
    lambda = alpha^2*(n + kappa) - n;     % facteur d'echelle (= 0 ici)
    % Poids : Wm pour la moyenne, Wc pour la covariance
    Wm = repmat(1/(2*(n + lambda)), 1, 2*n+1);
    Wc = Wm;
    Wm(1) = lambda/(n + lambda);                        % poids du point central (moyenne)
    Wc(1) = lambda/(n + lambda) + (1 - alpha^2 + beta); % poids du point central (covariance)

    % --- 4) Initialisation DEPUIS ZERO (identique a l'EKF) ---
    X = [z(1,1); z(2,1); 0; 0; 0];
    P = diag([sigma^2, sigma^2, 30^2, 30^2, 0.3^2]);

    x_est = zeros(2, N);  x_est(:,1) = X(1:2);

    % --- 5) Boucle de filtrage ---
    for k = 2:N
        % PREDICTION : points sigma de N(X,P) -> propages par f -> moyenne/cov
        [Xp, Pp] = transformation_unscented(@(s) f_etat(s, dt), X, P, Wm, Wc, lambda);
        Pp = Pp + Q;                      % + bruit de processus

        % CORRECTION : nouveaux points sigma autour de (Xp,Pp), passes dans h (mesure)
        [zp, S, Pxz] = transformation_unscented(@h_mesure, Xp, Pp, Wm, Wc, lambda);
        S = S + Rm;                       % + bruit de mesure
        K = Pxz / S;                      % gain de Kalman (covariance croisee / innovation)
        X = Xp + K*(z(:,k) - zp);         % correction par l'innovation
        P = Pp - K*S*K';                  % covariance a posteriori
        P = (P + P')/2;                   % symetrisation

        x_est(:,k) = X(1:2);
    end
end

function Xn = f_etat(X, dt)
%F_ETAT  Meme modele de virage que l'EKF (integration exacte sur dt).
    om = X(5);  c = cos(om*dt);  s = sin(om*dt);
    if abs(om) < 1e-6
        a = dt;  b = om*dt^2/2;           % limite omega -> 0
    else
        a = s/om;  b = (1-c)/om;
    end
    Xn = [ X(1) + a*X(3) - b*X(4);
           X(2) + b*X(3) + a*X(4);
           c*X(3) - s*X(4);
           s*X(3) + c*X(4);
           om ];
end

function zm = h_mesure(X)
%H_MESURE  Fonction de mesure : on observe uniquement la position (x, y).
    zm = X(1:2);
end

function chi = points_sigma(X, P, lambda)
%POINTS_SIGMA  Genere les 2n+1 points sigma (matrice n x (2n+1)).
%   Point central X, puis X +/- chaque colonne de L, ou L*L' = (n+lambda)*P
%   (L = "racine carree" de la covariance, par Cholesky).
    n = numel(X);
    L = cholesky_securisee((n + lambda)*P);
    chi = [X, X + L, X - L];
end

function [ym, Py, Pxy] = transformation_unscented(fun, X, P, Wm, Wc, lambda)
%TRANSFORMATION_UNSCENTED  Propage (X,P) a travers la fonction fun.
%   ym : moyenne de fun(X) ; Py : covariance de fun(X) ;
%   Pxy : covariance croisee entre X et fun(X).
    chi = points_sigma(X, P, lambda);     % 1) points sigma
    nb  = size(chi, 2);
    m   = numel(fun(chi(:,1)));           %    dimension de la sortie
    Y   = zeros(m, nb);
    for i = 1:nb
        Y(:,i) = fun(chi(:,i));           % 2) propagation NON lineaire de chaque point
    end
    ym  = Y*Wm(:);                        % 3) moyenne ponderee
    dY  = Y - ym;                         %    ecarts a la moyenne (sortie)
    dX  = chi - X;                        %    ecarts a la moyenne (entree)
    Py  = dY*diag(Wc)*dY';                % 4) covariances ponderees
    Pxy = dX*diag(Wc)*dY';
end

function L = cholesky_securisee(A)
%CHOLESKY_SECURISEE  Racine carree L (L*L' = A) avec regularisation.
%   Si A n'est pas definie positive (erreurs d'arrondi), on ajoute une petite
%   diagonale "jitter" croissante ; en dernier recours, valeurs propres ecretees.
    A = (A + A')/2;  n = size(A, 1);
    [L, p] = chol(A, 'lower');            % p = 0 si A est definie positive
    if p == 0, return; end
    jitter = 1e-9*max(1, trace(A)/n);
    for essai = 1:10
        [L, p] = chol(A + jitter*eye(n), 'lower');
        if p == 0, return; end
        jitter = 10*jitter;
    end
    [V, D] = eig(A);                      % repli ultime
    L = V*diag(sqrt(max(diag(D), 1e-12)));
end