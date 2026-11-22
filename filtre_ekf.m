function x_est = filtre_ekf(z, params)
%FILTRE_EKF  Filtre de Kalman ETENDU (EKF) DE ZERO, modele de virage.
%
%   Le filtre ne connait que : la 1ere mesure et le bruit sigma.
%   Il ne connait NI le rayon, NI omega, NI le centre du cercle.
%
%   Etat a 5 composantes : X = [ x ; y ; vx ; vy ; omega ]
%     (x,y)    position           (vx,vy)  vitesse en coordonnees cartesiennes
%     omega    vitesse de rotation de la vitesse (rad/s), a DECOUVRIR
%
%   Pourquoi "etendu" : le modele de virage est NON lineaire (cos, sin, produit
%   omega*v). On le LINEARISE a chaque pas autour de l'estimation courante,
%   grace a la Jacobienne F.

    if isfield(params, 'sigma_mes'),   sigma = params.sigma_mes;
    elseif isfield(params, 'sigma'),   sigma = params.sigma;
    else, error('Renseigner params.sigma_mes avant l''appel.');
    end
    dt = params.dt;                       % pas de temps (s)
    N  = size(z, 2);                      % nombre d'echantillons
    assert(size(z,1) == 2, 'z doit etre de taille 2 x N.');

    % --- 2) Reglages du bruit de processus (regles sur la graine 999) ---
    qa = 1e-2;                            % incertitude sur l'acceleration (m^2/s^3)
    qo = 1e-6;                            % incertitude sur la variation de omega
    if isfield(params, 'qa_ekf'), qa = params.qa_ekf; end   % surcharge possible
    if isfield(params, 'qo_ekf'), qo = params.qo_ekf; end   % (utilisee par tune)

    % --- 3) Matrices du modele ---
    Rm = sigma^2 * eye(2);                % bruit de mesure : x et y independants
    H  = [1 0 0 0 0;                      % on ne mesure que la position (x, y)
          0 1 0 0 0];
    Q1 = qa * [dt^3/3, dt^2/2;            % bruit d'acceleration (forme classique
               dt^2/2, dt    ];           %  position/vitesse), applique a chaque axe
    Q  = zeros(5);
    Q([1 3],[1 3]) = Q1;                  % axe x : composantes (x, vx)
    Q([2 4],[2 4]) = Q1;                  % axe y : composantes (y, vy)
    Q(5,5) = qo * dt;                     % omega peut deriver tres lentement

    % --- 4) Initialisation DEPUIS ZERO ---
    X = [z(1,1); z(2,1); 0; 0; 0];        % position = 1ere mesure ; vitesse=0 ; omega=0
    % P0 : incertitude initiale. Position : sigma (la mesure est bruitee de sigma).
    % Vitesse : tres large (30 m/s) et omega large (0.3 rad/s) : le filtre n'en
    % sait rien et doit se laisser corriger par les mesures.
    P = diag([sigma^2, sigma^2, 30^2, 30^2, 0.3^2]);

    x_est = zeros(2, N);                  % sortie : positions estimees
    x_est(:,1) = X(1:2);

    % --- 5) Boucle : prediction puis correction ---
    for k = 2:N
        % PREDICTION
        F = jacobienne_f(X, dt);          % linearisation du modele en X (5x5)
        X = f_etat(X, dt);                % l'etat avance via le modele NON lineaire
        P = F*P*F' + Q;                   % l'incertitude grandit (via F) + bruit Q

        % CORRECTION
        S = H*P*H' + Rm;                  % incertitude sur la mesure predite
        K = (P*H') / S;                   % gain de Kalman : confiance mesure vs modele
        X = X + K * (z(:,k) - H*X);       % on corrige par l'innovation (mesure - prevision)
        IKH = eye(5) - K*H;
        P = IKH*P*IKH' + K*Rm*K';         % forme de Joseph (reste definie positive)
        P = (P + P')/2;                   % symetrisation (securite numerique)

        x_est(:,k) = X(1:2);              % on ne renvoie que (x, y)
    end
end

function Xn = f_etat(X, dt)
%F_ETAT  Modele d'evolution : virage a omega constant (integration EXACTE sur dt).
%   La vitesse tourne de l'angle omega*dt, la position avance avec la vitesse moyenne.
    om = X(5);  c = cos(om*dt);  s = sin(om*dt);
    if abs(om) < 1e-6
        a = dt;  b = om*dt^2/2;           % limite omega -> 0 : mouvement rectiligne
    else
        a = s/om;  b = (1-c)/om;          % coefficients exacts de la rotation
    end
    Xn = [ X(1) + a*X(3) - b*X(4);        % x
           X(2) + b*X(3) + a*X(4);        % y
           c*X(3) - s*X(4);               % vx (rotation du vecteur vitesse)
           s*X(3) + c*X(4);               % vy
           om ];                          % omega suppose constant
end

function F = jacobienne_f(X, dt)
%JACOBIENNE_F  Derivees partielles de f_etat par rapport a chaque composante de X
%   (matrice 5x5). Ligne i, colonne j = derivee de la i-eme sortie / X(j).
    vx = X(3);  vy = X(4);  om = X(5);
    c = cos(om*dt);  s = sin(om*dt);
    if abs(om) < 1e-6
        a = dt;  b = om*dt^2/2;  da = -om*dt^3/3;  db = dt^2/2;
    else
        a = s/om;  b = (1-c)/om;
        da = (dt*c*om - s)/om^2;          % d(a)/d(omega)
        db = (dt*s*om - (1-c))/om^2;      % d(b)/d(omega)
    end
    F = eye(5);                           % par defaut : chaque etat se derive en lui-meme
    F(1,3) = a;  F(1,4) = -b;  F(1,5) = da*vx - db*vy;   % d(x+)/d(vx, vy, omega)
    F(2,3) = b;  F(2,4) = a;   F(2,5) = db*vx + da*vy;   % d(y+)/d(...)
    F(3,3) = c;  F(3,4) = -s;  F(3,5) = dt*(-s*vx - c*vy);
    F(4,3) = s;  F(4,4) = c;   F(4,5) = dt*(c*vx - s*vy);
end
