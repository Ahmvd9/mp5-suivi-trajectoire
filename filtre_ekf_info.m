function x_est = filtre_ekf_info(z, params)
%FILTRE_EKF_INFO , version INFORMEE (connait R, omega, cap, centre).
%
%   Etat a 5 composantes : X = [ x ; y ; v ; theta ; omega ]
%     v = vitesse lineaire, theta = cap (direction de la vitesse), omega = virage.
%   Modele (non lineaire : cos, sin, produit v*cos) :
%     x+ = x + v*cos(theta)*dt        y+ = y + v*sin(theta)*dt
%     v+ = v                          theta+ = theta + omega*dt      omega+ = omega

dt = params.dt;
assert(size(z,1) == 2, 'z doit etre de taille 2 x N.');
N  = size(z, 2);

% Lecture de sigma (sigma_mes de preference, sinon sigma, sinon erreur)
if isfield(params, 'sigma_mes'),   sigma = params.sigma_mes;
elseif isfield(params, 'sigma'),   sigma = params.sigma;
else, error('Renseigner params.sigma_mes avant l''appel.');
end
% omega initial suppose : params.omega0 si fourni, sinon la VRAIE valeur
if isfield(params, 'omega0'), omega0 = params.omega0; else, omega0 = params.omega; end

% Bruit de mesure et matrice d'observation
Rm = sigma^2 * eye(2);
H  = [1 0 0 0 0; 0 1 0 0 0];          % on ne mesure que (x, y)

% Bruit de processus Q (diagonale, petit : le mouvement est un cercle parfait)
q_pos = 1e-3;  q_v = 1e-4;  q_theta = 1e-6;  q_omega = 1e-7;
Q = diag([q_pos, q_pos, q_v, q_theta, q_omega]) * dt;

% INITIALISATION INFORMEE : c'est ici que le filtre "triche" (connait le cercle)
%   vitesse v0 = R*omega0 ; cap = tangent au cercle autour de l'ORIGINE
X = [ z(1,1);
    z(2,1);
    params.R * omega0;
    atan2(z(2,1), z(1,1)) + pi/2;
    omega0 ];
P = diag([sigma^2, sigma^2, 5^2, 0.5^2, 0.05^2]);   % incertitudes initiales

x_est = zeros(2, N);  x_est(:,1) = X(1:2);

for k = 2:N
    % PREDICTION
    F = jacobienne_f(X, dt);          % linearisation du modele en X
    X = f_etat(X, dt);                % propagation par le modele non lineaire
    P = F*P*F' + Q;

    % CORRECTION
    y_inn = z(:,k) - H*X;             % innovation (mesure - mesure predite)
    S = H*P*H' + Rm;                  % covariance de l'innovation
    K = (P*H') / S;                   % gain de Kalman
    X = X + K*y_inn;
    IKH = eye(5) - K*H;
    P = IKH*P*IKH' + K*Rm*K';         % forme de Joseph
    P = (P + P')/2;
    x_est(:,k) = X(1:2);
end
end

function Xn = f_etat(X, dt)
%F_ETAT  Evolution non lineaire (virage a vitesse constante), version Euler.
v = X(3);  theta = X(4);  omega = X(5);
Xn = [ X(1) + v*cos(theta)*dt;
    X(2) + v*sin(theta)*dt;
    v;
    theta + omega*dt;
    omega ];
end

function F = jacobienne_f(X, dt)
%JACOBIENNE_F  Jacobienne analytique df/dX (5x5), evaluee en X.
v = X(3);  theta = X(4);
F = eye(5);
F(1,3) =  cos(theta)*dt;   F(1,4) = -v*sin(theta)*dt;   % d(x+)/dv, d(x+)/dtheta
F(2,3) =  sin(theta)*dt;   F(2,4) =  v*cos(theta)*dt;   % d(y+)/dv, d(y+)/dtheta
F(4,5) =  dt;                                           % d(theta+)/domega
end