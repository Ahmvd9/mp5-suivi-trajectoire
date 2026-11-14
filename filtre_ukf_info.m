function x_est = filtre_ukf_info(z, params)
%FILTRE_UKF_INFO  , version INFORMEE. Meme etat que l'EKF info :
%   X = [ x ; y ; v ; theta ; omega ]. Memes initialisations "informees".

dt = params.dt;
assert(size(z,1) == 2, 'z doit etre de taille 2 x N.');
N  = size(z, 2);  n = 5;

if isfield(params, 'sigma_mes'),   sigma = params.sigma_mes;
elseif isfield(params, 'sigma'),   sigma = params.sigma;
else, error('Renseigner params.sigma_mes avant l''appel.');
end
if isfield(params, 'omega0'), omega0 = params.omega0; else, omega0 = params.omega; end

% Reglages identiques a l'EKF info (comparaison equitable)
Rm = sigma^2 * eye(2);
q_pos = 1e-3;  q_v = 1e-4;  q_theta = 1e-6;  q_omega = 1e-7;
Q  = diag([q_pos, q_pos, q_v, q_theta, q_omega]) * dt;

% Transformation unscented : etalement, poids
alpha = 1;  beta = 2;  kappa = 0;
lambda = alpha^2*(n + kappa) - n;
Wm = repmat(1/(2*(n + lambda)), 1, 2*n+1);  Wc = Wm;
Wm(1) = lambda/(n + lambda);
Wc(1) = lambda/(n + lambda) + (1 - alpha^2 + beta);

% Initialisation INFORMEE (v0 = R*omega0, cap tangent au cercle)
X = [ z(1,1); z(2,1); params.R*omega0; atan2(z(2,1), z(1,1)) + pi/2; omega0 ];
P = diag([sigma^2, sigma^2, 5^2, 0.5^2, 0.05^2]);

x_est = zeros(2, N);  x_est(:,1) = X(1:2);

for k = 2:N
    % PREDICTION : transformation unscented a travers f_etat
    [Xp, Pp] = transformation_unscented(@(s) f_etat(s, dt), X, P, Wm, Wc, lambda);
    Pp = Pp + Q;
    % CORRECTION : transformation unscented a travers h_mesure
    [zp, S, Pxz] = transformation_unscented(@h_mesure, Xp, Pp, Wm, Wc, lambda);
    S = S + Rm;
    K = Pxz / S;
    X = Xp + K*(z(:,k) - zp);
    P = Pp - K*S*K';
    P = (P + P')/2;
    x_est(:,k) = X(1:2);
end
end

function Xn = f_etat(X, dt)
%F_ETAT  Evolution non lineaire (virage a vitesse constante).
v = X(3);  theta = X(4);  omega = X(5);
Xn = [ X(1) + v*cos(theta)*dt;
    X(2) + v*sin(theta)*dt;
    v;
    theta + omega*dt;
    omega ];
end

function zm = h_mesure(X)
%H_MESURE  On observe uniquement la position (x, y).
zm = X(1:2);
end

function chi = points_sigma(X, P, lambda)
%POINTS_SIGMA  2n+1 points sigma (matrice n x (2n+1)) : X, X+L, X-L avec L*L' = (n+lambda)*P.
n = numel(X);
L = cholesky_securisee((n + lambda)*P);
chi = [X, X + L, X - L];
end

function [ym, Py, Pxy] = transformation_unscented(fun, X, P, Wm, Wc, lambda)
%TRANSFORMATION_UNSCENTED  Propage (X,P) a travers fun : moyenne, covariance, covariance croisee.
chi = points_sigma(X, P, lambda);
nb  = size(chi, 2);
m   = numel(fun(chi(:,1)));
Y   = zeros(m, nb);
for i = 1:nb, Y(:,i) = fun(chi(:,i)); end    % propagation de chaque point
ym  = Y*Wm(:);                               % moyenne ponderee
dY  = Y - ym;   dX = chi - X;
Py  = dY*diag(Wc)*dY';
Pxy = dX*diag(Wc)*dY';
end

function L = cholesky_securisee(A)
%CHOLESKY_SECURISEE  Racine carree avec regularisation (jitter) si A n'est pas definie positive.
A = (A + A')/2;  n = size(A, 1);
[L, p] = chol(A, 'lower');  if p == 0, return; end
jitter = 1e-9*max(1, trace(A)/n);
for essai = 1:10
    [L, p] = chol(A + jitter*eye(n), 'lower');
    if p == 0, return; end
    jitter = 10*jitter;
end
[V, D] = eig(A);
L = V*diag(sqrt(max(diag(D), 1e-12)));
end