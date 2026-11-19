%% benchmark_avance.m - MP5 : volet avance, test des LSTM / GRU sous Matlab
% Utilise uniquement les fonctions existantes de la Deep Learning Toolbox
% (lstmLayer, gruLayer, trainnet) : aucune architecture n'est implementee.
% Le filtre particulaire et le Deep Kalman Filter sont traites par la
% litterature (trackingPF indisponible sur notre licence).
%
% Meme scenario, memes graines de test (1 a 100) et memes metriques que
% run_experiments.m. Les reseaux sont entraines sur des CERCLES ALEATOIRES
% (rayon, vitesse, sens, centre, phase tires au hasard) : ils ne connaissent
% pas le cercle de test, comme nos EKF/UKF "de zero".
% Lancer d'abord avec cfg.quick = true (quelques minutes).

clear; clc; close all;
params = get_params();

%% 1. Parametres
cfg.quick      = false;         % true = test rapide, false = calcul complet
cfg.sigmas     = [1 5 20];      % niveaux de bruit testes (m)
cfg.nTrain     = 200;           % trajectoires d'entrainement par sigma
cfg.seedTrain0 = 1000;          % graines d'entrainement : 1001 a 1200
cfg.seedsVal   = 990:998;       % graines de validation (comme Wiener)
cfg.H          = 64;            % unites cachees du LSTM / GRU
cfg.epochs     = 150;           % on garde le reseau du meilleur score de validation
nMC = params.nMC;               % graines de test 1 a 100

if cfg.quick                    % test : sigma = 5 m, 5 graines, entrainement complet
    cfg.sigmas = 5; cfg.nTrain = 200; cfg.epochs = 60; nMC = 5;
end

archs = {'LSTM','GRU'}; nA = numel(archs);
sigmas = cfg.sigmas; nS = numel(sigmas);

rmse  = nan(nA,nS,nMC); eperm = nan(nA,nS,nMC);
tconv = nan(nA,nS,nMC); tcalc = nan(nA,nS,nMC);
NPAR  = nan(nA,nS);     TTRAIN = nan(nA,nS);

%% 2. Entrainement puis evaluation, pour chaque sigma
for j = 1:nS
    s = sigmas(j);
    fprintf('--- sigma = %g m ---\n', s);

    % Donnees d'entrainement et de validation (cercles aleatoires)
    [XT, YT] = genCercles(params, s, cfg.seedTrain0 + (1:cfg.nTrain));
    [XV, YV] = genCercles(params, s, cfg.seedsVal);

    for a = 1:nA
        % Entrainement
        rng(cfg.seedTrain0);
        tic; [net, NPAR(a,j)] = entrainerReseau(archs{a}, XT, YT, XV, YV, cfg);
        TTRAIN(a,j) = toc;
        fprintf('%s entraine en %.0f s (%d parametres)\n', archs{a}, TTRAIN(a,j), NPAR(a,j));

        % Echauffement (le 1er appel est toujours plus lent)
        [~, z0] = simulate_trajectory(s, 1);
        filtre_dl(z0, net, params.R);

        % Evaluation sur les memes graines que le volet classique
        for n = 1:nMC
            [x_true, z] = simulate_trajectory(s, n);
            tic; x_est = filtre_dl(z, net, params.R); tcalc(a,j,n) = toc;
            met = compute_metrics(x_true, x_est, s);
            rmse(a,j,n) = met.rmse;  tconv(a,j,n) = met.t_conv;
            eperm(a,j,n) = met.err_perm;
        end
        fprintf('%s : RMSE moyen = %.2f m\n', archs{a}, mean(rmse(a,j,:)));
    end
end

%% 3. Tableaux des reseaux seuls
ecrireTab(mean(rmse,3),            'tab_avance_rmse.csv',      archs, sigmas);
ecrireTab(max(rmse,[],3),          'tab_avance_pirecas.csv',   archs, sigmas);
ecrireTab(mean(eperm,3),           'tab_avance_eperm.csv',     archs, sigmas);
ecrireTab(mean(tconv,3,'omitnan'), 'tab_avance_tconv.csv',     archs, sigmas);
ecrireTab(1000*mean(tcalc,3),      'tab_avance_temps_ms.csv',  archs, sigmas);
ecrireTab(NPAR,                    'tab_avance_parametres.csv',archs, sigmas);
ecrireTab(TTRAIN,                  'tab_avance_entrainement_s.csv', archs, sigmas);

%% 4. Comparaison avec le volet classique (resultats_mp5.mat)
fic = 'resultats_mp5.mat';
if ~exist(fic,'file'), fic = fullfile('Résultats','resultats_mp5.mat'); end
if ~cfg.quick && exist(fic,'file')
    G = load(fic);
    [~, iS] = ismember(sigmas, G.params.sigmas);       % colonnes sigma communes
    noms  = [G.noms(:).', archs];
    R_all = cat(1, G.rmse(:,iS,:), rmse);
    E_all = cat(1, G.eperm(:,iS,:), eperm);
    C_all = cat(1, G.tcalc(:,iS,:), tcalc);
    ecrireTab(mean(R_all,3),       'tab_complet_rmse.csv',     noms, sigmas);
    ecrireTab(max(R_all,[],3),     'tab_complet_pirecas.csv',  noms, sigmas);
    ecrireTab(mean(E_all,3),       'tab_complet_eperm.csv',    noms, sigmas);
    ecrireTab(1000*mean(C_all,3),  'tab_complet_temps_ms.csv', noms, sigmas);
    disp(readtable('tab_complet_rmse.csv'));

    % Figure : RMSE en fonction de sigma, filtres du groupe + reseaux
    h = figure; hold on; grid on;
    for f = 1:numel(noms)
        if any(strcmp(noms{f}, archs)), st = '-^';
        elseif contains(noms{f}, 'info'), st = '--s';
        else, st = '-o'; end
        errorbar(sigmas, mean(R_all(f,:,:),3), std(R_all(f,:,:),0,3), st, 'LineWidth', 1.2);
    end
    set(gca,'XScale','log','YScale','log'); xticks(sigmas);
    xlabel('Bruit de mesure \sigma (m)'); ylabel('RMSE (m)');
    legend(noms, 'Location','northwest');
    title(sprintf('RMSE : filtres classiques et reseaux recurrents (%d realisations)', nMC));
    exportgraphics(h, 'fig_complet_rmse.png', 'Resolution', 200);
else
    disp(array2table(mean(rmse,3), 'RowNames', archs, ...
        'VariableNames', compose('sigma_%g', sigmas)));
end

save('resultats_avance.mat', 'cfg','archs','sigmas','rmse','eperm','tconv','tcalc','NPAR','TTRAIN');

%% ======================= Fonctions locales =======================
function x_est = filtre_dl(z, net, echelle)
    % Interface proche des filtres du groupe : z (2 x N) -> x_est (2 x N).
    % Reseau unidirectionnel : chaque estimation n'utilise que le passe (causal).
    Y = minibatchpredict(net, {(z/echelle).'}, UniformOutput=false);
    x_est = (Y{1}*echelle).';
end

function [X, Y] = genCercles(params, sig, graines)
    % Trajectoires circulaires aleatoires + bruit gaussien d'ecart-type sig.
    % Donnees en (temps x canaux), normalisees par le rayon nominal.
    n = numel(graines); X = cell(n,1); Y = cell(n,1);
    N = round(params.T/params.dt); t = (0:N-1)*params.dt;
    for i = 1:n
        rng(graines(i));
        R0  = 50 + 100*rand;                          % rayon 50 a 150 m
        w   = (0.05 + 0.1*rand)*(2*(rand>0.5)-1);     % 0,05 a 0,15 rad/s, deux sens
        c   = -100 + 200*rand(2,1);                   % centre dans [-100,100]^2
        phi = 2*pi*rand;                              % phase initiale
        P = c + R0*[cos(w*t+phi); sin(w*t+phi)];
        X{i} = ((P + sig*randn(size(P)))/params.R).';
        Y{i} = (P/params.R).';
    end
end

function [net, nPar] = entrainerReseau(arch, XT, YT, XV, YV, cfg)
    % Architecture minimale : entree (x,y) -> couche recurrente -> sortie (x,y)
    if strcmp(arch,'LSTM'), rec = lstmLayer(cfg.H, OutputMode="sequence");
    else,                   rec = gruLayer(cfg.H,  OutputMode="sequence"); end
    couches = [sequenceInputLayer(2); rec; fullyConnectedLayer(2)];
    opts = trainingOptions("adam", MaxEpochs=cfg.epochs, MiniBatchSize=32, ...
        InitialLearnRate=5e-3, GradientThreshold=1, Shuffle="every-epoch", ...
        ValidationData={XV, YV}, OutputNetwork="best-validation", ...
        Verbose=false, Plots="none");
    net  = trainnet(XT, YT, couches, "mse", opts);
    nPar = sum(cellfun(@numel, net.Learnables.Value));
end

function ecrireTab(M, fichier, noms, sigmas)
    T = array2table(M, 'VariableNames', compose('sigma_%g', sigmas));
    T = addvars(T, noms(:), 'Before', 1, 'NewVariableNames', 'Methode');
    writetable(T, fichier);
end
