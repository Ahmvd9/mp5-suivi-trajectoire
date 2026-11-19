# MP5 - Suivi de trajectoire bruitée

Traitement du signal, EPT DIC1. Mini-projet 5 : suivre une cible à trajectoire
circulaire à partir de mesures de position bruitées, et comparer quatre filtres.

## Problème
- Cible : cercle de rayon R = 100 m, vitesse angulaire ω = 0,1 rad/s
  (vitesse linéaire 10 m/s), durée T = 60 s, pas dt = 0,1 s (N = 600 échantillons).
- Mesures : position (x, y) + bruit gaussien d'écart-type σ = 1, 2, 5, 10 ou 20 m.
- Objectif : estimer la vraie position et comparer les filtres.

## Filtres comparés
| Filtre | Principe | Ce qu'il connaît au départ |
| --- | --- | --- |
| Kalman linéaire | modèle à vitesse constante | 1ère mesure, σ |
| Wiener (FIR) | filtre optimal appris sur une simulation | statistiques apprises sur le cercle |
| EKF | modèle de virage non linéaire, linéarisé | 1ère mesure, σ |
| UKF | même modèle, points sigma (sans Jacobienne) | 1ère mesure, σ |
| EKF info / UKF info | mêmes filtres, initialisés avec R, ω et le cap | rayon, ω, centre du cercle |

Les quatre filtres « de zéro » (Kalman, Wiener, EKF, UKF) sont les filtres principaux.
Les versions « info » servent à mesurer ce que vaut la connaissance a priori de la trajectoire.

## Lancer l'évaluation
1. Mettre tous les fichiers `.m` dans le même dossier Matlab.
2. Exécuter dans la fenêtre de commande : `run_experiments`
3. Compter quelques minutes (l'UKF est le plus lent).

Sorties : `resultats_mp5.mat`, tableaux `tab_*.csv`, figures `fig_*.png`.

## Protocole expérimental
- **Graine 999** : réglage des filtres (Kalman q, EKF/UKF qa et qo, apprentissage de Wiener).
- **Graines 990 à 998** : validation du choix de l'ordre L de Wiener.
- **Graines 1 à 100** : évaluation finale. Elles ne servent jamais au réglage.
- Tous les filtres reçoivent exactement les mêmes données pour un même σ et une même graine.
- Interface commune : `x_est = filtre_xxx(z, params)`, z de taille 2 x N, x_est de taille 2 x N.
  Le bruit de mesure est lu dans `params.sigma_mes` (erreur si absent).

## Métriques
RMSE de position, temps de convergence (erreur lissée sous 2σ), erreur en régime
permanent (derniers 20 % de la durée), temps de calcul, et RMSE du pire cas sur les graines.

## Paramètres retenus (réglés sur la graine 999)
| Paramètre | Valeur |
| --- | --- |
| Kalman q | 1 |
| Wiener, ordre L par σ (1, 2, 5, 10, 20) | 10, 20, 20, 50, 50 |
| EKF de zéro : qa, qo | 1e-2, 1e-6 |
| UKF de zéro : qa, qo | 1e-4, 1e-8 |

## Résultats : RMSE moyen (m) sur 100 graines
| Filtre | σ = 1 | 2 | 5 | 10 | 20 |
| --- | --- | --- | --- | --- | --- |
| Brut (mesures) | 1,42 | 2,83 | 7,08 | 14,15 | 28,31 |
| Kalman | 0,66 | 1,18 | 2,56 | 4,66 | 8,51 |
| Wiener | 0,89 | 1,33 | 3,04 | 4,71 | 8,69 |
| EKF | 0,41 | 0,75 | 1,68 | 3,22 | 6,53 |
| UKF | 0,33 | 0,63 | 1,52 | 3,02 | 5,97 |
| EKF info | 0,32 | 0,60 | 1,38 | 2,60 | 4,91 |
| UKF info | 0,32 | 0,60 | 1,38 | 2,60 | 4,93 |

Autres tableaux dans `Résultats/` : erreur en régime permanent, pire cas, temps de calcul,
temps de convergence.

## Limites
- La trajectoire est un cercle parfait : un bruit de processus très faible y est optimal.
  Sur une cible qui manœuvre, ce ne serait plus le cas.
- Le temps de convergence est nul presque partout : cette métrique ne distingue pas les filtres.
- Wiener est appris sur une seule réalisation : son réglage est sensible à l'ordre L.
- Le temps de calcul de Wiener inclut son apprentissage ; les temps sont mesurés sous MATLAB Online.
- À fort bruit (σ = 20), l'EKF de zéro est plus lent que l'UKF à identifier la vitesse
  angulaire sur certaines graines (voir `Résultats/fig_pire_graine.png`).

## Volet avancé (benchmark)
Conformément au sujet, aucune méthode avancée n'est implémentée par le groupe : on s'appuie
sur les fonctions Matlab existantes et sur la littérature.

| Méthode | Comment elle est évaluée |
| --- | --- |
| LSTM / GRU | testés sous Matlab avec la Deep Learning Toolbox (`lstmLayer`, `gruLayer`, `trainnet`) |
| Filtre particulaire (SIR) | littérature (`trackingPF` n'est pas disponible sur notre licence) |
| Deep Kalman Filter | littérature (Krishnan et al. 2015, KalmanNet 2022) |

Script `benchmark_avance.m` :
- mêmes graines de test (1 à 100) et mêmes métriques que `run_experiments` ;
- réseaux entraînés sur des **cercles aléatoires** (rayon, vitesse, sens, centre et phase tirés
  au hasard, graines 1001 à 1200, validation 990 à 998) : ils ne connaissent pas le cercle de test,
  comme les EKF/UKF « de zéro » ;
- niveaux de bruit testés : σ = 1, 5 et 20 m.

Lancer :
1. Prérequis : Deep Learning Toolbox (R2024a ou plus récent).
2. Exécuter `benchmark_avance` avec `cfg.quick = true` (test de quelques minutes).
3. Passer `cfg.quick = false` pour le calcul complet (environ 15 à 30 minutes).

Sorties : `tab_avance_*.csv`, `tab_complet_*.csv` (fusion avec `resultats_mp5.mat`),
`fig_complet_rmse.png`, `resultats_avance.mat`.

### Résultats du volet avancé : RMSE moyen (m) sur 100 graines
| Méthode | σ = 1 | 5 | 20 |
| --- | --- | --- | --- |
| Brut (mesures) | 1,42 | 7,08 | 28,31 |
| Kalman | 0,66 | 2,56 | 8,51 |
| UKF (de zéro) | 0,33 | 1,52 | 5,97 |
| LSTM | 2,71 | 4,65 | 9,76 |
| GRU | 1,89 | 4,96 | 10,35 |

- Les réseaux font mieux que les mesures brutes à σ = 5 et 20 m, mais **moins bien** à σ = 1 m :
  ils gardent une erreur résiduelle de quelques mètres, même quand le bruit est faible.
- Ils restent 1,6 à 3 fois moins précis que l'UKF, qui exploite le modèle physique du mouvement.
- Coût : environ 13 000 (GRU) à 17 000 (LSTM) paramètres appris, 1 à 2 minutes d'entraînement par σ,
  contre deux paramètres de bruit réglés pour l'EKF et l'UKF.

Figure : `Résultats/fig_complet_rmse.png`.

## Fichiers
| Fichier | Rôle |
| --- | --- |
| get_params.m | Paramètres communs (à ne pas modifier seul) |
| simulate_trajectory.m | Trajectoire vraie + mesures bruitées |
| compute_metrics.m | RMSE, temps de convergence, erreur en régime permanent |
| filtre_kf.m | Kalman linéaire (vitesse constante) |
| filtre_wiener.m | Filtre de Wiener FIR |
| filtre_ekf.m, filtre_ukf.m | EKF / UKF qui partent de zéro |
| filtre_ekf_info.m, filtre_ukf_info.m | EKF / UKF informés (connaissent le cercle) |
| tune_kf.m, tune_wiener.m, tune_ekf_ukf.m | Réglage des paramètres (graine 999) |
| run_kf_demo.m | Premier test du Kalman (graine 999, σ = 5) |
| run_experiments.m | Évaluation Monte Carlo finale (graines 1 à 100) |
| benchmark_avance.m | Volet avancé : test des LSTM / GRU (Deep Learning Toolbox) |
| Résultats/ | Résultats bruts, tableaux et figures de l'évaluation finale |

## Équipe
- Membre 1 : Cheikhna Ahmad Tidiane Diallo (simulation, Kalman, Wiener, évaluation commune)
- Membre 2 : Ndeye Yaram Diop (EKF, UKF)
- Membre 3 : Aminata Ndack Fall (volet avancé: Benchmark)
