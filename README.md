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

## Volet avancé
Conformément au sujet (benchmark sans implémentation), le volet avancé — filtre particulaire,
LSTM/GRU et Deep Kalman Filter (KalmanNet) — est traité à partir de la littérature
scientifique, dans le rapport. Aucun code n'est fourni pour ce volet.

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
| Résultats/ | Résultats bruts, tableaux et figures de l'évaluation finale |

## Équipe
- Membre 1 : Cheikhna Ahmad Tidiane Diallo (simulation, Kalman, Wiener, évaluation commune)
- Membre 2 : Ndeye Yarame Diop (EKF, UKF)
- Membre 3 : Aminata Ndack Fall (volet avancé : benchmark bibliographique)
