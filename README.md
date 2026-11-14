# MP5 - Suivi de trajectoire bruitee (EPT DIC1, Traitement du signal)

Suivi d'une cible en mouvement circulaire (R = 100 m, omega = 0,1 rad/s, dt = 0,1 s,
T = 60 s) a partir de mesures de position bruitees (sigma = 1, 2, 5, 10, 20 m).
Filtres compares : Kalman lineaire, Wiener (FIR), EKF, UKF.

## Lancer l'evaluation
Placer tous les fichiers .m dans le meme dossier Matlab, puis executer :

    run_experiments

Resultats : resultats_mp5.mat, tableaux (.csv) et figures (.png).

## Protocole experimental
- Graine 999 : reglage des filtres (Kalman q, EKF/UKF qa et qo, apprentissage Wiener).
- Graines 990 a 998 : validation du choix de l'ordre L de Wiener.
- Graines 1 a 100 : evaluation finale (jamais utilisees pour regler).
- Tous les filtres recoivent exactement les memes donnees ; sigma est passe par params.sigma_mes.

## Fichiers
| Fichier | Role |
| --- | --- |
| get_params.m | Parametres communs (a ne pas modifier seul) |
| simulate_trajectory.m | Trajectoire vraie + mesures bruitees |
| compute_metrics.m | RMSE, temps de convergence, erreur en regime permanent |
| filtre_kf.m | Kalman lineaire (vitesse constante) |
| filtre_wiener.m | Filtre de Wiener FIR |
| run_kf_demo.m | Premier test du Kalman (graine 999, sigma = 5) |
| filtre_ekf.m, filtre_ukf.m | EKF / UKF qui partent de zero (ne connaissent pas le cercle) |
| filtre_ekf_info.m, filtre_ukf_info.m | EKF / UKF informes (connaissent R, omega, cap) |
| tune_kf.m, tune_wiener.m, tune_ekf_ukf.m | Reglage des parametres (graine 999) |
| run_experiments.m | Evaluation Monte Carlo finale (graines 1 a 100) |

## Convention d'interface
Chaque filtre : `x_est = filtre_xxx(z, params)`, avec z de taille 2 x N et x_est 2 x N.
Le bruit de mesure est lu dans `params.sigma_mes` (erreur si absent).
## Resultats
Dossier `Résultats/` : resultats_mp5.mat (resultats bruts), tab_*.csv (tableaux), fig_*.png (figures).

## Equipe
- Membre 1 : Cheikhna Ahmad Tidiane Diallo(simulation, Kalman, Wiener, evaluation commune)
- Membre 2 : Ndeye Yaram Diop (EKF, UKF)
- Membre 3 : Aminata Ndack Fall [a completer]
