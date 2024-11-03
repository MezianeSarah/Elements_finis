# Projet d'Éléments Finis

Bienvenue dans le projet d'analyse par éléments finis, qui fait partie de  mes projets en Master 1 en Modélisation et Analyse Numérique à l'Université de Montpellier.

## Description

Ce projet se concentre sur la résolution de problèmes aux limites en utilisant la méthode des éléments finis, il explore l'application de cette méthode pour discrétiser des équations différentielles sur l'intervalle ouvert \( ]0,1[ \) et évaluer l'impact du maillage sur la précision des solutions obtenues.

## Contenu du Repository

Dans ce repository, vous trouverez les fichiers suivants :

- **non_uniforme.py** : Code Python pour la mise en œuvre de la méthode des éléments finis avec un maillage non uniforme.
- **uniforme.py** : Code Python pour la mise en œuvre de la méthode des éléments finis avec un maillage uniforme.
- **Elements_finis.pdf** : Document de mémoire qui présente le projet en détail, y compris les méthodes utilisées, les résultats obtenus et les analyses réalisées.

## Objectif

Les principaux objectifs de ce projet sont :
- Discrétiser un problème aux limites à l'aide de la méthode des éléments finis.
- Comparer l'impact des maillages uniformes et non-uniformes sur la précision des résultats.
- Vérifier que la fonction \( u(x) = \sin(\pi x^2) \) est la solution exacte du système discrétisé.

## Utilisation

Pour exécuter les scripts Python, assurez-vous d'avoir Python installé sur votre machine. Vous pouvez exécuter chaque fichier en ligne de commande comme suit :

```bash
python uniforme.py
