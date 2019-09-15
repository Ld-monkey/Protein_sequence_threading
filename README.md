# Protein sequence threading by double dynamic programming

## Introduction :
Développement d'un programme en python 3 permettant de prédire une structure tertiaire a partir d'une séquence d'acide aminée .

## Architecture du projet :
```bash
.
├── bin
├── data
│   └── 2019-09-10
│       ├── 2xri.fasta
│       ├── 2xri.fasta.txt
│       └── 2xri.pdb
├── doc
├── README.md
├── results
└── src
    └── main_threader.py

```

## But du projet :

Dans cette étude nous nous intéressons à la séquence protéique de la 2XRI fournit par la banque de donnée protein data bank (pdb).

## Usage :

Il est préférable d'avoir activé l'environnement conda pour permettre la reproductibilité des données.

```bash
python3 main_threader.py
```

## Resultats :

## Environnement conda :
Pour créer un environnement conda a partir du fichier yaml (conda_threader.yml) :
```bash
conda env create -f conda_threader.yml
```

Pour activer l'environnement conda :
```bash
conda active conda_threader
```

Pour désactiver l'environnement conda :
```bash
conda deactivate
```

Pour supprimer l'environnement :
```
conda env remove -n conda_threader
```

## Objectifs :

- [x] Ajouter un environnement conda avec un fichier yaml .
