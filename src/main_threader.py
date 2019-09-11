import pandas as pd
import math

# Message d'information apparait à l'ouverture du programme.
def begin_message():
    print("THREADER - Protein Sequence Threading Program")
    print("Built date : September 2019")
    print("usage : python3 main_threader.py\n")

"""
# Construire une classe Acide_aminé qui contiendrait toutes
les caractéristiques d'un AA => son nom et coordonnée x, y, z.
Mais attention on risque d'avoir des AA de meme nom.
Donc est-ce intelligent de construire plusieurs objets ou
une matrice.  bonne question

# Pour un ensemble d'objet (pas mal).
def matrice_distance(residue A, residue B):
  calcul euclidien

# Si matrice
def matrice_distance(dataframe, indice)
  calcul euclidien
"""

# Retourne la distance euclidienne entre 2 acides aminés.
def euclidean_distance(dataframe, index):
    # Coordonnée euclidienne de premier AA.
    xa = float(dataframe.iloc[index, 2])
    ya = float(dataframe.iloc[index, 3])
    za = float(dataframe.iloc[index, 4])

    # Coordonnée euclidienne du second AA.
    xb = float(dataframe.iloc[(index + 1), 2])
    yb = float(dataframe.iloc[(index + 1), 3])
    zb = float(dataframe.iloc[(index + 1), 4])

    # Application de la formule euclidienne.
    distance = math.sqrt(((xb - xa)**2 + (yb - ya)**2 + (zb - za)**2 ))
    return distance


if __name__ == '__main__':


    # Message du début de programme.
    begin_message()

    # Définition des clées à utiliser pour DataFrame.
    dict_seq_c_alpha = {'AA':[],
                        'num_list':[],
                        'x':[],
                        'y':[],
                        'z':[]}

    # L'index de la Dataframe correspondant a une suite AA dans un tableau.
    index_label = []

    # Ouverture du fichier pdb.
    with open("../data/2019-09-10/2xri.pdb", "r") as contents_pdb:
        numero_list = 0
        for line in contents_pdb:

            # Si c'est un atome avec un carbone alpha
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                numero_list += 1
                #dict_seq_c_alpha["num_list"] = numero_list

                # Ajout de l'index dans un tableau.
                index_label.append(numero_list)

                # Ajout dans un le dictionnaire.
                dict_seq_c_alpha['AA'].append(line[17: 20].strip())
                dict_seq_c_alpha['num_list'].append(line[22: 26].strip())
                dict_seq_c_alpha['x'].append(line[30: 38].strip())
                dict_seq_c_alpha['y'].append(line[38: 46].strip())
                dict_seq_c_alpha['z'].append(line[46: 54].strip())
            # Si on veut prendre seulement un nombre restreint d'acides aminés
            if numero_list == 10:
                break

    # Construction de la DataFrame.
    amino_acide_array = pd.DataFrame(data = dict_seq_c_alpha,
                                     index = index_label)

    # Affichage de la DataFrame.
    print(amino_acide_array)

    # Extration des données de la DataFrame.
    # Extraction de la premier ligne
    #print(amino_acide_array.iloc[0, [0, 2, 3, 4]])
    # Extraction de la deuxième ligne
    #print(amino_acide_array.iloc[1, [0, 2, 3, 4]])

    # Application de la fonction euclidienne.
    print("La distance euclidienne entre {} et {} = {} Angstrom.".format(
        amino_acide_array.iloc[0, 0],
        amino_acide_array.iloc[1, 0],
        euclidean_distance(amino_acide_array, index = 0)))

    # Faire une boucle pour construire un matrice de conctat.

