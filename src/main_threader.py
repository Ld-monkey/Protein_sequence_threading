
import pandas as pd
import numpy as np
import math

# Message d'information apparait à l'ouverture du programme.
def begin_message():
    print("THREADER - Protein Sequence Threading Program")
    print("Built date : September 2019")
    print("usage : python3 main_threader.py\n")

# Fonction qui ajoute une nouvelle ligne en output.
def newline():
    print("")

# Retourne la distance euclidienne entre 2 acides aminés.
def euclidean_distance(dataframe, index, position = 1):
    # Coordonnée euclidienne de premier AA.
    xa = float(dataframe.iloc[index, 2])
    ya = float(dataframe.iloc[index, 3])
    za = float(dataframe.iloc[index, 4])

    # Coordonnée euclidienne du second AA.
    xb = float(dataframe.iloc[(position), 2])
    yb = float(dataframe.iloc[(position), 3])
    zb = float(dataframe.iloc[(position), 4])

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

                # Ajout dans le dictionnaire.
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

    newline()

    # Faire une boucle pour construire un matrice de conctat.
    values_numpy = []
    index_aa = list(amino_acide_array['AA'])

    # Création d'une matrice de contact dans une dataframe.
    for i in range(0, len(amino_acide_array)):
        values_list = list()
        for y in range(0, len(amino_acide_array)):
            values_list.append(euclidean_distance(amino_acide_array, i, y))
        values_numpy.append(values_list)

    # Assemblage de la DataFrame.
    dataframe_aa = pd.DataFrame(data = values_numpy,
                                index = index_aa,
                                columns = index_aa)

    # Afficher la Dataframe (matrice de distance).
    print(dataframe_aa)
