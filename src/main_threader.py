import pandas as pd
import numpy as np
import math
import re

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
    distance = math.sqrt(( (xb - xa)**2 + (yb - ya)**2 + (zb - za)**2 ))
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
    print("Sortie du fichier .pdb :")
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
    print("La matrice de distance :")
    print(dataframe_aa)
    # Création de la colonne pour le dataframe.
    x = [i for i in np.arange(0.25, (0.25*31), 0.25, float)]

    print(amino_acide_array)

    pairwise_amino_acide = []

    # Extraire l'ensemble des couples AA pour notre sequence.
    for i in  range(0,len(amino_acide_array['AA'])):
        for y in  range(i,len(amino_acide_array['AA'])):
            pairwise_amino_acide.append([amino_acide_array.iloc[i,0],
                                          amino_acide_array.iloc[y,0]])

    newline()

    # Affiche la table de tout les couples possibles.
    print(pairwise_amino_acide)

    # Créer un liste d'expression regulière
    expr_regular_pairwise = []

    # Faire un tableau d'expressions régulières.
    for i in range(0, len(pairwise_amino_acide)):
        expr_regular = "^"+pairwise_amino_acide[i][0]+"\s[A-Z]{2}\s"+pairwise_amino_acide[i][1]+"\s[A-Z]{2}"
        expr_regular_pairwise.append([expr_regular])

    print(expr_regular_pairwise)
    #print(expr_regular_pairwise[0][0])

    # Lorsque dans un la comparaison de la ligne on trouve une paire
    # AA correspondant au notre liste de tous les couples possibles
    #==> on crée un dataframe associée.
    # Extraire le potentiel statistique pour chaque couple et le convertir
    # en dataframe.

    # Création de la colonne pour le dataframe.
    x = [i for i in np.arange(0.25, (0.25*31), 0.25, float)]

    # Création d'un dictionnaire avec les cles AA.
    potentiel_statistique_dict = dict()

    # On créer un dictionnaire de key car il ne peut exister qu'un seul
    # couple possible associé au nom de ce couple exemple key : ASNASN.
    with open("../data/2019-13-10/dope.par", "r") as dope_file:
        for line in dope_file:
            if re.search("[A-Z]{3}\sCA\s[A-Z]{3}\sCA", line):
                for i in range(0, len(expr_regular_pairwise)):
                    #print(expr_regular_pairwise[i][0])
                    if re.search(expr_regular_pairwise[i][0], line):
                        temporaire_list = list()
                        temporaire_list = line[13:].rsplit()
                        temporaire_name = pairwise_amino_acide[i][0]+pairwise_amino_acide[i][1]
                        #print(temporaire_name)
                        potentiel_statistique_dict[temporaire_name] = pd.DataFrame([temporaire_list], columns = x, index = ["E"])

    print(potentiel_statistique_dict.keys())

    # Récuperer le couple a partir de la matrice de contact
    #print(dataframe_aa.columns[0])
    #print(dataframe_aa.index[1])
    #print(dataframe_aa.iloc[0,0])
    #print(potentiel_statistique_dict.keys())

    #print(pairwise_amino_acide)

    potentiel_statistique_array = list()
    compteur_inutile = 0
    compteur_inutile_2 = 0

    for i in range(0, len(dataframe_aa.columns)):
        for y in range(0, len(dataframe_aa.index)):
            #print(dataframe_aa.index[i]+dataframe_aa.columns[y])
            temporaire_key = dataframe_aa.index[i]+dataframe_aa.columns[y]

            if temporaire_key not in potentiel_statistique_dict.keys():
                temporaire_key = dataframe_aa.index[y]+dataframe_aa.columns[i]
                #print("inverse")
                #print(temporaire_key)
            compteur_inutile += 1

            #print(temporaire_key)
            distance = dataframe_aa.iloc[i,y]
            #print(distance)

            # Si ce sont les memes acides aminés = 0.0
            if dataframe_aa.index[i] == dataframe_aa.columns[y]:
                potentiel_statistique_array.append(float(0.0))
            else :
                for z in range(0, len(potentiel_statistique_dict[temporaire_key].columns)):
                    if potentiel_statistique_dict[temporaire_key].columns[z] >= distance:
                        #print(abs(float(potentiel_statistique_dict[temporaire_key].iloc[0,z])))
                        potentiel_statistique_array.append(abs(float(potentiel_statistique_dict[temporaire_key].iloc[0,z])))
                        compteur_inutile_2 += 1
                        break
                    # cutoff
                    if distance > potentiel_statistique_dict[temporaire_key].columns[-1]:
                        potentiel_statistique_array.append(abs(float(potentiel_statistique_dict[temporaire_key].iloc[0,-1])))
                        compteur_inutile_2 +=1
                        break


    #print(compteur_inutile)
    #print(compteur_inutile_2)

    # On affiche l'ensemble du tableau.
    #print(potentiel_statistique_array)

    #print(index_label)

    potentiel_statistique_array = np.reshape(potentiel_statistique_array, (10,10))
    #print(potentiel_statistique_array)

    # Création de la low matrice.
    low_matrix_seq = pd.DataFrame(potentiel_statistique_array,
                                  columns = index_aa,
                                  index = index_label)

    # Affichage de la première low matrice.
    print("Première low_matrix :")
    print(low_matrix_seq)

    # Faire la première programmation dynamique et associer la dernière valeur
    # obtenue dans la hight matrix.

    # Création de l'index.
    i = 0
    y = 0
    final_score = low_matrix_seq.iloc[0,0]

    newline()

    # programmation dynamique
    while (y < (len(low_matrix_seq.columns)-1) and i < (len(low_matrix_seq.index)-1)):
        if low_matrix_seq.iloc[i, y+1] < low_matrix_seq.iloc[i+1, y+1] and low_matrix_seq.iloc[i, y+1] < low_matrix_seq.iloc[i+1, y]:
            #print("->")
            #print(i,y)
            final_score = low_matrix_seq.iloc[i, y+1]
            y += 1
            #print(i, y)
        elif low_matrix_seq.iloc[i+1, y+1] < low_matrix_seq.iloc[i, y+1] and low_matrix_seq.iloc[i+1, y+1] < low_matrix_seq.iloc[i+1, y]:
            #print("\>")
            final_score = low_matrix_seq.iloc[i+1, y+1]
            i += 1
            y += 1
            #print(i)
            #print(y)
        else:
            final_score = low_matrix_seq.iloc[i+1, y]
            #print("|")
            i +=1

    print("Le score final est :",final_score)

    # Création d'une matrice de haut niveau
    hight_matrix = np.full((10, 10), None)
    hight_matrix_seq = pd.DataFrame(hight_matrix,
                                    columns = index_aa,
                                    index = index_label)
    newline()
    print("La matrice de haut niveau (hight matrix)")

    # Ajouter la valeur dans la hight matrix.
    hight_matrix_seq.iloc[0,0] = final_score

    # Afficher la hight_matrix.
    print(hight_matrix_seq)

    # Une fois la premiere low matrice faire il faut
    # fois les autres propositions de low matrice a
    # réaliser a partir de la position initiale de la
    # high level.

    pos_x_hight_matrix = 0
    pos_y_hight_matrix = 0

    # A droite du point de départ
    print(hight_matrix_seq.columns[1], hight_matrix_seq.index[0])
    right_new_sequence = list(hight_matrix_seq.columns)

    # Applique la modification.
    right_new_sequence[hight_matrix_seq.index[0]-1] = hight_matrix_seq.columns[1]
    print(right_new_sequence)

    # En diagonale du point de départ
    print(hight_matrix_seq.columns[1], hight_matrix_seq.index[1])
    diago_new_sequence = list(hight_matrix_seq.columns)
    diago_new_sequence[hight_matrix_seq.index[1]-1] = hight_matrix_seq.columns[1]
    print(diago_new_sequence)

    # En bas du point de départ
    print(hight_matrix_seq.columns[0], hight_matrix_seq.index[1])
    botton_new_sequence = list(hight_matrix_seq.columns)
    botton_new_sequence[hight_matrix_seq.index[1]-1] = hight_matrix_seq.columns[0]
    print(botton_new_sequence)

    # Création de la low_matrice pour la sequence de droite.
    # (Faire un fonction pour la low_matrice)

    potentiel_statistique_array = list()

    for i in range(0, len(dataframe_aa.columns)):
        for y in range(0, len(dataframe_aa.index)):
            #print(dataframe_aa.index[i]+dataframe_aa.columns[y])
            temporaire_key = dataframe_aa.index[i]+right_new_sequence[y]


            if temporaire_key not in potentiel_statistique_dict.keys():
                temporaire_key = dataframe_aa.index[y]+right_new_sequence[i]
                #print("inverse")
                #print(temporaire_key)
            print(temporaire_key)
            distance = dataframe_aa.iloc[i,y]
            #print(distance)

            # Si ce sont les memes acides aminés = 0.0
            if dataframe_aa.index[i] == dataframe_aa.columns[y]:
                potentiel_statistique_array.append(float(0.0))
            else :
                for z in range(0, len(potentiel_statistique_dict[temporaire_key].columns)):
                    if potentiel_statistique_dict[temporaire_key].columns[z] >= distance:
                        #print(abs(float(potentiel_statistique_dict[temporaire_key].iloc[0,z])))
                        potentiel_statistique_array.append(abs(float(potentiel_statistique_dict[temporaire_key].iloc[0,z])))
                        break
                    # cutoff
                    if distance > potentiel_statistique_dict[temporaire_key].columns[-1]:
                        potentiel_statistique_array.append(abs(float(potentiel_statistique_dict[temporaire_key].iloc[0,-1])))
                        break


    potentiel_statistique_array = np.reshape(potentiel_statistique_array, (10,10))
    #print(potentiel_statistique_array)

    # Création de la low matrice.
    low_matrix_seq = pd.DataFrame(potentiel_statistique_array,
                                  columns = right_new_sequence,
                                  index = index_label)

    # Affichage de la première low matrice.
    print("2nd low_matrix :")
    print(low_matrix_seq)


    # Tout ce que j'ai fait est faux donc on repart de zero.
