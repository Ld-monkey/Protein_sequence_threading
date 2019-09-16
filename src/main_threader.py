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

    Asn_Asn = list()

    # Ouvrir le fichier dope.par qui contient les forces statiques.
    with open("../data/2019-13-10/dope.par", "r") as dope_file:
        for line in dope_file:
            # Selectionner les couples AA en fonction de leurs carbones alphas.
            if re.search("[A-Z]{3}\sCA\s[A-Z]{3}\sCA", line):
                if re.search("^ASN\s[A-Z]{2}\sASN\s[A-Z]{2}", line):
                    Asn_Asn = line[13:].rsplit()

    # Création de la colonne pour le dataframe.
    x = [i for i in np.arange(0.25, (0.25*31), 0.25, float)]

    # Création d'un dataframe de potentiel statistique pour le couple
    # ASN - CA - ASN - CA
    #print("Pour ASN - CA - ASN - CA :")
    Asn_Asn_dataframe = pd.DataFrame([Asn_Asn], columns = x, index = ["E"])
    #print(Asn_Asn_dataframe)

    # Association d'un couple AA en fonction de leur distance
    # pour obtenir leurs de potentiels statistes respectifs.

    # Example de distance récupérer pour un couple AA : example ASN - ASN.
    distance = 0.0

    # Extraire le potentiel statistique par rapport a la distance la plus proche.
    for i in range(0, len(Asn_Asn_dataframe.columns)):
        if Asn_Asn_dataframe.columns[i] >= distance:
            #print(Asn_Asn_dataframe.iloc[0,i])
            # Valeur absolue pour faciliter la programmation dynamique.
            #print(abs(float(Asn_Asn_dataframe.iloc[0,i])))
            break

    pairwise_amino_acide = []

    # Extraire l'ensemble des couples AA pour notre sequence.
    for i in  range(0,len(amino_acide_array['AA'])):
        for y in  range(i,len(amino_acide_array['AA'])):
            pairwise_amino_acide.append([amino_acide_array.iloc[i,0],
                                          amino_acide_array.iloc[y,0]])

    newline()

    # Affiche la table de tout les couples possibles.
    #print(pairwise_amino_acide)

    # Créer un liste d'expression regulière
    expr_regular_pairwise = []

    # Faire un tableau d'expressions régulières.
    for i in range(0, len(pairwise_amino_acide)):
        expr_regular = "^"+pairwise_amino_acide[i][0]+"\s[A-Z]{2}\s"+pairwise_amino_acide[i][1]+"\s[A-Z]{2}"
        expr_regular_pairwise.append([expr_regular])

    #print(expr_regular_pairwise)
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

    #print(potentiel_statistique_dict.keys())

    # Par exemple pour le coucple LEU SER = LEUSER
    #print(potentiel_statistique_dict.get('LEUSER'))

    # Ensuite en fonction qui créer un low matricre (dataframe) associé a une
    # une séquence donnée et qui pioche toutes les informations


    # Extraire le potentiel statistique par rapport a la distance la plus proche.
    for i in range(0, len(potentiel_statistique_dict['ASNASN'].columns)):
        if potentiel_statistique_dict['ASNASN'].columns[i] >= distance:
            #print(potentiel_statistique_dict['ASNASN'].iloc[0,i])
            #print(abs(float(potentiel_statistique_dict['ASNASN'].iloc[0,i])))
            break

    #low_matrix = pd.Dataframe(data = , index = , columns = )
    # Pour constuire la low-matrix : pour chaque couples AA recupérer
    #la distance et la convertir en potentiel statistique.

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
    print("Première low_matrice :")
    print(low_matrix_seq)

    # Atention erreur grossière peut etre que dans le potentil statistique:
    # ASN -> LEU c'est pas la meme chose que LEU - > ASN a verifier !!!.
    # car on a developper un algorithme qui concidère que c'est la meme chose.
