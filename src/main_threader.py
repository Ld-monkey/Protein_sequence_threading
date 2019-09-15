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
    print("Pour ASN - CA - ASN - CA :")
    Asn_Asn_dataframe = pd.DataFrame([Asn_Asn], columns = x, index = ["E"])
    print(Asn_Asn_dataframe)

    # Association d'un couple AA en fonction de leur distance
    # pour obtenir leurs de potentiels statistes respectifs.

    # Example de distance récupérer pour un couple AA : example ASN - ASN.
    distance = 0.0

    # Extraire le potentiel statistique par rapport a la distance la plus proche.
    for i in range(0, len(Asn_Asn_dataframe.columns)):
        if Asn_Asn_dataframe.columns[i] >= distance:
            print(Asn_Asn_dataframe.iloc[0,i])
            # Valeur absolue pour faciliter la programmation dynamique.
            print(abs(float(Asn_Asn_dataframe.iloc[0,i])))
            break

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

    #print(expr_regular_pairwise)
    #print(expr_regular_pairwise[0][0])

    # Lorsque dans un la comparaison de la ligne on trouve une paire
    # AA correspondant au notre liste de tous les couples possibles
    #==> on crée un dataframe associée.
    # Extraire le potentiel statistique pour chaque couple et le convertir
    # en dataframe.

    # Exemple de string contennant les potentiels statistique.
    string_file = "ALA CA ALA CA 10.00 10.00 10.00 10.00 10.00 10.00 10.00 -1.64 1.02 -0.00 -1.23 -0.62 -0.63 0.33 0.58 0.48 0.10 -0.33 -0.07 -0.30 -0.25 -0.16 0.02 0.03 -0.08 0.01 -0.02 -0.08 -0.12 -0.02"+"ALA CA LEU CA 10.00 10.00 10.00 10.00 10.00 10.00 4.41 0.89 -0.94 -1.33 -0.11 -0.57 -0.21 -0.47 0.08 0.20 0.12 -0.11 -0.16 -0.30 -0.06 0.01 -0.07 0.01 0.06 0.01 -0.00 -0.04 -0.03 -0.02"

    print(pairwise_amino_acide[0][0]+
          pairwise_amino_acide[0][1])

    # Création de la colonne pour le dataframe.
    x = [i for i in np.arange(0.25, (0.25*31), 0.25, float)]

    # Création d'un dictionnaire avec les cles AA.
    potentiel_statistique_dict = dict()

    compteur_inutile = 0

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
                        compteur_inutile +=1

    print(compteur_inutile)

    #print(potentiel_statistique_dict.keys())

    # Par exemple pour le coucple LEU SER = LEUSER
    print(potentiel_statistique_dict.get('LEUSER'))

    # Ensuite en fonction qui créer un low matricre (dataframe) associé a une
    # une séquence donnée et qui pioche toutes les informations 
