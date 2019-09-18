import pandas as pd
import numpy as np
import math
import re
import sys

def begin_message():
    """ Méthode qui affiche différentes informations au démarrage du programme."""

    print("THREADER - Protein Sequence Threading Program")
    print("Built date : September 2019")
    print("usage : python3 main_threader.py\n")

def newline():
    """ Méthode qui ajoute une nouvelle ligne en output."""
    print("")

class ParserPdb_Dope:
    """
    Cette classe permet d'ouvrir le fichier pdb et dope.par.
    """
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file

    def read_pdb_file_to_dataframe(self):
        """ Méthode qui ouvre un fichier pdb et retourne sous forme de dataframe \
        les informations nécessaires.
        """
        dict_seq_c_alpha = {'AA':[],
                            'num_list':[],
                            'x':[],
                            'y':[],
                            'z':[]}

        with open(self.pdb_file, "r") as contents_pdb:
            numero_list = 0
            for line in contents_pdb:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    numero_list += 1
                    dict_seq_c_alpha['AA'].append(line[17: 20].strip())
                    dict_seq_c_alpha['num_list'].append(line[22: 26].strip())
                    dict_seq_c_alpha['x'].append(line[30: 38].strip())
                    dict_seq_c_alpha['y'].append(line[38: 46].strip())
                    dict_seq_c_alpha['z'].append(line[46: 54].strip())
                    if numero_list == 10:
                        break

            dataframe = pd.DataFrame(data = dict_seq_c_alpha,
                                     index = np.arange(1, numero_list+1, 1))
            return dataframe

class AcideAmine:
    """
    Classe représentant l'ensemble des méthoes utilisés pour les acides aminés.
    """
    def __init__(self, dataframe):
        self.dataframe = dataframe

    def euclidean_distance(self, dataframe, index, position = 1):
        """Méthode aui retourne la distance euclidienne entre 2 acides aminés."""
        xa = float(self.dataframe.iloc[index, 2])
        ya = float(self.dataframe.iloc[index, 3])
        za = float(self.dataframe.iloc[index, 4])

        # Coordonnée euclidienne du second AA.
        xb = float(self.dataframe.iloc[(position), 2])
        yb = float(self.dataframe.iloc[(position), 3])
        zb = float(self.dataframe.iloc[(position), 4])

        distance = math.sqrt(( (xb - xa)**2 + (yb - ya)**2 + (zb - za)**2 ))
        return distance

    def create_distance_matrix(self):
        """
        Méthode qui permet de créer un matrice de distance (= matrice de contact).
        """
        values_numpy = []
        index_and_columns_names = list(self.dataframe['AA'])

        for i in range(0, len(self.dataframe)):
            values_list = []
            for y in range(0, len(self.dataframe)):
                values_list.append(self.euclidean_distance(self.dataframe, i, y))
            values_numpy.append(values_list)

        dataframe_aa = pd.DataFrame(data = values_numpy,
                                    index = index_and_columns_names,
                                    columns = index_and_columns_names)
        return dataframe_aa

def pairwise_amino_acide(dataframe):
    """ Méthode qui revoie l'ensemble des couples AA possibles. """
     # Création de la colonne pour le dataframe.
    pairwise_amino_acide = []

    # Extraire l'ensemble des couples AA pour notre sequence.
    for i in  range(0,len(dataframe['AA'])):
        for y in  range(i,len(dataframe['AA'])):
            pairwise_amino_acide.append([dataframe.iloc[i,0],
                                         dataframe.iloc[y,0]])

    return pairwise_amino_acide

def create_regular_expression(pairwise_amino_acide):
    """
    Méthode qui retourne une liste contenant l'ensemble des
    expressions régulières a partir de tout les couples AA possibles.
    """
    expr_regular_pairwise = []

    for i in range(0, len(pairwise_amino_acide)):
        expr_regular = "^"+pairwise_amino_acide[i][0]+"\s[A-Z]{2}\s"+pairwise_amino_acide[i][1]+"\s[A-Z]{2}"
        expr_regular_pairwise.append([expr_regular])

    return expr_regular_pairwise

def create_potentiel_stat_from_dope(dope_file_par,
                                    pairwise_aa,
                                    expr_regular_pairwise):
    """
    Ouvre le fichier avec les potentiels statistiques pour extraire les
    informations nécessaire et les retroune dans un dictionnaire.
    """
    potentiel_statistique_dict = dict()
    x = [i for i in np.arange(0.0, (0.25*30), 0.25, float)]

    with open(dope_file_par, "r") as dope_file:
        for line in dope_file:
            if re.search("[A-Z]{3}\sCA\s[A-Z]{3}\sCA", line):
                for i in range(0, len(expr_regular_pairwise)):
                    if re.search(expr_regular_pairwise[i][0], line):
                        temporaire_list = list()
                        temporaire_list = line[13:].rsplit()
                        temporaire_name = pairwise_aa[i][0]+pairwise_aa[i][1]
                        potentiel_statistique_dict[temporaire_name] = pd.DataFrame([temporaire_list], columns = x, index = ["E"])

    return potentiel_statistique_dict

def create_empty_hight_matrix(lenght_x, lenght_y, columns_names):
    """ Méthode qui créé et retourne une matrice de haut niveau."""
    hight_matrix = np.full((lenght_x, lenght_y), None)
    hight_matrix_seq = pd.DataFrame(hight_matrix,
                                    columns = columns_names,
                                    index = np.arange(1,
                                                      len(columns_names)+1,
                                                      1))
    return hight_matrix_seq

def create_empty_low_matrix(lenght_x, lenght_y, columns_names):
    """ Méthode qui créé et retourne une matrice de haut niveau."""
    low_matrix = np.full((lenght_x, lenght_y), float(20),)
    low_matrix_seq = pd.DataFrame(low_matrix,
                                  columns = columns_names,
                                  index = np.arange(1,
                                                    len(columns_names)+1,
                                                    1))
    return low_matrix_seq

def find_potentiel_statistique(n, m, p, q, distance, potentiel, temp_key):
    if ( n == p and m == q):
        return 0.0
    if distance > 30*0.25:
        return potentiel[temp_key].iloc[0,-1]
    elif distance == 0.0:
        return 10.0
    else :
        for i in range(0, len(potentiel[temp_key].columns)):
            if potentiel[temp_key].columns[i] >= distance:
                rst = float(potentiel[temp_key].iloc[0,i])
                break

    return rst

def built_low_matrix(pos_x_low_matrix,
                     pos_y_low_matrix,
                     aa_fixed_hight_matrix,
                     pos_x_hight_matrix_fixed_position,
                     pos_y_hight_matrix_fixed_position,
                     ca_fixed_hight_matrix,
                     low_matrix,
                     potentiel_statistique,
                     matrice_distance):

    # Ajout d'une cles
    temps_key = aa_fixed_hight_matrix+low_matrix.columns[pos_y_low_matrix]

    # Calculer une distance.
    distance = matrice_distance.iloc[pos_y_hight_matrix,
                                     pos_x_low_matrix]

    energie = find_potentiel_statistique(aa_fixed_hight_matrix,
                                         ca_fixed_hight_matrix,
                                         low_matrix.columns[pos_y_low_matrix],
                                         low_matrix.index[pos_x_low_matrix],
                                         distance,
                                         potentiel_statistique,
                                         temps_key)

    low_matrix.iloc[pos_x_low_matrix, pos_y_low_matrix] = energie
    return low_matrix

def dynamique_programming(low_matrix, x=1, y=1):
    i = x
    while (low_matrix.iloc[i+1, y+1] != 20 or
           low_matrix.iloc[i, y+1] != 20 or
           low_matrix.iloc[i+1, y] != 20):
        if float(low_matrix.iloc[i, y+1]) < float(low_matrix.iloc[i+1, y+1]) and float(low_matrix.iloc[i, y+1]) < float(low_matrix.iloc[i+1, y]):
            final_score = low_matrix.iloc[i, y+1]
            y += 1
        elif float(low_matrix.iloc[i+1, y+1]) < float(low_matrix.iloc[i, y+1]) and float(low_matrix.iloc[i+1, y+1]) < float(low_matrix.iloc[i+1, y]):
            final_score = low_matrix.iloc[i+1, y+1]
            i += 1
            y += 1
        else:
            final_score = low_matrix.iloc[i+1, y]
            i +=1
    return y, i

if __name__ == '__main__':

    # Message du début de programme.
    begin_message()

    try:
        link_pdb_file = str(sys.argv[1])
        link_dope_file = str(sys.argv[2])
    except IndexError:
        sys.exit("Erreur d'argument : example d'utilisation\n"+
                 "python main_threader.py ../data/pdb/2xri.pdb"+
                 "../data/dope/dope.par")

    pdb = ParserPdb_Dope(link_pdb_file)
    # Extraction des information du pdb.
    pdb_file_dataframe = pdb.read_pdb_file_to_dataframe()

    print("Informations du fichier .pdb :")
    print(pdb_file_dataframe)

    newline()

    all_acide_amine = AcideAmine(pdb_file_dataframe)

    # Création de la matrice de distance ou matrice de contact.
    distance_matrix = all_acide_amine.create_distance_matrix()
    print(distance_matrix)

    # Création de tout les couples AA possibles.
    all_pairwises_aa_list = pairwise_amino_acide(pdb_file_dataframe)

    # Création de l'ensemblde des expressions régulières.
    expr_regular = create_regular_expression(all_pairwises_aa_list)

    # Création d'un dictionnaire stockant l'ensemble des potentiels statistiques.
    pot_stat_dict = create_potentiel_stat_from_dope(link_dope_file,
                                                    all_pairwises_aa_list,

                                                    expr_regular)

    # Création d'une matrice de haut niveau
    print("La matrice de haut niveau (hight matrix)")
    hight_matrix = create_empty_hight_matrix(len(distance_matrix.columns),
                                             len(distance_matrix.columns),
                                             distance_matrix.columns)
    print(hight_matrix)

    # Suite.
    newline()

    aa_hight_matrix = hight_matrix.columns
    ca_hight_matrix = hight_matrix.index

    pos_x_hight_matrix = 0
    pos_y_hight_matrix = 0

    AA = aa_hight_matrix[pos_y_hight_matrix]
    CA = ca_hight_matrix[pos_x_hight_matrix]

    pos_x_low_matrix = 0
    pos_y_low_matrix = 0

    print("Création de la low_matrice.")
    # Création de la low matrice
    low_matrix = create_empty_low_matrix(len(distance_matrix.columns),
                                         len(distance_matrix.columns),
                                         distance_matrix.columns)

    # variable car low matrice
    pos_x_low_matrix = 0
    pos_y_low_matrix = 0

    pos_x_initiale = pos_x_low_matrix
    pos_y_initiale = pos_y_low_matrix

    # Ne bouge pas tant que la low matrice na pas finit

    pos_x_hight_matrix = 0
    pos_y_hight_matrix = 0
    low_matrix = built_low_matrix(pos_x_low_matrix,
                                  pos_y_low_matrix,
                                  AA,
                                  pos_x_hight_matrix,
                                  pos_y_hight_matrix,
                                  CA,
                                  low_matrix,
                                  pot_stat_dict,
                                  distance_matrix)


    # Modification des variable des postion low_matrice
    # On va vers la droite :
    for i in range(0, 5):
        pos_y_low_matrix +=1
        low_matrix = built_low_matrix(pos_x_low_matrix, pos_y_low_matrix, AA, pos_x_hight_matrix,
                                      pos_y_hight_matrix, CA, low_matrix, pot_stat_dict, distance_matrix)
        pos_x_low_matrix +=1
        low_matrix = built_low_matrix(pos_x_low_matrix, pos_y_low_matrix, AA, pos_x_hight_matrix,
                                      pos_y_hight_matrix, CA, low_matrix, pot_stat_dict, distance_matrix)
        temps = pos_y_low_matrix
        pos_y_low_matrix -= 1
        low_matrix = built_low_matrix(pos_x_low_matrix, pos_y_low_matrix, AA, pos_x_hight_matrix,
                                      pos_y_hight_matrix, CA, low_matrix, pot_stat_dict, distance_matrix)

        pos_x_low_matrix, pos_y_low_matrix = dynamique_programming(low_matrix,
                                                               pos_x_initiale,
                                                               pos_y_initiale)

    print(low_matrix)
