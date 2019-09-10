def begin_message():
    print("THREADER - Protein Sequence Threading Program")
    print("Built date : September 2019")
    print("usage : python3 main_threader.py\n")

if __name__ == '__main__':
    # Message du début de programme.
    begin_message()

    # Déclaration des variables.
    dict_seq_c_alpha = {}

    # Création d'une liste de dictionnaire.
    list_dict_seq_c_a = list()

    # Ouverture du fichier pdb.
    with open("../data/2019-09-10/2xri.pdb", "r") as contents_pdb:
        numero_list = 0
        for line in contents_pdb:

            # Si c'est un atome avec un carbone alpha
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                numero_list += 1
                dict_seq_c_alpha["num_list"] = numero_list
                dict_seq_c_alpha["Resid_name"] = line[17: 20].strip()
                dict_seq_c_alpha["Resid_seq_num"] = line[22: 26].strip()
                dict_seq_c_alpha["x"] = line[30: 38].strip();
                dict_seq_c_alpha["y"] = line[38: 46].strip();
                dict_seq_c_alpha["z"] = line[46: 54].strip();
                list_dict_seq_c_a.append(dict_seq_c_alpha.copy())
            # Si on veut prendre seulement un nombre restreint d'acides aminés
            if numero_list == 10:
                break

    # Affichage de la liste de dictionnaire.
    for i in range(0, len(list_dict_seq_c_a)):
        print(list_dict_seq_c_a[i])
