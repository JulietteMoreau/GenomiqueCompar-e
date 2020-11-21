#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 09:30:17 2020

@author: annecharlottemichel
"""

import os as os
import pandas as pd

liste_bacterie = os.listdir('./prot')
nb_bact = len(liste_bacterie)

# Création du dictionnaire
dict_final = {}
clefs = liste_bacterie
correspondance = {}

for i in range(0, nb_bact):
    dict_final[liste_bacterie[i][:-3]]={}
    correspondance[liste_bacterie[i][:-3]] = ''
    
for i in range(0, nb_bact):
    for j in range(i+1, nb_bact):

        bact_X = liste_bacterie[i][:-3]
        bact_Y = liste_bacterie[j][:-3]
        
        fichier = pd.read_csv('./reciproque/R-' +bact_X +'-' +bact_Y +'.txt', 
                              sep = ' ', header = None)
        
        clefs_bact_X = list(fichier[0])
        clefs_bact_Y = list(fichier[1])
        
        if correspondance[bact_X] == '' :
            n = 0
            nom = ''
            carac = fichier.loc[0,0][n]
            while carac != '_':
                nom = nom + carac
                n += 1
                carac = fichier.loc[0,0][n]
            correspondance[bact_X] = nom
            
        elif correspondance[bact_Y] == '' :
            n = 0
            nom = ''
            carac = fichier.loc[0,1][n]
            while carac != '_':
                nom = nom + carac
                n += 1
                carac = fichier.loc[0,1][n]
            correspondance[bact_Y] = nom
        
        for k in range(0, len(fichier)):
            if clefs_bact_X[k] not in dict_final[bact_X].keys():
                dict_final[bact_X][fichier.loc[k,0]] = [fichier.loc[k,0], fichier.loc[k,1]]
            else :
                dict_final[bact_X][fichier.loc[k,0]].append(fichier.loc[k,1])
            
            if clefs_bact_Y[k] not in dict_final[bact_Y].keys():
                dict_final[bact_Y][fichier.loc[k,1]] = [fichier.loc[k,1], fichier.loc[k,0]]
            else :
                dict_final[bact_Y][fichier.loc[k,1]].append(fichier.loc[k,0])
    
del fichier
del clefs
del clefs_bact_X
del clefs_bact_Y
del bact_Y
del bact_X


# Sélection du core génome
bact = liste_bacterie[0][:-3]
nb_core_genome = 0

correspondance = pd.DataFrame.from_dict(correspondance, orient = 'index')

correspondance = correspondance.reset_index()
        
# Parcours les gènes du génome de la 1ère bactérie
for j in dict_final[bact].keys():
    
    # Vérification si le gène a bien 20 orthologues
    if len(dict_final[bact][j]) == nb_bact : 
        core_genome = True
        k = 1
        
        # Parcours de la liste des orthologues pour vérifier que ces gènes
        # appartiennent au core génome
        while core_genome == True and k < nb_bact :
            gene = dict_final[bact][j][k]
            nom = ''
            
            # Récupération du nom de la bactérie 
            genome = ''
            n = 0
            carac = gene[n]
            while carac != '_':
                genome += carac
                n += 1
                carac = gene[n]
            
            for l in range(0, nb_bact): 
                if genome == correspondance.iloc[l,1] : 
                    nom = correspondance.iloc[l,0]
                    break
            
            # Vérification que les listes contiennent les mêmes éléments
            if (gene in dict_final[nom].keys()
                and sorted(dict_final[nom][gene]) == sorted(dict_final[bact][j])):
                k += 1
                
            else :
                core_genome = False
        print(k)  
        
        if core_genome == True :
            nb_core_genome += 1
                    
            





        

