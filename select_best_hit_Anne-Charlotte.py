#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 13:15:46 2020

@author: annecharlottemichel
"""

import os as os
import pandas as pd


liste_fichiers = os.listdir('./blast_outputs_d')
liste_fichiers.remove('.DS_Store')
nb_fichiers = len(liste_fichiers)

os.makedirs('best_hits')

for n in range(0, nb_fichiers):
    table = open('table.csv', 'w')
    table.write("query id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tgaps\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tquery length\tsubject length\n")
    
    fichier = open('./blast_outputs_d/'+liste_fichiers[n], 'r')
    liste_lignes = fichier.readlines()
    for ligne in liste_lignes:
        if ligne[0]!='#':
            table.writelines(ligne)
            
    table.close()
    
    table = pd.read_csv('table.csv', sep ='\t')
    
    Best_hits = table[['query id', 'subject id', '% identity', 'evalue', 
                    'alignment length', 'query length']]
    
    del table
    os.remove('table.csv')
    
    Best_hits = Best_hits[Best_hits['% identity'] > 95]
    Best_hits = Best_hits[Best_hits['evalue'] < 1e-10]
    Best_hits['couverture'] = Best_hits['alignment length']/ Best_hits['query length'] 
    Best_hits = Best_hits[Best_hits['couverture'] > 0.75]
    Best_hits = Best_hits[Best_hits['query id'] != Best_hits['subject id']]
    
    counts = Best_hits['query id'].value_counts()
    
    i = 0
    while counts[i] > 1 and i < len(counts) :
        liste_lignes = list(Best_hits[Best_hits['query id'] == counts.index[i]].index)
        liste_lignes = liste_lignes[1:len(liste_lignes)]
        Best_hits = Best_hits.drop(index=liste_lignes)
        i += 1
    
    Best_hits = Best_hits[['query id', 'subject id']]
    
    Best_hits.to_csv('./best_hits/'+liste_fichiers[n][:-3]+'.txt', sep=' ', index=False, header=False)




