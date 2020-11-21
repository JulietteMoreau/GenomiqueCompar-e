#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 14:10:50 2020

@author: annecharlottemichel
"""

import os as os
import pandas as pd

liste_fichiers = os.listdir('./best_hits')
liste_bacterie = os.listdir('./prot')

#os.makedirs('reciproque')

for i in range (0, len(liste_bacterie)):
    print(i)
    for j in range(i+1, len(liste_bacterie)):
        bact_X = liste_bacterie[i][:-3]
        bact_Y = liste_bacterie[j][:-3]
        table_1 = pd.read_csv('./best_hits/'+bact_X+'-vs-'+bact_Y+'.txt', 
                              sep =' ', header=None)
        table_2 = pd.read_csv('./best_hits/'+bact_Y+'-vs-'+bact_X+'.txt', 
                              sep =' ', header = None)
        
        Reciproque = pd.merge(table_1, table_2, how = 'inner', 
                              left_on = [0, 1], right_on=[1, 0], 
                              suffixes = ("_table_1", "_table_2"))
        Reciproque = Reciproque[['0_table_1', '1_table_1']]
        Reciproque.to_csv('./reciproque/R-'+bact_X+'-' + bact_Y +'.txt', 
                          sep=' ', index=False, header=False)

