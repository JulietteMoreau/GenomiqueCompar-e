# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 17:46:46 2020

@author: julie
"""

import os as os
import pandas as pd


######################BEST HITS


liste_fichiers = os.listdir('./Blast_output')
nb_fichiers = len(liste_fichiers)

#os.makedirs('best_hits')

for n in range(0,nb_fichiers):
    table=open('table.csv','w')
    table.write("query id\tsubject id\t% identity\talignment length\tmismatches\tgap opens\tgaps\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\tquery length\tsubject length\n")
    fichier = open('./Blast_output/'+liste_fichiers[n], 'r')
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
    
    counts=Best_hits['query id'].value_counts()
    
    i=0
    while counts[i]>1 and i<len(counts):
        liste_lignes=list(Best_hits[Best_hits['query id']==counts.index[i]].index)
        liste_lignes=liste_lignes[1:len(liste_lignes)]
        Best_hits=Best_hits.drop(index=liste_lignes)
        i+=1
        
    Best_hits=Best_hits[['query id', 'subject id']]
    
    Best_hits.to_csv('./best_hits/'+liste_fichiers[n][:-3]+'.txt', sep=' ', index=False, header='False')
    
    
################RECIPROQUE

liste_fichiers = os.listdir('./best_hits')
liste_bacteries=os.listdir('./prot')

#os.makedirs('reciproque')

for i in range(0,len(liste_bacteries)):
    print(i)
    for j in range(i+1,len(liste_bacteries)):
        bact_X=liste_bacteries[i][:-3]
        bact_Y=liste_bacteries[j][:-3]
        table1=pd.read_csv('./best_hits/'+bact_X+'-vs-'+bact_Y+'.txt', sep=' ', header=None)
        table2=pd.read_csv('./best_hits/'+bact_Y+'-vs-'+bact_X+'.txt', sep=' ', header=None)
        
        Reciproque=pd.merge(table1, table2, how="inner", left_on=[0,1], right_on=[1,0], suffixes=('_table_1', '_table_2'))
        Reciproque=Reciproque[['0_table_1','1_table_1']]
        Reciproque.to_csv('./reciproque/R-'+bact_X+'-'+bact_Y+'.txt', sep=' ', index=False, header=False)
        
        
        
#################CORE GENOME

liste_bacteries=os.listdir('./prot')
liste_reciproque=os.listdir('./reciproque')

dict_final={}
clefs=liste_bacteries

for i in range(0,len(liste_bacteries)):
    dict_final[liste_bacteries[i][:-3]]={}
    
for i in range(0,len(liste_bacteries)):
    for j in range(i+1,len(liste_bacteries)):
        bact_X=liste_bacteries[i][:-3]
        bact_Y=liste_bacteries[j][:-3]
        fichier=pd.read_csv('./reciproque/R-'+bact_X+'-'+bact_Y+'.txt', sep=' ', header=None)
        clefs_bact_X=list(fichier[0])
        clefs_bact_Y=list(fichier[1])
        for k in range(0, len(fichier)):
            if clefs_bact_X[k] not in dict_final[bact_X].keys():
                dict_final[bact_X][fichier.loc[k,0]]=[fichier.loc[k,0],fichier.loc[k,1]]
            else:
                dict_final[bact_X][fichier.loc[k,0]].append(fichier.loc[k,1])

            if clefs_bact_Y[k] not in dict_final[bact_Y].keys():
                dict_final[bact_Y][fichier.loc[k,1]]=[fichier.loc[k,1],fichier.loc[k,0]]
            else:
                dict_final[bact_Y][fichier.loc[k,1]].append(fichier.loc[k,0])
                
del fichier

core_genome=[]     
           
for cle in range(0,len(dict_final[liste_bacteries[0][:-3]])):
    traceur=0
    cles_gene=list(dict_final[liste_bacteries[0][:-3]].keys())
    
    if len(dict_final[liste_bacteries[0][:-3]][cles_gene[cle]])==21:
        
        for genome in range(1,len(dict_final)):
            for gene in dict_final[liste_bacteries[0][:-3]][cles_gene[cle]]:
                
                cles_gene_genome=list(dict_final[liste_bacteries[genome][:-3]].keys())
                     
                if gene in dict_final[liste_bacteries[genome][:-3]].keys():
                    
                    if gene not in dict_final[liste_bacteries[genome][:-3]][gene]:
                        traceur+=1
                    
        if traceur==0:
            core_genome.append(dict_final[liste_bacteries[0][:-3]][cles_gene[cle]])
            
        
taille_core_genome=len(core_genome)       
print(taille_core_genome)


#dict_final[liste_bacteries[genome][:-3]][cles_gene_genome[cle]]
