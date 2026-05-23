import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn import datasets
from mpl_toolkits.mplot3d import Axes3D

def vector(file, list):
    codon_dict = {}
    database = []
    code_part = ""
    for line in file:
        line = line.strip()
        if line[0] != '>':
            code_part += line
        if line[0] == '>' and code_part != "":
            if len(code_part) % 3 != 0:
                code_part = code_part[:len(code_part)-len(code_part)%3]
            for item in list:
                codon_dict[item] = 0
            freq_list = []
            RNA_seq = code_part.replace('T', 'U')
            for i in range(0, len(RNA_seq), 3):
                codon = RNA_seq[i:i+3]
                if codon in list:
                    codon_dict[codon] += 1
            for item in list:
                codon_dict[item] = codon_dict[item]/(len(RNA_seq)/3)
                codon_dict[item] = round(codon_dict[item], 5)
                freq_list.append(codon_dict[item])
            database.append(freq_list)
            code_part = ""
        

    if len(code_part) % 3 != 0:
        code_part = code_part[:len(code_part)-len(code_part)%3]
    for item in list:
        codon_dict[item] = 0
    freq_list = []
    RNA_seq = code_part.replace('T', 'U')
    for i in range(0, len(RNA_seq), 3):
        codon = RNA_seq[i:i+3]
        if codon in list:
            codon_dict[codon] += 1  
    for item in list:
        codon_dict[item] = codon_dict[item]/(len(RNA_seq)/3)
        codon_dict[item] = round(codon_dict[item], 5)
        freq_list.append(codon_dict[item])
    database.append(freq_list)

    return database

  
        
            

codon_file = open('data/codons.txt', 'r')
human_file = open('data/homosapien_genes.txt', 'r')
ecoli_file = open('data/e.coli_genes.txt', 'r')

codon_list = []
for line in codon_file:
    line = line.strip('\n')
    codon_list.append(line)

human_data = vector(human_file, codon_list)
ecoli_data = vector(ecoli_file, codon_list)
#print(human_data)
print('\n')
#print(ecoli_data)

data = human_data + ecoli_data
print(data)

human_len = len(human_data)
ecoli_len = len(ecoli_data)
#print(human_len)
#print(ecoli_len)

pca= PCA(n_components=50)



data_pca = pca.fit_transform(data)



pc_x = 1
pc_y = 5
pc_z =8
 
x = data_pca[:, pc_x-1]
y = data_pca[:, pc_y-1]
z = data_pca[:, pc_z-1]


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(x[:human_len], y[:human_len], c='blue', edgecolor='k')
ax.scatter(x[human_len:], y[human_len:], c='red', edgecolor='k')

ax.scatter(x[:human_len], y[:human_len], z[:human_len], c='blue', edgecolor='k')
ax.scatter(x[human_len:], y[human_len:], z[human_len:], c='red', edgecolor='k')



ax.set_xlabel('Principal Component 1')
ax.set_ylabel('Principal Component 2')
ax.set_zlabel('Principal Component 3')
ax.set_title('PCA of Codon Usage')
plt.show()


codon_file.close()
human_file.close()
ecoli_file.close()