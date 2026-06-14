
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn import datasets
from mpl_toolkits.mplot3d import Axes3D
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import train_test_split
import random
from sklearn.neighbors import KNeighborsClassifier



#function that creates a vector of codon frequencies for each gene in the file
def vector(file, list):

    codon_dict = {}    # Dictionary to store codon frequencies
    database = []    # Stores frequency vectors for all genes
    code_part = ""

    for line in file:
        line = line.strip()
        if line[0] != '>':
            code_part += line
        if line[0] == '>' and code_part != "":
            if len(code_part) % 3 != 0:
                code_part = code_part[:len(code_part)-len(code_part)%3]   #makes sure the code part can be divided by 3
            for item in list:
                codon_dict[item] = 0
            freq_list = []   #store codon frequencies for current gene
            RNA_seq = code_part.replace('T', 'U')
            for i in range(0, len(RNA_seq), 3):
                codon = RNA_seq[i:i+3]
                if codon in list:
                    codon_dict[codon] += 1
            for item in list:
                codon_dict[item] = codon_dict[item]/(len(RNA_seq)/3)
                codon_dict[item] = round(codon_dict[item], 5)
                freq_list.append(codon_dict[item])
            database.append(freq_list)  #add frequency vector for current gene to database
            code_part = ""  #reset code part for next gene
        
#same process for the last gene in the file
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

def confirm(start, end, size):
    big_list = []
    for i in range(size):
        list = []
        for j in range(64):
            list.append(0)
        for u in range(0,64):
            list[u] = random.randrange(0, 10)
        for o in range(start, end):
            list[o] = random.randrange(10,50)
        total = sum(list)
        list = [x/total for x in list]
            # print(list)
        big_list.append(list)
    return big_list
        

check = True
#open files
codon_file = open('data/codons.txt', 'r')
human_file = open('data/homosapien_genes1.txt', 'r')
ecoli_file = open('data/e.coli_genes1.txt', 'r')

#create list of codons from codon file
codon_list = []
for line in codon_file:
    line = line.strip('\n')
    codon_list.append(line)

#call vector function for both human and ecoli files and print results
#print(vector(human_file, codon_list))
if check == False:
    human_data = vector(human_file, codon_list)
    ecoli_data = vector(ecoli_file, codon_list)
else:
    human_data = confirm(0, 32, 1000)
    ecoli_data = confirm(32, 64, 1000)
#print('\n')


data = human_data + ecoli_data
#print(data)

human_len = len(human_data)
ecoli_len = len(ecoli_data)
#print(human_len)
#print(ecoli_len)

pca= PCA(n_components=50)



data_pca = pca.fit_transform(data)



pc_x = 1
pc_y = 2
pc_z = 3
 
x = data_pca[:, pc_x-1]
y = data_pca[:, pc_y-1]
z = data_pca[:, pc_z-1]


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#ax.scatter(x[:human_len], y[:human_len], c='blue', edgecolor='k')
#ax.scatter(x[human_len:], y[human_len:], c='red', edgecolor='k')

ax.scatter(x[:human_len], y[:human_len], z[:human_len], c='blue', edgecolor='k')
ax.scatter(x[human_len:], y[human_len:], z[human_len:], c='red', edgecolor='k')



ax.set_xlabel('Principal Component %d' % pc_x)
ax.set_ylabel('Principal Component %d' % pc_y)
ax.set_zlabel('Principal Component %d' % pc_z)
ax.set_title('PCA of Codon Usage')
plt.show()

labels = []
for i in range(int(human_len)):
    labels.append('human')
for i in range(int(ecoli_len)):
    labels.append('ecoli')

train_data, test_data, train_labels, test_labels = train_test_split(data, labels, test_size=0.2, random_state=42)
accurate_count = 0


lda = LDA(n_components=1)
lda.fit(train_data, train_labels)

predictions = lda.predict(test_data)


for i in range (len(predictions)):
    if predictions[i] == test_labels[i]:
        accurate_count += 1
sim_per = accurate_count/len(test_labels)*100

print('The accuracy percentage between the predictions and the test labels is: ' + str(sim_per) + '%')
print('The predictions are:' )
print(predictions)
print('The actual labels are:')
print(test_labels)


#close files
codon_file.close()
human_file.close()
ecoli_file.close()