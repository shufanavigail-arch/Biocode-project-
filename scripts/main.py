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

  
        
            
#open files
codon_file = open('data/codons.txt', 'r')
human_file = open('data/homosapien_genes.txt', 'r')
ecoli_file = open('data/e.coli_genes.txt', 'r')

#create list of codons from codon file
codon_list = []
for line in codon_file:
    line = line.strip('\n')
    codon_list.append(line)

#call vector function for both human and ecoli files and print results
print(vector(human_file, codon_list))
print('\n')
print(vector(ecoli_file, codon_list))

#close files
codon_file.close()
human_file.close()
ecoli_file.close()