import os
import numpy as np

output_path = 'outputs'
species=[]

for filename in os.listdir(output_path):
    species.append(filename.strip('.log'))

file = open("Reactions.dat","r")
output = open("scratch/ReactionMatrix.dat","w")
for line in file:
    tmp=line.strip('\n').split()
    reaction=np.zeros(len(species))
    sign = -1
    for i in range(0,len(tmp)):
        for j in range(0,len(species)):
            if str(tmp[i]) == str(species[j]):
                reaction[j]=sign*int(tmp[i-1])
            if str(tmp[i]) == '>>>':
                sign=1
    tmp2=''
    tmp2=str(reaction).strip('[').strip(']')
    output.write(str(tmp2)+'\n')

file.close()
output.close()