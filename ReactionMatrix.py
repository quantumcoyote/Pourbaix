import os
import numpy as np
import sys

output_path = 'outputs'
species=[]

for filename in os.listdir(output_path):
    species.append(filename.strip('.out'))

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

output = open("scratch/concentrations.dat","w")
for i in range(0,len(species)):
    if str(species[i])=='H2O':
      output.write(str(species[i]) + '  55.6 \n')
    if str(species[i])=='OH':
      output.write(str(species[i]) + '  '+str(sys.argv[1])+' \n')
    if str(species[i])!='OH' and str(species[i])!='H2O':
      output.write(str(species[i])+'  1.0 \n')

output.close()
