import os
import math
import sys

output_path = 'scratch/'
species=[]
concentrations=[]
H=[]
G=[]
T=float(sys.argv[1])
P=float(sys.argv[2])
Reactions_G=[]

#
# Obtain the concentrations
#

for filename in os.listdir(output_path):
    if "concentrations.dat" in filename:
        file=open("scratch/"+str(filename),"r")
        for line in file:
            tmp=line.split()
            concentrations.append(tmp[1])

#
# Obtain the Enthalpy and Gibbs Free Energy
#

for filename in os.listdir(output_path):
    if "Energy" in filename:
        file=open("scratch/"+str(filename),"r")
        for line in file:
            species.append(filename.strip(".dat").replace("Energy_",""))
            tmp = line.strip('\n').split()
            H.append(tmp[0])
            G.append(tmp[1])
        file.close()

#
# Apply the reaction matrix
#

plot=open("scratch/plot.dat","a")

file_labels = open("Reactions.dat","r")
for filename in os.listdir(output_path):
    if "ReactionMatrix.dat" in filename:
        file=open("scratch/"+str(filename),"r")
        for line in file:
            tmp = line.strip('\n').split()
            stdstate=1.0
            Tot_G=0.0
            Tot_H=0.0
            for i in range(0,len(tmp)):
                if(float(tmp[i])) < 0.0:
                    stdstate = stdstate * (1/(float(concentrations[i])**(-1.0*float(tmp[i]))))
                if(float(tmp[i])) > 0.0:
                    stdstate = stdstate * ((float(concentrations[i])**(float(tmp[i]))))
                Tot_H=Tot_H+(float(tmp[i])*float(H[i]))
                Tot_G=Tot_G+(float(tmp[i])*float(G[i]))
            stdstate=0.008314462175*T*math.log(stdstate)
            Tot_G=Tot_G+stdstate
            Reactions_G.append(Tot_G)
            #print(str(file_labels.readline()).strip('\n'),' dG=',Tot_G,' kJ/mol')
            plot.write(str(T)+' '+str(sys.argv[3])+' '+str(Tot_G)+'\n')
        file.close()
file_labels.close()

#print('')
#print('Minimum Energy Reaction')
#print(Reactions_G.index(min(Reactions_G))+1)
