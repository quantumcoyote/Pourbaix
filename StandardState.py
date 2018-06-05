import os

output_path = 'scratch/'
species=[]
H=[]
G=[]

for filename in os.listdir(output_path):
    if "Energy" in filename:
        file=open("scratch/"+str(filename),"r")
        for line in file:
            species.append(filename.strip(".dat").replace("Energy_",""))
            tmp = line.strip('\n').split()
            H.append(tmp[0])
            G.append(tmp[1])
        file.close()

for filename in os.listdir(output_path):
    if "ReactionMatrix.dat" in filename:
        file=open("scratch/"+str(filename),"r")
        for line in file:
            tmp = line.strip('\n').split()
            for i in range(0,len(tmp)):
                print(float(tmp[i])*float(H[i]))
                print(float(tmp[i])*float(H[i]))
        file.close()