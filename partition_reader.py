fichier = "partitions/partition_7.txt"

with open(fichier) as file:
    lines = file.readlines()
    output = []
    for line in lines:
        line_list = line.split()
        output.append(list(map(int,line_list)))
    P=output

for part in P:
    n = len(part)
    for i in range(n):
        part.append(841-part[i])
        
with open("partitions/partition_7.txt","w") as file:
    for part in P:
        part.sort()
        for i in part:
            file.write(str(i)+" ")
        file.write("\n")