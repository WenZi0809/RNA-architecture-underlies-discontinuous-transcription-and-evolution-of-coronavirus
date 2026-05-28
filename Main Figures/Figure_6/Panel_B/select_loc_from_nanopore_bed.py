from sys import argv

pyname, bed_file, sam_file, out_sam_file = argv

#reads_id_test = []
reads_id = set()

with open(bed_file,"r")as bed:
    for line in bed:
        line = line.strip().split()
        name = line[0]
        map_loc = line[-1].split("-")
        if (len(map_loc) < 3):
            continue
        if ((map_loc[0] == 'leader') & (map_loc[1] == 'ORF1a') & (map_loc[2] == '>N')):
            reads_id.add(name)
            #reads_id_test.append(name)

#with open("name_test.txt","w")as f:
#    for i in reads_id_test:
#        f.write(i+"\n")

sam = open(sam_file,"r")
out_sam = open(out_sam_file,"w")

for line in sam:
    if line.startswith("@"):
        out_sam.write(line)
    else:
        line_split = line.strip().split()
        seq_id = line_split[0]
        if seq_id in reads_id:
            out_sam.write(line)

sam.close()
out_sam.close()