in_name = "zTest.tad"
out_name = "zTestMerged.tad"


with open(in_name, "r") as infile,  open(out_name, "w") as outfile:
    for line in infile:
        if line.startswith("@"):
            outfile.write(line);
        elif line.startswith("scaffold_"):
            split = line.split()
            if(len(split)!=27):
                print "ERROR!!:", line, split
                
            
            if len(split[2])==1:
                counts = map(int, split[3:])
                total = str(sum(counts))
            else:
                counts = [ map(int, c.split(",")) for c in split[3:] ]
                each_total = [0]*len(split[2])
                for i in counts:
                    each_total =  [x+y for x,y in zip(each_total, i)]
                total = ",".join(map(str,each_total))
            
            labels = split[0:3]
            new_line = "\t".join(labels) + "\n"
            outfile.write(new_line)
