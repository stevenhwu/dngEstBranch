in_name = "bwa_filtered.tad"
out_name = "bwa_filtered_merged.tad"
sum_file = "bwa_filtered_sum.tad"

with open(in_name, "r") as infile,  open(out_name, "w") as outfile, open(sum_file, "w") as sumfile:
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
                sum_count = total
            else:
                counts = [ map(int, c.split(",")) for c in split[3:] ]
                each_total = [0]*len(split[2])
                for i in counts:
                    each_total =  [x+y for x,y in zip(each_total, i)]
                total = ",".join(map(str,each_total))
                sum_count = sum(each_total)
                
            new_line = "\t".join(split[0:3]) + "\t" + total + "\n"
            outfile.write(new_line)
            
            sum_line = "\t".join(split[0:2]) + "\t" + str(sum_count) + "\n"
            sumfile.write(sum_line);
