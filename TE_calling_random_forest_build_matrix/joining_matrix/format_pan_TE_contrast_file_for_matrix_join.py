B73_bed = []
with open("B73_allContrasts_filteredTE_13Jan19.bed") as bedfile:
    for TE_line in bedfile:
        # Split each line.
        TE_line = TE_line.strip().split()
        TE_name = TE_line[3].split(";")[0]
        TE_coord = TE_line[0] +":"+ TE_line[1] + "-" + TE_line[2]
        TE_bed_fmt = TE_name + "\t" + TE_coord
        B73_bed.append(TE_bed_fmt)
        
with open('B73_TE_fmt.bed', 'w') as f:
    for item in B73_bed:
        f.write("%s\n" % item)
        
W22_bed = []
with open("W22_allContrasts_filteredTE_13Jan19.bed") as bedfile:
    for TE_line in bedfile:
        # Split each line.
        TE_line = TE_line.strip().split()
        TE_name = TE_line[3].split(";")[0]
        TE_coord = TE_line[0] +":"+ TE_line[1] + "-" + TE_line[2]
        TE_bed_fmt = TE_name + "\t" + TE_coord
        W22_bed.append(TE_bed_fmt)
        
with open('W22_TE_fmt.bed', 'w') as f:
    for item in W22_bed:
        f.write("%s\n" % item)
        

Mo17_bed = []
with open("Mo17_allContrasts_filteredTE_13Jan19.bed") as bedfile:
    for TE_line in bedfile:
        # Split each line.
        TE_line = TE_line.strip().split()
        TE_name = TE_line[3].split(";")[0]
        TE_coord = TE_line[0] +":"+ TE_line[1] + "-" + TE_line[2]
        TE_bed_fmt = TE_name + "\t" + TE_coord
        Mo17_bed.append(TE_bed_fmt)
        
with open('Mo17_TE_fmt.bed', 'w') as f:
    for item in Mo17_bed:
        f.write("%s\n" % item)
        
        
PH207_bed = []
with open("PH207_allContrasts_filteredTE_13Jan19.bed") as bedfile:
    for TE_line in bedfile:
        # Split each line.
        TE_line = TE_line.strip().split()
        TE_name = TE_line[3].split(";")[0]
        TE_coord = TE_line[0] +":"+ TE_line[1] + "-" + TE_line[2]
        TE_bed_fmt = TE_name + "\t" + TE_coord
        PH207_bed.append(TE_bed_fmt)
        
with open('PH207_TE_fmt.bed', 'w') as f:
    for item in PH207_bed:
        f.write("%s\n" % item)        
