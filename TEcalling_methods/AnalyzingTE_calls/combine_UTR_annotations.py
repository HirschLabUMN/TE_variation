"""This script will take a 5' or 3' UTR annotation file and collapse any multiple UTR annotations that exist for the same transcript. Output will be gff format file with collapsed UTR annotations     
"""
import sys
import getopt
import operator
import re
import tempfile
#usage etc setup taken from Meesh's script                                                                                                                                                              
def usage():
    print("""\n                                                                                                                                                                                         
        python3 edit_gene_name.py                                                                                                                                                                       
            -g or gff_file gff file with original version of names                                                                                                                                      
            -o or output_file  Name of new gff file that will be created                                                                                                                                
        \n""")
#add in 2 arguements to get lists of genes that do and don't overlap with                                                                                                                               
try:
    opts, args = getopt.getopt(sys.argv[1:], "g:o:h", ["gff_file=",
                                                         "output_file=",
                                                         "help"])

except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-g", "--gff_file"):
        gff_file = arg
    elif opt in ("-o", "--output_file"):
        output_file = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"

prev_utr_start = 0
prev_utr_end = 0
prev_utr_name = "name"
prev_utr_strand = ""
prev_utr_chr = 1
#get the name of gene in the last line of the file                                                                                                                                                      
f_read = open(gff_file, "r")
last_line = f_read.readlines()[-1]
f_read.close()
last_line = last_line.strip("\n")
last_line_name = last_line.split("\t")[8]
#get n lines in a file                                                                                                                                                                                  
n = len(open(gff_file).readlines( ))
print(last_line_name)
print(n)
count = 0

with open(output_file, "w") as output,\
     open(gff_file, "r") as orig_gff:
    for line in orig_gff:
        if line.startswith("#"):
            count += 1
            continue
        else:
            count += 1
            print(count)
            line = line.strip("\n")
            fields = line.split("\t")
            curr_utr_start = int(fields[3])
            curr_utr_end = int(fields[4])
            curr_utr_name = fields[8]
            curr_utr_strand = fields[6]
            curr_utr_chr = fields[0]
            #check if UTR gene name is the same                                                                                                                                                         
            if curr_utr_name == last_line_name: #if we're at the end of the file write out the prev gene info                                                                                           
                #this should save final lines that have multiple utr entries                                                                                                                            
                #print out prev_gene info                                                                                                                                                               
                print(str(prev_utr_chr)+"\t"+fields[1]+"\t"+fields[2]+"\t"+str(prev_utr_start)+"\t"+str(prev_utr_end)+"\t"+fields[5]+"\t"+prev_utr_strand+"\t"+fields[7]+"\t"+prev_utr_name, file = out\
put)
                prev_utr_start = int(fields[3])
                prev_utr_end = int(fields[4])
                prev_utr_name = fields[8]
                prev_utr_strand = fields[6]
                prev_utr_chr = fields[0]
                #count += 1                                                                                                                                                                             
            if count == n: #if we're at the last line                                                                                                                                                   
                #this isn't working quite as expected, but for at most 4 files I'll hand edit where necessary                                                                                           
                if curr_utr_name == prev_utr_name:
                #if the gene name is the same as the previous one the update previous UTR line end, and start if necessary                                                                              
                    if curr_utr_end > prev_utr_end:
                        prev_utr_end = int(fields[4])
                    if prev_utr_start > curr_utr_start:
                        #if for some reason the current UTR annotation starts before the previous one, change it                                                                                        
                        prev_utr_start = int(fields[3])
                    prev_utr_chr = fields[0]
                    prev_utr_strand = fields[6]
                    prev_utr_name = fields[8]
                    print(str(prev_utr_chr)+"\t"+fields[1]+"\t"+fields[2]+"\t"+str(prev_utr_start)+"\t"+str(prev_utr_end)+"\t"+fields[5]+"\t"+prev_utr_strand+"\t"+fields[7]+"\t"+prev_utr_name, file =\
 output)
                else:
                    print(line, file = output)

            if curr_utr_name == prev_utr_name:
                #if the gene name is the same as the previous one the update previous UTR line end, and start if necessary                                                                              
                if curr_utr_end > prev_utr_end:
                    prev_utr_end = int(fields[4])
                #else:                                                                                                                                                                                  
                #    continue                                                                                                                                                                           
                if prev_utr_start > curr_utr_start:
                    #if for some reason the current UTR annotation starts before the previous one, change it                                                                                            
                    prev_utr_start = int(fields[3])
                prev_utr_name = fields[8]
                #count += 1                                                                                                                                                                             
            else:
                #if the prev_utr_name = name (a check for the first line of the gff file), then write current line infor the prev_utr variables                                                         
                if prev_utr_name == "name":
                    prev_utr_start = int(fields[3])
                    prev_utr_end = int(fields[4])
                    prev_utr_name = fields[8]
                    prev_utr_strand = fields[6]
                    prev_utr_chr = fields[0]
                    #count += 1                                                                                                                                                                         
                else:
                    #if the name of the current utr isn't the same as the previous utr, write out prev utr information to new file                                                                      
                    print(str(prev_utr_chr)+"\t"+fields[1]+"\t"+fields[2]+"\t"+str(prev_utr_start)+"\t"+str(prev_utr_end)+"\t"+fields[5]+"\t"+prev_utr_strand+"\t"+fields[7]+"\t"+prev_utr_name, file =\
 output)
                    #output.write(prev_utr_name)                                                                                                                                                        
                    #update prev_gene infor with current line info                                                                                                                                      
                    prev_utr_start = int(fields[3])
                    prev_utr_end = int(fields[4])
                    prev_utr_name = fields[8]
                    prev_utr_strand = fields[6]
                    prev_utr_chr = fields[0]
                    #count += 1                                                                                                                                                                         


output.close()
orig_gff.close()

#end of script 
