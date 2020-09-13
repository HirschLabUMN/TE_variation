#script to remove duplicate genes in files with TE location to genes in base pairs                                                                                                                      
#if a TE is present more than 1x, save smallest (closest to gene) value                                                                                                                                 

import pandas as pd
import os
import sys
#will need pandas for data frame manipulation                                                                                                                                                           

input_file = sys.argv[1]
output_file = sys.argv[2]

orig_file = pd.read_csv(input_file, sep = "\t", header = None)
df_orig_file = pd.DataFrame(orig_file)
df_orig_file.columns = ["TE_name","location_10kbspan","distance_to_gene"]

#df.groupby('A', as_index=False).max()                                                                                                                                                                  
#group by TE name and location and save only the smallest value                                                                                                                                         
newdf_nodup = df_orig_file.groupby(['TE_name','location_10kbspan']).min()
#print(newdf_nodup)                                                                                                                                                                                     
newdf_nodup.to_csv(output_file, sep = '\t', header = True, index = True)


