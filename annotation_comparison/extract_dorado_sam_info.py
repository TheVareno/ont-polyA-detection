
import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import os 
import re
import pysam # type: ignore


#pd.set_option('display.max_columns', None)
#pd.set_option("max_rows", None)

def set_working_directory():
    target_dir = '/home/hi68ren/Dokumente/ProjektModul/Implementation/main_data/ecoli/annotation_comparison'
    if os.getcwd != target_dir : 
        os.chdir(target_dir)


set_working_directory()

sam_data = pysam.AlignmentFile("dorado_output/ecoli_dorado_calls_polya.sam", "r")

print(sam_data)

#sam_data = pd.read_csv('dorado_output/ecoli_dorado_calls_polya.sam', sep='\t', comment='@', header=None)

# Filter rows containing 'pt:i:number' pattern
#rows_with_pti = sam_data[sam_data.apply(lambda row: row.astype(str).str.contains(r'\bpt:i\b').any(), axis=1)]


# rows_with_pti.to_csv('rows_with_pt.csv', index=True)

