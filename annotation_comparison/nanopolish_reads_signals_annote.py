"""
status : paths revised. 

description : 
get nanopolish estimated read_ids. 
since nnps performs mapping in first hand then polyA length estimation. 

"""

import numpy as np  # type: ignore
import pandas as pd # type: ignore

import matplotlib.pyplot as plt 
import ont_fast5_api.fast5_interface    
import random as rnd
from csv import writer
import os 
from read5.Reader import read



def setup_working_directory():
    target_dir = '/home/hi68ren/Dokumente/ProjektModul/Implementation/main_data/ecoli/annotation_comparison'
    if os.getcwd() != target_dir:
        os.chdir(target_dir)



def load_data():
    nnps_polya_df = pd.read_csv("nanopolish_polya/ecoli_polya_all_reads.tsv", sep='\t')
    with open("annotations/used_readids.txt") as readids:
        used_read_ids = [line.rstrip() for line in readids]
    return nnps_polya_df, used_read_ids



def categorize_read_ids(nnps_read_ids, used_read_ids):
    cmn_read_ids = [read_id for read_id in nnps_read_ids if read_id in used_read_ids]
    nnps_uniq_read_ids = [read_id for read_id in nnps_read_ids if read_id not in used_read_ids]
    
    print(f"Common read IDs: {len(cmn_read_ids)}")
    print(f"Unique Nanopolish read IDs: {len(nnps_uniq_read_ids)}")
    
    return cmn_read_ids, nnps_uniq_read_ids


# main fast5 
# ecoli_ivv = read("FAX28269_36c48ee6_b042d0cd_0.fast5") 


def get_random_read_nnps(nnps_uniq_read_ids):
    return rnd.choice(nnps_uniq_read_ids)


def get_random_read(ecoli_ivv):
    all_reads = ecoli_ivv.getReads()
    return rnd.choice(all_reads)


# append to table
def add_coordinate(read_id:str, start:int, end:int, alt_start:any, alt_end:any):
    
    alt_start = None if alt_start == 'n' else alt_start
    alt_end = None if alt_end == 'n' else alt_end

    new_row = [read_id, start, end, alt_start, alt_end]

    with open('coordinates_nnps.csv', 'a', newline='') as csv:
        writer_obj = writer(csv)
        writer_obj.writerow(new_row)
       

# add to used read ids 
def add_readid(read_id:str):
    with open('used_readids_nnps.txt', 'a') as txt:
        txt.write(read_id + '\n') 
             

def plot_signal(read_signal, read_id, pA_coordinate=None, save=False):
    plt.figure(figsize=(20, 8))
    plt.plot(read_signal, color='green', linewidth=0.8)
    
    if pA_coordinate:
        for i in pA_coordinate:
            plt.axvline(x=i, color='red', ls='--', lw=0.8)
    
    plt.title(f'Raw Signal of read: {read_id}')
    plt.xlabel('X Values')
    plt.ylabel('Signal Intensity')
    plt.grid(True)
    
    if save:
        plt.savefig(f'plots_nnps/{read_id}.png', bbox_inches='tight')
    
    plt.show()


def convert_ls_to_txt(ls:list , path_to_file:str) :
    with open(path_to_file, 'a') as txt: 
        for item in ls: 
            txt.write(item + '\n')




def main():
    setup_working_directory()
    
    nnps_polya_df, used_read_ids = load_data()
    nnps_read_ids = nnps_polya_df['readname']
    cmn_read_ids, nnps_uniq_read_ids = categorize_read_ids(nnps_read_ids, used_read_ids)
    
    print('success')

    #ecoli_ivv = read("FAX28269_36c48ee6_b042d0cd_0.fast5")
    #read_id = get_random_read_nnps(nnps_uniq_read_ids)
    
    #print(f'The current chosen read ID: {read_id}')
    #read_signal = ecoli_ivv.getZNormSignal(read_id)
    
    #plot_signal(read_signal, read_id)
    
    #pA_start = int(input("polyA start coordinate: "))
    #pA_alt_start = input("polyA alternative start coordinate: ")
    #pA_end = int(input("polyA end coordinate: "))
    #pA_alt_end = input("polyA alternative end coordinate: ")
    
    #pA_coordinate = (pA_start, pA_end)
    #add_coordinate(read_id=read_id, start=pA_start, end=pA_end, alt_start=pA_alt_start, alt_end=pA_alt_end)
    #add_readid(read_id=read_id)
    
    #plot_signal(read_signal, read_id, pA_coordinate, save=True)



if __name__ == "__main__":
    main()