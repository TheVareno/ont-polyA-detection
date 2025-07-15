
"""
Author : Hadi Vareno 
Email : mohammad.noori.vareno@uni-jena.de
Github:  https://github.com/TheVareno
"""

from read5.Reader import read # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import ont_fast5_api.fast5_interface     # type: ignore
import random as rnd
from csv import writer
import os 
import pandas as pd # type: ignore
import seaborn as sns # type: ignore


def set_directory():
    if os.getcwd() != '/home/hi68ren/Dokumente/ProjektModul/Implementation/main_data/ivt' : 
        os.chdir('/home/hi68ren/Dokumente/ProjektModul/Implementation/main_data/ivt')
    else : 
        pass


def read_ivt_fast5(): 
    ivt_reads = read('./FAO07649_853a2ebb_0.fast5') 
    return ivt_reads


def get_random_read_ids(input_tsv: str, output_file: str, n: int = 100):
    
    # Read the Nanopolish TSV file into a DataFrame
    df = pd.read_csv(input_tsv, sep='\t')

    sampled_read_ids = df['readname'].sample(n=n, random_state=42)

    with open(output_file, 'w') as file:
        for read_id in sampled_read_ids:
            file.write(f"{read_id}\n")


def process_read_id(input_file: str, used_file: str):
    
    if os.stat(input_file).st_size == 0:
        print("No read IDs left to process.")
        return None

    # Read remaining read_ids 
    with open(input_file, 'r') as file:
        read_ids = file.readlines()

    current_read_id = read_ids[0].strip()

    with open(used_file, 'a') as used_file_obj:
        used_file_obj.write(f"{current_read_id}\n")

    # Remove the processed read_id from the input file
    with open(input_file, 'w') as file:
        file.writelines(read_ids[1:])  # Write back the remaining IDs

    return current_read_id


def add_coordinate(read_id:str, start:int, end:int):
    
    new_row = [read_id, start, end]
    
    new_row_str = ','.join(str(e) for e in new_row)
    
    with open('./annotations/coordinates.txt', 'a', newline='') as coord_file:
        coord_file.write(f'{new_row_str}\n')

def plot_signals(read_signal, next_read_id): 
    
    plt.figure(figsize=(20, 8))
    plt.plot(read_signal, color='#4cbb17', linewidth=0.8)
    plt.title(f'Raw Signal of read : {next_read_id}')
    plt.xlabel('Data Points')
    plt.ylabel('Normalized Signal Intensity')
    plt.grid(True)
    plt.pause(10)

    # next plot with boundaries
    pA_start = int(input("polyA start coordinate: "))
    pA_end = int(input("polyA end coordinate: "))

    pA_coordinate = (pA_start, pA_end) 

    add_coordinate(read_id=next_read_id, start= pA_start, end= pA_end)

    plt.figure(figsize=(20, 8))
    plt.plot(read_signal, color='#4cbb17', linewidth=0.8)

    for i in pA_coordinate:
        plt.axvline(x = i, color='#ff2c2c', ls='--', lw=1.0)

    plt.title(f'Raw Signal of read : {next_read_id}', fontweight='bold', fontsize=16)
    plt.xlabel('Data Points', fontweight='bold')
    plt.ylabel('Normalized Signal Intensity', fontweight='bold')
    plt.grid(True)
    plt.savefig(f'./annotations/plots/{next_read_id}.png', bbox_inches = 'tight')
    plt.show()


def set_style():     
    sns.set_style('darkgrid')




def main(): 
    
    set_directory()
    
    set_style()
    
    # get from nanopolish 
    
    # get_random_read_ids('./nanopolish/ivt_polyA_nanopolish_output.tsv', './annotations/to_annotate_read_ids.txt', 100)
    
    ivt = read_ivt_fast5() 
    next_read_id = process_read_id('./annotations/to_annotate_read_ids.txt', './annotations/used_read_ids.txt')
    
    if next_read_id:
        print(f"Processed read ID: {next_read_id}")
 
    read_signal = ivt.getZNormSignal(next_read_id)

    plot_signals(read_signal, next_read_id)


    



if __name__ == '__main__' : 
    main()


    
    