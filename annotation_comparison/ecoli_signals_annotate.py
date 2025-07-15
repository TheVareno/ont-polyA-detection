
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


if os.getcwd() != '/home/hi68ren/Dokumente/ProjektModul/Implementation/main_data/ecoli/annotation_comparison' : 
    os.chdir('/home/hi68ren/Dokumente/ProjektModul/Implementation/main_data/ecoli/annotation_comparison')
else : 
    pass


ecoli_ivv = read("first_fast5/FAX28269_36c48ee6_b042d0cd_0.fast5") 


def get_random_read():
    all_reads = ecoli_ivv.getReads()
    return rnd.choice(all_reads)


# append to table
def add_coordinate(read_id:str, start:int, end:int, alt_start:any, alt_end:any):
    
    alt_start = None if alt_start == 'n' else alt_start
    alt_end = None if alt_end == 'n' else alt_end

    new_row = [read_id, start, end, alt_start, alt_end]

    with open('annotations/coordinates.csv', 'a', newline='') as csv:
        writer_obj = writer(csv)
        writer_obj.writerow(new_row)
       

# add to used read ids 
def add_readid(read_id:str):
    with open('annotations/used_readids.txt', 'a') as txt:
        txt.write(read_id + '\n') 
           

#! read id random [1, 3999)
read_id = get_random_read()
print(f'The current chosen read id : {read_id}')


#! main signal  
read_signal = ecoli_ivv.getZNormSignal(read_id)

# show first plot     
plt.figure(figsize=(20, 8))
plt.plot(read_signal, color='green', linewidth=0.8)
plt.title(f'Raw Signal of read : {read_id}')
plt.xlabel('X Values')
plt.ylabel('Signal Intensity')
plt.grid(True)
plt.show()


# next plot with boundaries
pA_start = int(input("polyA start coordinate: "))
pA_alt_start = input("poly alternative start coordinate: ")

pA_end = int(input("polyA end coordinate: "))
pA_alt_end = input("polya alternative end coordinate:")

pA_coordinate = (pA_start, pA_end) 

# append to csv file as new row
add_coordinate(read_id=read_id, start= pA_start, end= pA_end, alt_start=pA_alt_start, alt_end=pA_alt_end)
add_readid(read_id=read_id)

# second plot with bounderies 
plt.figure(figsize=(20, 8))
plt.plot(read_signal, color='green', linewidth=0.8)

for i in pA_coordinate:
    plt.axvline(x = i, color='red', ls='--', lw=0.8)

plt.title(f'Raw Signal of read : {read_id}')
plt.xlabel('X Values')
plt.ylabel('Signal Intensity')
plt.grid(True)

plt.savefig(f'plots/{read_id}.png', bbox_inches = 'tight')

plt.show()
    

