
from ont_fast5_api.fast5_interface import get_fast5_file

def print_all_raw_data():
    
    fast5_filepath = "/Users/hadivareno/Documents/Projects/ProjektModul/Implementation/main_data/annotation_comparison/ecoli/first_fast5/FAX28269_36c48ee6_b042d0cd_0.fast5"
    
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read in f5.get_reads():
            raw_data = read.get_raw_data()
            print(read.read_id, raw_data)


print_all_raw_data()