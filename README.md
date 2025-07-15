
# Segmentation of ONT Raw Signal of Data for Detecting Poly(A) Region  

## Overview 

- Identifying poly Adenosine (poly(A)) regions within raw signal data from Oxford Nanopore Technology (ONT) direct RNA sequencing (dRNA-seq) is crucial for
  understanding RNA biology, especially poly(A) tail length, which plays a regulatory role in gene expression. 

- In this work, we first collect 100 random reads to manually annotate the poly(A) start and end positions. Then we develop a segmentation model based on a Hidden Markov Model (HMM). 

- We next evaluate the performance of our model together with two existing tools, TailfindR and Nanopolish polya, to find out how accurate each tool is in detecting the poly(A)region versus human annotations, which are considered as ground truth during this study.

- While TailfindR outperformed Nanopolish polya in terms of segmentation accuracy, our own model performed the best consistently, particularly on difficult reads where other tools struggled.


## Data Availability

Due to privacy constraints and the large volume of raw data, **sequencing reads and processed datasets are not included** in this repository.

If you are interested in data access or collaboration, please contact me directly.


## Contact

Feel free to reach out with any questions or collaboration ideas:

- 📧 mohammad.noori.vareno@uni-jena.de  
- 📧 hadivareno@gmail.com
