# gRNA
Guide RNA desigh for CRISPR/Cas9 genome editing


This Python package identifies candidate guide RNA for CRISPR/Cas9 genome editing based on the NGG structure of gRNA, then aligns those sequences to the reference genome and scores them based on the alignment. Next, the package reads the provided gff3 file and annotates the gRNA based on this information. The final output includes the location, strand, annotation, and score of each candidate gRAN.


Users' Guide

Double-click on the gRNA program to run it, and wait for the software interface to appear as shown below (Figure 1):

![image](https://github.com/user-attachments/assets/d3a4e114-693b-4091-9dea-aee8036dd4d3)


1. lick the "Target Sequence File" button in the software, and select a Fasta format file, such as the test file "target_seq.Fasta" located in the test_data folder (as shown in Figure 2).

2. Click the "Reference Sequence File" button in the software, and select a Fasta format file, such as the test file "reference_seq.Fasta" located in the test_data folder (as shown in Figure 2).

3. Click the "Target Sequence Annotation" button in the software, and select a GFF3 format file, such as the test file "target_gene.GFF3" located in the test_data folder (as shown in Figure 2). (If you don’t have the file, you can skip this step; the annotation output will be marked as "NA").

4. Click the "Select Output File" button, choose the output file location and name. For example, name the result as "output" and save it to the test_data folder (as shown in Figure 2).


