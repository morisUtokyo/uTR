## Real data

The fasta file includes the following DNA sequences:

- SAND12(control,BAFME) hg38_dna range=chr8:118366813-118366928 in the 4th intron of SAMD12, pattern=(AAAAT)23

- SAND12(case,BAFME) Tandem repeat in the 4th intron of SAND12 found in Patient II-1 in family F6115 (Supplementary Figure 6, Ishiura, H. et al. Expansions of intronic TTTCA and TTTTA repeats in benign adult familial myoclonic epilepsy. Nat Genet 50, 581–590 (2018))  pattern=(ATTTT)221(ATTTC)221(ATTTT)82

- SAND12(case,BAFME) Tandem repeat in the 4th intron of SAND12 found in (Supplementary Figure 6, Ishiura, H. et al. Expansions of intronic TTTCA and TTTTA repeats in benign adult familial myoclonic epilepsy. Nat Genet 50, 581–590 (2018))  pattern=(ATTTT)613(ATTTC)320(ATTTT)5(ATTTC)130

- RFC1(control,CANVAS) hg38_dna range=chr4:39348425-39348483 pattern=(AAAAG)11

- KAZN(control) hg38_dna range=chr1:14883297-14883426 pattern=(AAAG)6(AG)11(AAAG)20

- ZNF37A(control) hg38_dna range=chr10:38112731-38112826 pattern=(CTTTT)12(CTTGT)3(CTTTT)2

To apply uTR, use

    bash test.sh
