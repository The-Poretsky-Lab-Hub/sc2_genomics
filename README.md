# sc2_genomics
Scripts and such supporting SARS-CoV-2 genomic analyses

translate_nt_substitution.R: Input a SNP (or column of SNPs) and recieve the amino acid change it creates. Returns "Synonymous" if the SNP is in a protein-coding region but does not confer an AA change. Returns "Non-coding" if the SNP is in an untranslated region. This script uses the SARS-CoV-2 reference genome, but should work for any reference with non-overlapping ORFs. 
