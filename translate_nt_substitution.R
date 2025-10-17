library(tidyverse)
library(Biostrings)


#Load a reference annotation in gff3 format (.gff or .gff3)
read_gff3 <- function(gff) {
  read_tsv(gff,
    comment = "#",
    col_names = c("seqid", "source", "type", "start", "end",
                  "score", "strand", "phase", "attrb"),
    col_types = "ccciicccc"
  )
}

#Import SARS-CoV-2 example
annot <- read_gff3("NC_045512_Hu-1.gff")

#Also import SARS-CoV-2 reference FASTA
ref <- readDNAStringSet("NC_045512.2.fna")[[1]]

#Find start and end positions of genes in reference
get_gene_positions <- function(annot) {
  filter(annot, type == "gene") %>%
    select(start, end, attrb) %>%
    mutate(
      gene = str_remove(
        map(
          str_split(
            .$attrb,
            pattern = ";"),
           3),
        pattern = "Name=")) %>%
    select(gene, start, end)
}

#Get the gene positions of the example reference:
gene_positions <- get_gene_positions(annot)

###Use these examples for testing###
example_snps <- c("C26060T",
                  "C241T",
                  "T670G",
                  "C897A",
                  "C2790T",
                  "G21624C",
                  "C21711T",
                  "T22928C",
                  "T22942A",
                  "G22992A",
                  "T23005A",
                  "C25000T")
####################################

#Translate a character vector of SNPs, where
#gene_positions should be a dataframe created by using get_gene_positions() on
#an imported GFF3-format annotation file, and
#ref is a DNAString imported with Biostrings::readDNAStringSet()
translate_snps <- function(snps, gene_positions, ref) {
  snp_info <- tibble(snp = snps) %>%
    mutate(loc = as.integer(str_sub(snps, 2, -2))) %>%
    mutate(gene = map_chr(loc, function(l) {
      grange <- gene_positions %>%
        filter(l >= start & l <= end)
      if (nrow(grange) == 0) {
        "None"
      } else {
        grange$gene[1]
      }
    })) %>%
    left_join(gene_positions, by = "gene") %>%
    mutate(codon_rf = (loc - start) %% 3) %>%
    mutate(aa_loc = (loc - codon_rf - start) / 3 + 1) %>%
    mutate(ref_codon = pmap(list(loc, codon_rf), function(l, c) {
      s <- l - c
      e <- s + 2
      if (is.na(s)) {
        DNAString("")
      } else {
        (subseq(ref, start = s, end = e))
      }
    })) %>%
    mutate(ref_aa = map(ref_codon, translate)) %>%
    mutate(new_codon = pmap(list(ref_codon, codon_rf, snp), function(r, c, s) {
      if (length(r) == 0) {
        DNAString("")
      } else {
        new_codon <- r
        new_codon[c + 1] <- str_sub(s, start = -1)
        return(new_codon)
      }
    })) %>%
    mutate(new_aa = map(new_codon, translate)) %>%
    mutate(aa_sub = pmap(list(ref_aa, aa_loc, new_aa), function(r, a, n) {
      if (r != n) {
        paste0(as.character(r), as.character(a), as.character(n))
      } else if (length(r) == 0) {
        "Non-coding"
      } else {
        "Synonymous"
      }
    })) %>%
    select(aa_sub) %>%
    unlist() %>%
    as.character()
  return(snp_info)
}

test <- data.frame(example_snps)
test_translated <- test %>%
  mutate(translated = translate_snps(example_snps, gene_positions, ref))
