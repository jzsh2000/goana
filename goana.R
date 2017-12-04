suppressMessages(library(tidyverse))
library(stringr)
suppressMessages(library(magrittr))
suppressMessages(library(glue))
library(limma)
suppressMessages(library(Homo.sapiens))

gene_list_file = commandArgs(trailingOnly = TRUE)
gene_list = read_lines(gene_list_file)

gene_list = gene_list[-(1:which(str_detect(gene_list, '^>')))]
cat(glue('The number of input genes is {length(gene_list)}\n\n'))

# map gene symbols to entrez ids
gene_entrezid = suppressMessages(mapIds(Homo.sapiens,
                       keys = gene_list,
                       keytype = "SYMBOL",
                       column = "ENTREZID"))

# create a clean data_frame with two columns (gene symbol and entrez id)
gene_name2id = enframe(gene_entrezid,
                       name = 'symbol', value = 'entrezid')

if (sum(is.na(gene_name2id$entrezid)) > 0) {
    cat('unmatched genes:\n')
    gene_name2id %>%
        filter(is.na(entrezid)) %>%
        pull(symbol)
}

gene_name2id %<>%
    dplyr::filter(!is.na(entrezid))

cat(glue('The number of matched genes is {nrow(gene_name2id)}\n\n'))

# run GO term enrichment analysis using `goana` function in limma package
goana_res = goana(gene_name2id$entrezid,
                  species = 'Hs')

cat(glue('Write output csv file to {gene_list_file}.goana.csv\n\n'))

# Add detailed gene names for each GO term and write result to local file
# Thresholds:
# * Only use biological pathway
# * Only use GO term with more than 10 genes
# * FDR < 0.01
go_eg_fulltable = toTable(org.Hs.egGO2ALLEGS)
as_data_frame(goana_res) %>%
    dplyr::mutate(GO_id = rownames(goana_res)) %>%
    dplyr::select(GO_id, everything()) %>%
    dplyr::filter(Ont == 'BP', N >= 10) %>%
    dplyr::mutate(FDR = p.adjust(P.DE, method = 'fdr')) %>%
    dplyr::arrange(P.DE) %>%
    dplyr::filter(FDR < 0.01) %>%
    dplyr::mutate(gene_list = map_chr(GO_id, function(id) {
        full_table = filter(go_eg_fulltable, go_id == id)
        gene_name_list = gene_name2id %>%
            dplyr::filter(entrezid %in% full_table$gene_id) %>%
            dplyr::pull(symbol)
        paste(unique(gene_name_list), collapse = ', ')
    })) %>%
    write_csv(paste0(gene_list_file, '.goana.csv'))