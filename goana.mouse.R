#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    cat('Rscript goana.R geneset.txt [output.csv]\n')
    q(save = 'no')
}

gene_list_file = args[1]

if (file.exists(gene_list_file)) {
    suppressMessages(library(tidyverse))
    library(stringr)
    suppressMessages(library(magrittr))
    suppressMessages(library(glue))
    library(progress)
    library(limma)
    suppressMessages(library(org.Mm.eg.db))
    suppressMessages(library(crayon))

    gene_list = read_lines(gene_list_file)
    output_file = paste0(gene_list_file, '.goana.csv')
} else {
    q(save = 'no')
}

if (length(args) >= 2 && dir.exists(dirname((args[2])))) {
    output_file = args[2]
}

gene_list = gene_list[-(1:which(str_detect(gene_list, '^>')))]

# map gene symbols to entrez ids
gene_entrezid = suppressMessages(mapIds(org.Mm.eg.db,
                       keys = gene_list,
                       keytype = "SYMBOL",
                       column = "ENTREZID"))

# create a clean data_frame with two columns (gene symbol and entrez id)
gene_name2id = enframe(gene_entrezid,
                       name = 'symbol', value = 'entrezid')

glue('
    The number of matched genes is \\
    {bold(sum(!is.na(gene_name2id$entrezid)))}/{bold(length(gene_list))}
    ')

if (sum(is.na(gene_name2id$entrezid)) > 0) {
    gene_unmatched = gene_name2id %>%
        filter(is.na(entrezid)) %>%
        pull(symbol)
    if (length(gene_unmatched) > 5) {
        glue('Unmatch genes: {paste(gene_unmatched[1:5], collapse=", ") ...}')
    } else {
        glue('Unmatch genes: {paste(red(gene_unmatched), collapse=", ")}')
    }
}

gene_name2id %<>%
    dplyr::filter(!is.na(entrezid))

# run GO term enrichment analysis using `goana` function in limma package
goana_res = goana(gene_name2id$entrezid,
                  species = 'Mm')

glue('Top gene ontology terms:')
goana_res[, c('Term', 'P.DE')] %>%
    rownames_to_column() %>%
    arrange(P.DE) %>%
    column_to_rownames() %>%
    head()

glue('Write output csv file to {output_file}')

# Add detailed gene names for each GO term and write result to local file
# Thresholds:
# * Only use biological pathway
# * Only use GO term with more than 10 genes
# * FDR < 0.01
go_eg_fulltable = toTable(org.Mm.egGO2ALLEGS)
goana_res %<>%
    rownames_to_column('GO_id') %>%
    as_data_frame() %>%
    dplyr::filter(Ont == 'BP', N >= 10) %>%
    dplyr::mutate(FDR = p.adjust(P.DE, method = 'fdr')) %>%
    dplyr::arrange(P.DE) %>%
    dplyr::filter(FDR < 0.01)

pb <- progress_bar$new(total = nrow(goana_res))
goana_res %>%
    dplyr::mutate(gene_list = map_chr(GO_id, function(id) {
        full_table = filter(go_eg_fulltable, go_id == id)
        gene_name_list = gene_name2id %>%
            dplyr::filter(entrezid %in% full_table$gene_id) %>%
            dplyr::pull(symbol)
        pb$tick()
        paste(unique(gene_name_list), collapse = ', ')
    })) %>%
    write_csv(output_file)
