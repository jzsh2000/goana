library(tidyverse)
library(limma)
library(Homo.sapiens)

# Hallmark gene set (apoptosis)
# Genes mediating programmed cell death (apoptosis) by activation of caspases.
# See: <http://software.broadinstitute.org/gsea/msigdb>
HALLMARK_APOPTOSIS <- c(
    'ADD1', 'AIFM3', 'ANKH', 'ANXA1', 'APP', 'ATF3', 'AVPR1A', 'BAX', 'BCAP31',
    'BCL10', 'BCL2L1', 'BCL2L10', 'BCL2L11', 'BCL2L2', 'BGN', 'BID', 'BIK',
    'BIRC3', 'BMF', 'BMP2', 'BNIP3L', 'BRCA1', 'BTG2', 'BTG3', 'CASP1', 'CASP2',
    'CASP3', 'CASP4', 'CASP6', 'CASP7', 'CASP8', 'CASP9', 'CAV1', 'CCNA1', 'CCND1',
    'CCND2', 'CD14', 'CD2', 'CD38', 'CD44', 'CD69', 'CDC25B', 'CDK2', 'CDKN1A',
    'CDKN1B', 'CFLAR', 'CLU', 'CREBBP', 'CTH', 'CTNNB1', 'CYLD', 'DAP', 'DAP3',
    'DCN', 'DDIT3', 'DFFA', 'DIABLO', 'DNAJA1', 'DNAJC3', 'DNM1L', 'DPYD', 'EBP',
    'EGR3', 'EMP1', 'ENO2', 'ERBB2', 'ERBB3', 'EREG', 'ETF1', 'F2', 'F2R', 'FAS',
    'FASLG', 'FDXR', 'FEZ1', 'GADD45A', 'GADD45B', 'GCH1', 'GNA15', 'GPX1', 'GPX3',
    'GPX4', 'GSN', 'GSR', 'GSTM1', 'GUCY2D', 'H1F0', 'HGF', 'HMGB2', 'HMOX1',
    'HSPB1', 'IER3', 'IFITM3', 'IFNB1', 'IFNGR1', 'IGF2R', 'IGFBP6', 'IL18',
    'IL1A', 'IL1B', 'IL6', 'IRF1', 'ISG20', 'JUN', 'KRT18', 'LEF1', 'LGALS3',
    'LMNA', 'LPPR4', 'LUM', 'MADD', 'MCL1', 'MGMT', 'MMP2', 'NEDD9', 'NEFH',
    'PAK1', 'PDCD4', 'PDGFRB', 'PEA15', 'PLAT', 'PLCB2', 'PMAIP1', 'PPP2R5B',
    'PPP3R1', 'PPT1', 'PRF1', 'PSEN1', 'PSEN2', 'PTK2', 'RARA', 'RELA', 'RETSAT',
    'RHOB', 'RHOT2', 'RNASEL', 'ROCK1', 'SAT1', 'SATB1', 'SC5DL', 'SLC20A1',
    'SMAD7', 'SOD1', 'SOD2', 'SPTAN1', 'SQSTM1', 'TAP1', 'TGFB2', 'TGFBR3',
    'TIMP1', 'TIMP2', 'TIMP3', 'TNF', 'TNFRSF12A', 'TNFSF10', 'TOP2A', 'TSPO',
    'TXNIP', 'VDAC2', 'WEE1', 'XIAP'
)

# map gene symbols to entrez ids
gene_entrezid = mapIds(Homo.sapiens,
                       keys = HALLMARK_APOPTOSIS,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

# create a clean data_frame with two columns (gene symbol and entrez id)
gene_name2id = enframe(gene_entrezid,
                       name = 'symbol', value = 'entrezid') %>%
    dplyr::filter(!is.na(entrezid))

# run GO term enrichment analysis using `goana` function in limma package
goana_res = goana(gene_name2id$entrezid,
                  species = 'Hs')

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
    write_csv('goana_res.csv')
