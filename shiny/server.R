#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(magrittr)
library(stringr)
library(limma)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(DT)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

    get_species_dat <- reactive({
        if (input$species == 'human') {
            list(
                eg_db = org.Hs.eg.db,
                go_eg = toTable(org.Hs.egGO2ALLEGS),
                species = 'Hs'
            )
        } else {
            list(
                eg_db = org.Mm.eg.db,
                go_eg = toTable(org.Mm.egGO2ALLEGS),
                species = 'Mm'
            )
        }
    }) %>% debounce(1000)

    get_gene_list <- reactive({
        gene_list = str_split(input$gene, '\\n')[[1]] %>%
            map_chr(str_trim) %>%
            str_subset('.')

        if (length(gene_list) == 0) {
            NULL
        } else {
            gene_list_id = mapIds(get_species_dat()$eg_db,
                                  keys = gene_list,
                                  keytype = 'SYMBOL',
                                  column = 'ENTREZID')
            gene_name2id = enframe(gene_list_id,
                                   name = 'symbol', value = 'entrezid') %>%
                filter(!is.na(entrezid))

            if (nrow(gene_name2id) == 0) {
                NULL
            } else {
                gene_name2id
            }
        }

    }) %>% debounce(1000)

    observeEvent(input$example, {
        updateTextAreaInput(session, 'gene',
                            value = paste(readLines('geneset.txt')[-c(1:2)],
                                          collapse = '\n'))
    })

    observeEvent(input$clear, {
        updateTextAreaInput(session, 'gene', value = '')
    })

    output$goana_res <- renderDataTable({

        gene_list = get_gene_list()
        if (is.null(gene_list)) {
            datatable(data_frame())
        } else {
            withProgress(
                message = 'GO term analysis',
                detail = 'run goana',
                value = 0.1, {
                    goana_res = goana(get_gene_list()$entrezid,
                                      species = get_species_dat()$species)

                    incProgress(detail = "calculate FDR",
                                amount = 0.3)
                    goana_res = as_data_frame(goana_res) %>%
                        dplyr::mutate(GO_id = rownames(goana_res)) %>%
                        dplyr::select(GO_id, everything()) %>%
                        dplyr::filter(Ont == 'BP', N >= 10) %>%
                        dplyr::mutate(FDR = p.adjust(P.DE, method = 'fdr')) %>%
                        dplyr::arrange(P.DE) %>%
                        dplyr::filter(FDR < 0.01) %>%
                        dplyr::slice(1:50)

                    incProgress(detail = "add gene names",
                                amount = 0.3)
                    goana_res %<>%
                        dplyr::mutate(gene_list = map_chr(GO_id, function(id) {
                            full_table = filter(get_species_dat()$go_eg, go_id == id)
                            gene_name_list = get_gene_list() %>%
                                dplyr::filter(entrezid %in% full_table$gene_id) %>%
                                dplyr::pull(symbol)
                            paste(unique(gene_name_list), collapse = ', ')
                        })) %>%
                        dplyr::mutate(P.DE = signif(P.DE, digits = 3),
                                      FDR = signif(P.DE, digits = 3))

                    setProgress(1)
                    goana_res %<>%
                        datatable(
                            extension = "Buttons",
                            options = list(
                                dom = 'Bfrtip',
                                buttons = list(
                                    'copy',
                                    list(
                                        extend = 'collection',
                                        buttons = c('csv', 'excel'),
                                        text = 'Download'
                                    ),
                                    list(
                                        extend = 'colvis',
                                        columns = 3:8
                                    )
                                ),
                                columnDefs = list(list(visible = FALSE,
                                                       targets = c(3:6, 8)))
                            )
                        )
            })

            goana_res
        }
    })

})
