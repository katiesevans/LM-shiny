
library(tidyverse)
library(linkagemapping)
library(shiny)
library(DT)
library(plotly)
library(fst)

# setwd
# setwd("~/Desktop/LM-shiny/")
# load RIAIL regressed phenotype data to get possible condition/traits
load("data/N2xCB4856cross_full2.Rda")
allRIAILsregressed <- fst::read_fst("data/allRIAILsregressed.fst")
source("query_genes.R")

# Define UI for application
ui <- fluidPage(
   
   # Application title
   titlePanel("Linkage Mapping Analysis - Andersen Lab"),
   
   # no sidebar layout
   shiny::wellPanel(
       # input well
       shiny::fluidRow(
           # conditions
           column(3, shiny::selectInput(inputId = "drug_input", label = "Condition:", choices = sort(unique(allRIAILsregressed$condition)), selected = NULL)),
           # set
           column(2, shiny::radioButtons(inputId = "set_input", label = "RIAIL set:", choices = c(1, 2), selected = 2, inline = TRUE)),
           # trait
           column(3, shiny::selectInput(inputId = "trait_input", label = "Trait:", choices = unique(allRIAILsregressed$trait))),
           # qtl marker
           column(3, shiny::uiOutput("qtl"))
       ),
       # code download button
       shiny::downloadButton(outputId = "code_download", label = "Get code")
   ),
      
  # Main panel for plots
  mainPanel(width = 12,
     # uiOutput("allplots")
     shiny::tabsetPanel(type = "tabs",
                        shiny::tabPanel("QTL Analysis: Condition", shiny::plotOutput("genplot", height = "800px")),
                        shiny::tabPanel("QTL Analysis: Trait",   
                                        # h3("Trait:"),
                                        # shiny::selectInput(inputId = "trait_input", label = NULL, choices = unique(allRIAILsregressed$trait)),
                                        shiny::uiOutput("allplots")),
                        shiny::tabPanel("eQTL Overlap",
                                        # h3(shiny::textOutput("trait_name")), 
                                        # h3("Trait:"),
                                        # shiny::selectInput(inputId = "trait_input2", label = NULL, choices = unique(allRIAILsregressed$trait)),
                                        # shiny::uiOutput("qtl2"),
                                        plotly::plotlyOutput("eqtlplot", height = "500px"),
                                        DT::dataTableOutput("eqtl_data")),
                        shiny::tabPanel("Candidate Genes",
                                        # shiny::uiOutput("qtl"),
                                        shiny::uiOutput("candidate_genes")),
                        shiny::tabPanel("Help"))
  )
)


# Define server logic 
server <- function(input, output) {
    
    # load mapping data per condition
    loaddata <- shiny::reactive({

        cond <- input$drug_input
        
        # load data
        annotatedmap <- fst::read_fst(glue::glue("data/drug_data/{cond}-GWER.chromosomal.annotated.fst"))
        
    })
    
    # plot alllodplots for each condition (or multiple conditions)
    output$genplot <- shiny::renderPlot({
        
        # get drug
        cond <- input$drug_input
        strainset <- input$set_input

        # load data
        annotatedmap <- loaddata()
        
        # plot all traits
        newmap <- annotatedmap %>%
            na.omit() %>%
            dplyr::filter(set == strainset)%>%
            arrange(desc(trait)) %>%
            dplyr::mutate(n2res = ifelse(eff_size < 0, "yes", "no"),
                          ci_l_pos = as.numeric(ci_l_pos),
                          ci_r_pos = as.numeric(ci_r_pos),
                          pos = as.numeric(pos),
                          lod = as.numeric(lod),
                          condition = cond) %>%
            dplyr::mutate(trait = stringr::str_split_fixed(trait, paste0(cond, "."), 2)[,2])
        
        newmap$trait <- factor(newmap$trait, levels = unique(newmap$trait))
            
        
        # get chromosome lengths
        chr_lens <- data.frame(chr = c("I", "II", "III", "IV", "V", "X"), 
                               start = rep(1,6), 
                               end = c(14972282,15173999,13829314,17450860,20914693,17748731),
                               condition = newmap$condition[1],
                               trait = newmap$trait[1])
        
        # plot
        ggplot(newmap)+
            aes(x=pos/1E6, y=trait)+
            theme_bw() +
            viridis::scale_fill_viridis(name = "LOD") + 
            viridis::scale_color_viridis(name = "LOD") +
            geom_segment(aes(x = ci_l_pos/1e6, y = trait, xend = ci_r_pos/1e6, yend = trait, color = lod), size = 2, alpha = 1) +
            geom_segment(data=chr_lens, aes(x =  start/1e6, xend = end/1e6, y = trait, yend=trait), color='transparent', size =0.1) +
            geom_point(aes(fill = lod, shape = n2res), color = "black",size = 3, alpha = 1)+
            scale_shape_manual(values = c("yes" = 24, "no" = 25)) +
            xlab("Genomic position (Mb)") + ylab("") +
            guides(shape = FALSE) +
            theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                  axis.ticks.y = element_blank(),
                  legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10),
                  legend.key.size = unit(.75, "cm"),
                  panel.grid.major.x = element_line(),
                  panel.grid.major.y = element_line(),
                  panel.grid.minor.y = element_blank(),
                  axis.text.y = element_text(size = 10, face = "bold", color = "black"),
                  axis.title.x = element_text(size=12, face="bold", color= "black"),
                  axis.title.y = element_blank(),
                  strip.text.x = element_text(size=12, face="bold", color="black"),
                  strip.text.y = element_text(size=12, face="bold", color="black", angle = 0),
                  strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
                  plot.title = element_text(size=12, face="bold")) +
            facet_grid(condition ~ chr, scales = "free_x", space = "free")
        
        
    })
   
    # plot linkage defaults (LOD, pxg, riail pheno)
    output$allplots <- shiny::renderUI({
        
        # define variables
        trt <- input$trait_input
        cond <- input$drug_input
        strainset <- input$set_input

        # load data
        annotatedmap <- loaddata()

        # filter data
        traitmap <- annotatedmap %>%
            dplyr::filter(trait == glue::glue("{cond}.{trt}"),
                          set == strainset)
        
        pheno <- allRIAILsregressed %>%
            dplyr::filter(condition == cond,
                          trait == trt)
        
        drugcross <- linkagemapping::mergepheno(N2xCB4856cross_full2, pheno, strainset)
        
        #########
        # Plots #
        #########
        
        # to make barplot a shiny, need to split into separate plots...
        output$barplot <- plotly::renderPlotly({
            # riail pheno plot
            riailset <- drugcross$pheno %>%
                dplyr::filter(set == strainset)
            phenodf <- pheno %>%
                dplyr::filter(strain %in% c(riailset$strain, "N2", "CB4856")) %>%
                dplyr::mutate(strain_fill = dplyr::case_when(strain == "N2" ~ "n2",
                                                             strain == "CB4856" ~ "cb",
                                                             TRUE ~ "RIL")) %>%
                dplyr::group_by(strain, trait) %>%
                dplyr::mutate(avg_phen = mean(phenotype, na.rm = T)) %>%
                dplyr::distinct(strain, trait, .keep_all = T) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(norm_pheno = ((avg_phen - min(avg_phen, na.rm = T)) / (max(avg_phen, na.rm = T) - min(avg_phen, na.rm = T)))) %>%
                dplyr::arrange(norm_pheno)
            
            phenodf$strain <- factor(phenodf$strain, levels = unique(phenodf$strain))
            
            # plot
            barplot <- phenodf %>%
                ggplot2::ggplot(.) +
                ggplot2::aes(x = strain, y = norm_pheno, fill = strain_fill, color = strain_fill) +
                ggplot2::geom_bar(stat = "identity") +
                ggplot2::scale_fill_manual(values = c("n2" = "orange", "cb" = "blue", "RIL" = "grey")) +
                ggplot2::scale_color_manual(values = c("n2" = "orange", "cb" = "blue", "RIL" = "grey")) +
                theme_bw(15) +
                theme(axis.ticks.x = element_blank(),
                      axis.text.x = element_blank(),
                      panel.grid = element_blank(),
                      axis.title = element_text(face = "bold", color = "black"),
                      axis.text.y = element_text(face = "bold", color = "black"),
                      legend.position = "none") +
                labs(x = "Strain", y = "Normalized phenotype")
            
            # plotly
            plotly::ggplotly(barplot + aes(text = glue::glue("Strain: {strain}")), tooltip = "text")
        })
        
        # lod plot
        output$lodplot <- shiny::renderPlot({
            lodplot <- linkagemapping::maxlodplot(traitmap) +
                theme(plot.title = element_blank())
            
            if(nrow(traitmap %>% na.omit()) > 0) {
                pxgplot <- linkagemapping::pxgplot(drugcross, traitmap) +
                    theme(plot.title = element_blank(),
                          panel.grid = element_blank())
            } else {
                pxgplot <- ggplot2::ggplot(traitmap) +
                    geom_blank() +
                    geom_text(x = 0.5, y = 0.5, label = "No QTL")
            }
            
            cowplot::plot_grid(lodplot, pxgplot, nrow = 2, align = "v", axis = "lr")
        })
        
        
        # peaks data 
        output$peaks <- DT::renderDataTable({
            traitmap %>%
                na.omit() %>%
                dplyr::arrange(chr, pos) %>%
                dplyr::select(marker, lod, var_exp, eff_size, ci_l_pos, ci_r_pos) %>%
                dplyr::mutate(lod = round(lod, digits = 4),
                              var_exp = round(var_exp, digits = 4),
                              eff_size = round(eff_size, digits = 4))
        })
        
        # what to show
        tagList(
            h3("QTL plots:"),
            plotly::plotlyOutput("barplot", height = "300px"),
            shiny::plotOutput("lodplot", height = "600px"),
            h3("QTL peaks:"),
            DT::dataTableOutput("peaks")
        )
      
    })
    
    get_eQTL <- shiny::reactive({
        
        # first, get QTL
        qtl_marker <- input$whichqtl

        # define variables
        trt <- input$trait_input
        cond <- input$drug_input
        strainset <- input$set_input

        # load data
        annotatedmap <- loaddata()
        
        # filter data
        peaks <- annotatedmap %>%
            dplyr::filter(trait == glue::glue("{cond}.{trt}"),
                          set == strainset) %>%
            na.omit() %>%
            dplyr::filter(marker == qtl_marker)
        
        # made for multiple peaks but I think just one at a time is good.
        all_eQTL <- NULL
        for(i in 1:nrow(peaks)) {
            test <- eqtlmap %>%
                dplyr::filter(chr == peaks$chr[i],
                              ci_l_pos < peaks$ci_r_pos[i],
                              ci_r_pos > peaks$ci_l_pos[i]) %>%
                dplyr::mutate(QTL = peaks$pos[i])
            all_eQTL <- rbind(all_eQTL, test)
        }
        
        # combine to get gene names if available
        newprobes <- all_probe_info %>% 
            dplyr::select(probe:wbgene) %>% 
            dplyr::distinct() %>%
            dplyr::group_by(probe) %>%
            dplyr::mutate(all_genes = paste(unique(gene), collapse = ", "),
                          all_genes2 = paste(unique(wbgene), collapse = ", "),
                          all_genes3 = paste(all_genes, all_genes2, sep = ", ")) %>%
            dplyr::select(probe, gene = all_genes3) %>%
            dplyr::mutate(gene = gsub(", NA", "", gene)) %>%
            dplyr::distinct()
        all_eQTL <- all_eQTL %>%
            dplyr::left_join(newprobes, 
                             by = c("trait" = "probe")) %>%
            # color by distant or local (within 1 Mb of peak)
            dplyr::mutate(class = dplyr::case_when(probe_chr != chr ~ "diff_chr",
                                                   (probe_start + (probe_stop - probe_start)/2) - pos > 1e6 ~ "distant",
                                                   TRUE ~ "cis"))
        
    })
    
    # text output for eQTL trait
    output$trait_name <- renderText({
        input$trait_input
    })
    
    # Show eQTL overlap
    output$eqtlplot <- plotly::renderPlotly({
        
        all_eQTL <- get_eQTL()
        
        # factor to order chr
        all_eQTL$chr <- factor(all_eQTL$chr, levels = c("I", "II", "III", "IV", "V", "X"))
        all_eQTL$probe_chr <- factor(all_eQTL$probe_chr, levels = c("X", "V", "IV", "III", "II", "I"))
        
        # plot eQTL peaks
        tsize <- 12
        eplot <- ggplot(all_eQTL) +
            aes(x = pos / 1e6, y = lod, color = class, size = var_exp) +
            # geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 1, ymax=1.01), color='transparent', fill='transparent', size =0.1, inherit.aes =  F) +
            geom_point() +
            facet_grid(~chr, scales = "free", space = "free") +
            theme_bw(tsize) +
            labs(x = "QTL position (Mb)", y = "LOD") +
            scale_color_manual(values = c("cis" = "grey", "distant" = "yellow", "diff_chr" = "red"), name = "Class") +
            scale_size_continuous(range = c(0.1,2), guide = "none") +
            theme(panel.grid = element_blank(),
                  legend.position = "right",
                  axis.title = element_text(face = "bold"),
                  axis.text.y = element_text(face = "bold", color = "black"),
                  axis.text.x = element_text(face = "bold", color = "black", angle = 90),
                  strip.text = element_text(face = "bold")) +
            scale_alpha(guide = "none") +
            geom_vline(aes(xintercept = QTL/1e6), linetype = "dashed", color = "blue")
        
        plotly::ggplotly(eplot +
                             aes(text = glue::glue("Probe: {trait}\n Gene: {gene}\n Probe_pos: {probe_chr}:{round(probe_start/1e6, digits = 3)} Mb")), 
                             tooltip = "text")

    })
    
    # eQTL dataframe
    output$eqtl_data <- DT::renderDataTable({
        
        all_eQTL <- get_eQTL()
        
        # clean up and print
        all_eQTL %>%
            dplyr::mutate(probe_pos = paste0(probe_chr, "_", probe_start))%>%
            dplyr::select(probe = trait, eQTL_pos = marker, lod, var_exp, eff_size, ci_l_marker, ci_r_marker, probe_pos, class)
        
    })
    
    # get list of QTL to choose from
    output$qtl <- shiny::renderUI({
        
        # get inputs
        trt <- input$trait_input
        cond <- input$drug_input
        strainset <- input$set_input
        # interval <- input$intervals
        
        annotatedmap <- loaddata()
        
        # filter data
        traitmap <- annotatedmap %>%
            dplyr::filter(trait == glue::glue("{cond}.{trt}"),
                          set == strainset) %>%
            na.omit()
        
        # output
        tagList(
            # h3("Select QTL:"),
            shiny::selectInput(inputId = "whichqtl", label = "Select QTL:", choices = unique(traitmap$marker))
        )
    })
    
    # get list of QTL to choose from - for eqtl plot
    # output$qtl2 <- shiny::renderUI({
    #     
    #     # get inputs
    #     trt <- input$trait_input
    #     cond <- input$drug_input
    #     strainset <- input$set_input
    #     # interval <- input$intervals
    #     
    #     annotatedmap <- loaddata()
    #     
    #     # filter data
    #     traitmap <- annotatedmap %>%
    #         dplyr::filter(trait == glue::glue("{cond}.{trt}"),
    #                       set == strainset) %>%
    #         na.omit()
    #     
    #     # output
    #     tagList(
    #         h3("Select QTL:"),
    #         shiny::selectInput(inputId = "whichqtl2", label = NULL, choices = unique(traitmap$marker))
    #     )
    # })
    
    # candidate gene function dataframe output
    output$candidate_genes <- shiny::renderUI({
        
        # first, get QTL
        qtl_marker <- input$whichqtl
        
        # get inputs
        trt <- input$trait_input
        cond <- input$drug_input
        strainset <- input$set_input
        # interval <- input$intervals
        
        annotatedmap <- loaddata()
        
        # filter data
        traitmap <- annotatedmap %>%
            dplyr::filter(trait == glue::glue("{cond}.{trt}"),
                          set == strainset) %>%
            na.omit()
        
        # get the region from the marker
        markerdf <- traitmap %>%
            dplyr::filter(marker == qtl_marker) %>%
            dplyr::mutate(reg = paste0(chr, ":", ci_l_pos, "-", ci_r_pos))
        
        # next, run query genes
        # first element: dataframe,second element: text
        test <- query_genes(markerdf$reg)
        
        # organize data and add wormbase links
        df <- test[[1]] %>%
            dplyr::group_by(wbgene) %>%
            dplyr::mutate(go_term = paste(go_term, collapse = "; "),
                          go_name = paste(go_name, collapse = "; "),
                          go_description = paste(go_description, collapse = "; ")) %>%
            dplyr::distinct() %>%
            dplyr::select(-go_term, -go_description, -gene_class_description, -gene_id, -go_annotation) %>%
            dplyr::rename(GO_term = go_name) %>%
            dplyr::ungroup() %>%
            # dplyr::mutate(WormBase = paste0("https://wormbase.org/species/c_elegans/gene/", wbgene)) %>%
            dplyr::mutate(wbgene = paste0("<a href='",paste0("https://wormbase.org/species/c_elegans/gene/", wbgene),"'>",wbgene,"</a>"))

        output$dataframe <- DT::renderDataTable(escape = FALSE, {
            df
        })
        
        # might need to structure this text
        output$gene_text <- shiny::renderUI({

            # print as bullets
            tags$div(
                tags$ul(
                    tags$li(test[[2]][1]),
                    tags$li(test[[2]][2]),
                    tags$li(test[[2]][3]),
                    tags$li(test[[2]][4]),
                    tags$li(test[[2]][5]),
                    tags$li(test[[2]][6]),
                    tags$li(test[[2]][7])
                )
            )
        })
        
        # what to print to tab
        tagList(
            shiny::uiOutput("gene_text"),
            shiny::h3("List of genes:"),
            DT::dataTableOutput("dataframe")
        )

    })
    
    # download report
    output$code_download <- shiny::downloadHandler(
        
        # Dynamic file name -- WHY IS THIS ONLY PARTIALLY WORKING?
        filename = function() {
            # load inputs
            trt <- input$trait_input
            cond <- input$drug_input
            strainset <- as.numeric(input$set_input)
            glue::glue("{gsub('-', '', Sys.Date())}_{input$drug_input}_{input$trait_input}_set{input$set_input}_report.html")
            },
        
        content = function(file) {
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            tempReport <- file.path(tempdir(), "report.Rmd")
            file.copy("report.Rmd", tempReport, overwrite = TRUE)

            # load inputs
            trt <- input$trait_input
            cond <- input$drug_input
            strainset <- as.numeric(input$set_input)
            qtl_marker <- input$whichqtl
            
            # load data
            annotatedmap <- loaddata() %>%
                dplyr::filter(set == strainset)
            load("data/N2xCB4856cross_full2.Rda")
            allRIAILsregressed <- fst::read_fst("data/allRIAILsregressed.fst")
            eqtlmap <- data.table::fread("data/eqtlmap.tsv")
            all_probe_info <- data.table::fread("data/all_probe_info.tsv.gz")
            gene_annotations <- data.table::fread("data/gene_annotations.tsv.gz")
            
            # Knit the document - use local environment to keep all the ^ above variables ^
            rmarkdown::render(tempReport, output_file = file)
            
        }
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)

