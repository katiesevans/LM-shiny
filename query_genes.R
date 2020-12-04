library(linkagemapping)

# load gene descriptions
gene_annotations <- fst::read_fst("data/gene_annotations.fst")
db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "data/strain_vcf.db")

# db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "https://github.com/katiesevans/LM-shiny/blob/main/data/strain_vcf.db?raw=true")


# vcf <- fst::read_fst("data/CB4856.smaller.vcf.fst")
# data("eQTLpeaks")
# data("probe_info")

# Look for genes in interval
# new update 20201013 - simpler (better??) eQTL predictions
query_genes <- function(region, GO = NULL, strain = "CB4856") {
    
    # Region
    chrom <- stringr::str_split_fixed(region, ":", 2)[,1]
    left_pos <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,1])
    right_pos <- as.numeric(stringr::str_split_fixed(stringr::str_split_fixed(region, ":", 2)[,2], "-", 2)[,2])
    
    # how many genes are in the interval?
    # all_genes <- cegwas2::query_vcf(region, impact = "ALL", samples = strain)
    query <- RSQLite::dbSendQuery(db, glue::glue("SELECT * FROM {strain} WHERE ((CHROM = '{chrom}') AND (POS >= {left_pos}) AND (POS <={right_pos}))"))
    all_genes <- RSQLite::dbFetch(query) 
    t1 <- (glue::glue("There are {length(unique(all_genes$gene_id))} genes in the interval {region}"))
    
    # filter eqtl to > 1% VE
    eQTLpeaks2 <- eQTLpeaks %>%
        dplyr::filter(var_exp >= 0.01)
    
    # how many eQTL map to this region?
    all_eQTL <- eQTLpeaks2 %>%
        dplyr::filter(chr == chrom,
                      ci_l_pos < right_pos,
                      ci_r_pos > left_pos)
    t2 <- (glue::glue("There are {nrow(all_eQTL)} eQTL ({length(unique(all_eQTL$trait))} traits) with VE > 5% that map to {region}"))
    
    # all eQTL probes
    all_eQTL_probes <- probe_info %>%
        dplyr::filter(probe%in% all_eQTL$trait)
    
    ##############################
    # if wbgene is NA - try to fix
    ##############################
    
    # filter na
    # filter na
    na1 <- all_eQTL_probes %>%
        dplyr::group_by(probe) %>%
        dplyr::mutate(num_na = sum(is.na(wbgene))/length(wbgene)) %>%
        dplyr::filter(num_na == 1)
    
    unique_probes <- paste(unique(na1$probe), collapse = ",")
    t3 <- (glue::glue("There are {nrow(na1)} genes with an eQTL that need to be hand curated: {unique_probes}"))
    
    ##################################
    
    # which of the eQTL are overlapping with genes in interval?
    eQTL_outside_CI <- all_eQTL_probes %>%
        dplyr::filter(!wbgene %in% all_genes$gene_id)
    t4 <- (glue::glue("There are {nrow(all_eQTL)-length(unique(eQTL_outside_CI$wbgene))-nrow(na1)} genes in the region with an eQTL, {length(unique(eQTL_outside_CI$wbgene))} genes outside the region with an eQTL, and {nrow(na1)} unknown"))
    
    # Total genes of interest:
    t5 <- (glue::glue("There are at least {length(unique(all_genes$gene_id)) + length(unique(eQTL_outside_CI$wbgene))} total genes of interest."))
    
    # how many of the genes in interval have variation?
    vars <- all_genes %>%
        dplyr::mutate(GT = ifelse(a1 == REF, "ref", "alt")) %>%
        dplyr::filter(GT == "alt")
    
    # genes with protein coding vars
    proteincode <- vars %>%
        dplyr::filter(impact %in% c("MODERATE", "HIGH"))
    t6 <- (glue::glue("There are {length(unique(vars$gene_id))}/{length(unique(all_genes$gene_id))} genes in interval with genetic variation, {length(unique(proteincode$gene_id))}/{length(unique(vars$gene_id))} have protein-coding variation"))
        
    # return final dataframe with all info (might be off, only has 133 instead of 134?)
    total_genes <- gene_annotations %>%
        dplyr::filter(wbgene %in% c(all_genes$gene_id, eQTL_outside_CI$wbgene)) %>%
        dplyr::mutate(inside_CI = ifelse(wbgene %in% all_genes$gene_id, T, F),
                      eqtl = ifelse(wbgene %in% all_eQTL_probes$wbgene, T, F),
                      vars = ifelse(wbgene %in% vars$gene_id, T, F),
                      pc_vars = ifelse(wbgene %in% proteincode$gene_id, T, F),
                      go_annotation = NA)
    
    distinct <- total_genes %>%
        dplyr::distinct(wbgene, .keep_all = T)
    
    # pc alone
    pc <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == F))
    
    # pc + eQTL
    pc_eqtl <- nrow(distinct %>% dplyr::filter(pc_vars == T, eqtl == T))
    
    # eQTL alone
    e <- nrow(distinct %>% dplyr::filter(pc_vars == F, eqtl == T))
    
    t7 <- (glue::glue("There are {pc + pc_eqtl + e} genes with protein-coding variation and/or an eQTL (top priority)"))
    
    text_element <- c(t1, t2, t3, t4, t5, t6, t7)
    
    return(list(total_genes, text_element))
}    
