# CB vcf flat file
library(cegwas2)
vcf <- data.table::fread("~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/CB4856-20170531.snpeff.vcf.gz", skip = 77)


filtered <- vcf %>%
    dplyr::filter(FILTER == "PASS") %>%
    dplyr::mutate(genotype = stringr::str_split_fixed(CB4856, ":", 2)[,1],
                  annotation = stringr::str_split_fixed(INFO, "ANN=", 2)[,2]) %>%
    dplyr::select(chrom = `#CHROM`, pos = POS, ref = REF, alt = ALT, genotype, annotation) %>%
    dplyr::mutate(annotation = ifelse(annotation == "", NA, annotation))
annotated <- filtered %>%
    tidyr::separate(annotation, into = c("allele", "allele_annotation", "impact", "gene_name", "gene_id", "feature_type", "feature_id", "transcript_biotype",
                                         "rank", "hgvs.c", "hgvs.p", "cDNA.pos/cDNA.length", "CDS.pos/CDS.length", "AA.pos/AA.length", "distance",
                                         "errors", "warnings/info"), sep = "\\|", remove = FALSE)

# annotation: "Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'"

vcf <- "~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/CB4856-20170531.snpeff.vcf.gz"
vcf_header_query <- glue::glue("bcftools view -h {vcf}")
vcf_header <- readr::read_lines(suppressWarnings(pipe(vcf_header_query)))
info_set <- stringr::str_match(vcf_header, '##INFO=<ID=([^,>]+).*Type=([^,>]+).*Description=\"([^,>]+)\"') # nolint
colnames(info_set) <- c("g", "ID", "Type", "Description")
info_columns <- purrr::discard(info_set[, 2], is.na)
info_column_types <- purrr::discard(info_set[, 3], is.na)
"ANN" %in% info_columns
format_set <- stringr::str_match(vcf_header, '##FORMAT=<ID=([^,>]+).*Type=([^,>]+).*Description=\"([^,>]+)\"') # nolint
colnames(format_set) <- c("g", "ID", "Type", "Description")
format_columns <- c(purrr::discard(format_set[, 2], is.na), "TGT")
format_column_types <- c(purrr::discard(format_set[, 3], is.na), "TGT")

elegans_gff <- get_db(table = "wormbase_gene")
query = "V:5000000-5005000"

base_header <- c("CHROM",
                 "POS",
                 "REF",
                 "ALT",
                 "FILTER")
ann_header <- c("allele",
                "effect",
                "impact",
                "gene_name",
                "gene_id",
                "feature_type",
                "feature_id",
                "transcript_biotype",
                "exon_intron_rank",
                "nt_change",
                "aa_change",
                "cdna_position_or_len",
                "protein_position",
                "distance_to_feature",
                "error",
                "extra")
info = c()
format = c("TGT")
info_query <- paste0(info, collapse = "\\t%")
format_query <- paste0(format, collapse = "!%")
query_string <- glue::glue("'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%{info_query}[\\t%{format_query}]\\n'")


command <- paste("bcftools",
                 "query",
                 "--samples CB4856",
                 "-f",
                 glue::glue("'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%{info_query}[\\t%{format_query}]\\n'"),
                 "/Users/katieevans/Dropbox/Andersenlab/Reagents/WormReagents/_SEQ/WI/WI-20170531/vcf/WI.20170531.snpeff.vcf.gz",
                 ">",
                 "/Users/katieevans/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/CB4856-20170531-cegwas2.snpeff.vcf.gz")

debug(query_vcf)
query_vcf("sqst-5", samples = "CB4856", impact = "ALL")

# read in vcf from command line...
result <- try(dplyr::tbl_df(data.table::fread("/Users/katieevans/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/CB4856-20170531-cegwas2.snpeff.vcf.gz",
                                              col.names = c(base_header,
                                                            "INFO",
                                                            "CB4856"),
                                              sep = "\t",
                                              na.strings = c("", "."))),
              silent = FALSE)
tsv <- result %>%
    dplyr::mutate(REF = ifelse(REF == TRUE, "T", REF), # T nucleotides are converted to 'true'
                  ALT = ifelse(ALT == TRUE, "T", ALT))
tsv <- dplyr::mutate(tsv, ANN = "")
tsv <- tsv %>%
    dplyr::mutate(ANN = stringr::str_split(ANN, ",")) %>%
    tidyr::unnest(ANN) %>%
    {
        if (!is.null(ann_header))
            tidyr::separate(., ANN, into = ann_header, sep = "\\|") %>%
            dplyr::select(dplyr::one_of(c(base_header, ann_header)), dplyr::everything(), -extra)
        else
            dplyr::select(., dplyr::one_of(c(base_header)), dplyr::everything())
    } %>%
    dplyr::mutate(query = query, region = region) %>%
    dplyr::select(CHROM, POS, query, region, dplyr::everything())

## here
if ("impact" %in% names(tsv)) {
    tsv <- tsv[tsv$impact %in% impact, ]
} else if (impact != "ALL" | !is.na(impact)) {
    warning("Warning: No ANN column specified; Variants will not be filtered on impact.")
}
tsv <- tidyr::gather_(tsv, "SAMPLE", "FORMAT_COLUMN", samples)  %>%
    tidyr::separate(FORMAT_COLUMN,
                    into = format,
                    sep = "\\!",
                    convert = TRUE,
                    remove = T) %>%
                    {
                        if ("DP" %in% format) dplyr::mutate(., DP = as.integer(ifelse( (DP == ".") | is.na(DP), 0, DP))) else .
                    } %>%
                    {
                        if ("TGT" %in% format)
                            tidyr::separate(.,
                                            TGT,
                                            sep = "\\/|\\|",
                                            into = c("a1", "a2"),
                                            convert = T,
                                            remove = T)
                        else
                            .
                    } %>%
                    {
                        if ("GT" %in% format)
                            tidyr::separate(.,
                                            GT,
                                            sep = "\\/|\\|",
                                            into = c("g1", "g2"),
                                            convert = T,
                                            remove = T) %>%
                            dplyr::mutate_at(as.integer, .vars = c("g1", "g2")) %>%
                            # Why this weird way? It's faster.
                            dplyr::mutate(genotype = as.integer(rowSums(.[, c("g1", "g2")])))
                        else
                            .
                    }
# Fix NAs and convert columns
tsv <- tsv %>%
    dplyr::mutate_all(function(x) { ifelse((x == ".") | (x == ""), NA, x)}) %>%
    # Info columns
    dplyr::mutate_at(.vars = info_columns[info_column_types == "Integer" & info_columns %in% info],
                     as.integer) %>%
    dplyr::mutate_at(.vars = info_columns[info_column_types == "Float" & info_columns %in% info],
                     as.numeric) %>%
    dplyr::mutate_at(.vars = info_columns[info_column_types == "Flag" & info_columns %in% info],
                     as.logical) %>%
    # Format Columns
    dplyr::mutate_at(.vars = format_columns[format_column_types == "Integer" & format_columns %in% format],
                     as.integer) %>%
    dplyr::mutate_at(.vars = format_columns[format_column_types == "Float" & format_columns %in% format],
                     as.numeric) %>%
    dplyr::mutate_at(.vars = format_columns[format_column_types == "Flag" & format_columns %in% format],
                     as.logical)
column_order <- c("CHROM",
                  "POS",
                  "REF",
                  "ALT",
                  "SAMPLE",
                  "FILTER",
                  "FT",
                  "a1",
                  "a2",
                  "g1",
                  "g2",
                  "genotype",
                  "query",
                  "region")
column_order_use <- c(column_order[column_order %in% names(tsv)],
                      names(tsv)[!names(tsv) %in% column_order])

tsv <- tsv %>%
    dplyr::select_at(.vars = column_order_use) %>%
    dplyr::arrange(CHROM, POS)

# Remove ANN field if its empty
if (!("ANN" %in% info_columns)) {
    tsv <- tsv %>% dplyr::select(-ANN)
}


# remake fake crossobject with less strains
pheno <- fst::read_fst("~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/allRIAILsregressed.fst") %>%
    na.omit()
length(unique(pheno$strain))

linkagemapping::load_cross_obj("N2xCB4856cross_full")

View(newcross$pheno)

gene_annotations <- fst::read_fst("~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/gene_annotations.fst")
gene_annotations2 <- gene_annotations %>%
    dplyr::group_by(wbgene) %>%
    dplyr::mutate(go_term = paste(go_term, collapse = "; "),
                  go_name = paste(go_name, collapse = "; "),
                  go_description = paste(go_description, collapse = "; ")) %>%
    dplyr::distinct()
fst::write_fst(gene_annotations2, "~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/gene_annotations2.fst", 100)



#########
# make CB vcf into an sqlite db
library(RSQLite)
vcf <- fst::read_fst("~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/CB4856.smaller.vcf.fst")
db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = "~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/strain_vcf.db")
RSQLite::dbWriteTable(db, "CB4856", vcf)
# RSQLite::dbDisconnect(db)
dbListTables(db)
dbListFields(db, "CB4856")
dbReadTable(db, "CB4856")

res <- RSQLite::dbSendQuery(db, "SELECT * FROM CB4856 WHERE ((CHROM = 'I') AND (POS > 3000))")
RSQLite::dbFetch(res)
