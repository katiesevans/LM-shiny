library(fst)

# read in vcf from command line...
result <- try(dplyr::tbl_df(data.table::fread("/Users/katieevans/Dropbox/AndersenLab/LabFolders/Katie/projects/genome/CB4856-20170531-cegwas2.snpeff.vcf.gz",
                                              col.names = c("CHROM", "POS", "REF", "ALT", "FILTER",
                                                            "INFO",
                                                            "CB4856"),
                                              sep = "\t",
                                              na.strings = c("", "."))),
              silent = FALSE)
tsv <- result %>%
    dplyr::mutate(REF = ifelse(REF == TRUE, "T", REF), # T nucleotides are converted to 'true'
                  ALT = ifelse(ALT == TRUE, "T", ALT))

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

new <- tsv %>%
    dplyr::mutate(SAMPLE = "CB4856") %>%
    tidyr::separate(CB4856, into = c("a1", "a2"), sep = "/") %>%
    tidyr::separate(INFO, into = ann_header, sep = "\\|") %>%
    dplyr::select(dplyr::one_of(c(base_header, ann_header)), dplyr::everything(), -extra) %>%
    dplyr::select(CHROM, POS, REF, ALT, SAMPLE, FILTER, a1, a2, dplyr::everything())

info_columns <- c('INDEL' , 'IDV' , 'IMF' , 'DP' , 'VDB' , 'RPB' , 'MQB' , 'BQB' , 'MQSB' , 'SGB' , 'MQ0F' , 'AD' , 'ICB' , 'HOB' , 'MQ' , 'AC' , 'NS' , 'AN' , 'ANN' , 'LOF' , 'NMD')
info_column_types <- c('Flag' , 'Integer' , 'Float' , 'Integer' , 'Float' , 'Float' , 'Float' , 'Float' , 'Float' , 'Float' , 'Float' , 'Integer' , 'Float' , 'Float' , 'Integer' , 'Integer' , 'Integer' , 'Integer' , 'String' , 'String' , 'String')
format_columns <- c('PL' , 'FT' , 'HP' , 'DP' , 'SP' , 'AD' , 'ADF' , 'ADR' , 'GT' , 'TGT')
format_column_types <- c('Integer' , 'String' , 'String' , 'Integer' , 'Integer' , 'Integer' , 'Integer' , 'Integer' , 'String' , 'TGT')
info <- NULL
format <- 'TGT'

# Fix NAs and convert columns
new2 <- new %>%
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

fst::write_fst(new2, "~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/CB4856.vcf.fst", 100)

test <- new2 %>%
    dplyr::group_by(FILTER) %>%
    dplyr::summarise(num = n())

smaller_vcf <- new2 %>%
    dplyr::select(CHROM, POS, REF, ALT, a1, a2, effect:gene_id)

allRIAILsregressed <- fst::read_fst("~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/allRIAILsregressed.fst")
