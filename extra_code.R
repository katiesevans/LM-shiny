setwd("~/Desktop/LM-shiny/")
# load RIAIL regressed phenotype data to get possible condition/traits
# load("data/N2xCB4856cross_full2.Rda")
allRIAILsregressed <- data.table::fread(glue::glue("data/allRIAILsregressed.tsv.gz"))


load(glue::glue("data/set{strainset}/{cond}-GWER.chromosomal.annotated.Rda"))
readr::write_tsv(annotatedmap, glue::glue("data/set{strainset}/{cond}-GWER.chromosomal.annotated.tsv"))

system.time(load(glue::glue("data/set{strainset}/{cond}-GWER.chromosomal.annotated.Rda"))) #5.867s
system.time(readr::read_tsv(glue::glue("data/set{strainset}/{cond}-GWER.chromosomal.annotated.tsv"))) #7.338s
system.time(data.table::fread(glue::glue("data/set{strainset}/{cond}-GWER.chromosomal.annotated.tsv"))) #1.65s


test <- data.table::fread(glue::glue("data/set{strainset}/{cond}-GWER.chromosomal.annotated.tsv"))

# change flx-250 to flx.250 etc. for all set1
# get all conditions with -:
dash <- grep("-", unique(test$condition), value = T)
# change to "." to load data
dot <- gsub("-", "\\.", grep("-", unique(test$condition), value = T))

for(d in grep("-", unique(test$condition), value = T)) {
    # load data
    annotatedmap <- data.table::fread(glue::glue("data/set1/{d}-GWER.proximal.annotated.tsv"))
    
    d2 <- gsub("-", "\\.", d)
    
    annotatedmap$trait <- gsub(d, d2, annotatedmap$trait)
    
    # save as tsv
    readr::write_tsv(annotatedmap, glue::glue("data/set1/{d2}-GWER.proximal.annotated.tsv"))
}

allRIAILsregressed <- data.table::fread(glue::glue("data/allRIAILsregressed.tsv")) %>%
    dplyr::select(round:phenotype)
allRIAILsregressed$condition <- gsub("-", "\\.", allRIAILsregressed$condition)

readr::write_tsv(allRIAILsregressed, "data/allRIAILsregressed.tsv")

load("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/RIAIL_data/allRIAILsPCregressed.Rda")



# set 2: chromosomal to proximal
for(d in unique(allRIAILsregressed$condition)) {
    # load data
    map <- data.table::fread(glue::glue("data/set2/{d}-GWER.chromosomal.annotated.tsv")) %>%
        dplyr::filter(!grepl("PC", trait))
    
    readr::write_tsv(map, glue::glue("data/set2/{d}-GWER.chromosomal.annotated.tsv"))
    
    # re-annotate
    map2 <- map %>%
        dplyr::select(marker:iteration)
    
    # phenotype
    drugcross <- linkagemapping::mergepheno(N2xCB4856cross_full2, allRIAILsregressed %>% 
                                                dplyr::filter(condition == d))
    
    annotatedmap <- linkagemapping::annotate_lods(map2, drugcross, cutoff = "proximal")
    
    # save as tsv
    readr::write_tsv(annotatedmap, glue::glue("data/set2/{d}-GWER.proximal.annotated.tsv"))
}


# set 1: proximal to chromosomal
for(d in unique(allRIAILsregressed$condition)) {
    # load data
    map <- data.table::fread(glue::glue("data/set1/{d}-GWER.proximal.annotated.tsv")) %>%
        dplyr::filter(!grepl("PC", trait))
    
    readr::write_tsv(map, glue::glue("data/set1/{d}-GWER.proximal.annotated.tsv"))
    
    # re-annotate
    map2 <- map %>%
        dplyr::select(marker:iteration)
    
    # phenotype
    drugcross <- linkagemapping::mergepheno(N2xCB4856cross_full2, allRIAILsregressed %>% 
                                                dplyr::filter(condition == d))
    
    annotatedmap <- linkagemapping::annotate_lods(map2, drugcross, cutoff = "chromosomal")
    
    # save as tsv
    readr::write_tsv(annotatedmap, glue::glue("data/set1/{d}-GWER.chromosomal.annotated.tsv"))
}

# combine chrom and proximal and set1 set2 to one dataframe for quicker times?
setwd("~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/")
d <- "HT115"
# for(d in unique(allRIAILsregressed$condition)) {
    df1 <- data.table::fread(glue::glue("data/set1/{d}-GWER.proximal.annotated.tsv")) %>%
        dplyr::mutate(set = 1,
                      annotation = "proximal")
    df2 <- data.table::fread(glue::glue("data/set2/{d}-GWER.proximal.annotated.tsv"))  %>%
        dplyr::mutate(set = 2,
                      annotation = "proximal")
    df3 <- data.table::fread(glue::glue("data/set1/{d}-GWER.chromosomal.annotated.tsv"))  %>%
        dplyr::mutate(set = 1,
                      annotation = "chromosomal")
    df4 <- data.table::fread(glue::glue("data/set2/{d}-GWER.chromosomal.annotated.tsv"))  %>%
        dplyr::mutate(set = 2,
                      annotation = "chromosomal")
    alldata <- df1 %>%
        dplyr::bind_rows(df2, df3, df4)
    readr::write_tsv(alldata, glue::glue("data/drug_data/{d}-GWER.annotated.tsv"))
# }

# why didn'tit work for zinc?20201013 - ran out of space on computer i think lol
test2 <- data.table::fread("data/drug_data/daunorubicin-GWER.annotated.tsv.gz")
saveRDS(test, "data/drug_data/daunorubicin-GWER.annotated.rds")
unique(test$annotation)
unique(test$set)

eqtl_probes2 <- eqtl_probes %>%
    dplyr::mutate(gene = paste0(gene_id, ",", gene_name)) %>%
    dplyr::select(probe, gene) %>%
    tidyr::separate_rows(gene, sep = ",") %>%
    dplyr::mutate(gene = ifelse(gene == "", NA, gene)) %>%
    na.omit() %>%
    dplyr::distinct()

ann2 <- gene_annotations %>%
    dplyr::mutate(gene = paste0(gene_name, ",", gene_id)) %>%
    dplyr::select(-gene_name, -gene_id) %>%
    tidyr::separate_rows(gene, sep = ",") %>%
    dplyr::distinct()

all_probe_info <- eqtl_probes2 %>%
    dplyr::left_join(ann2)
save(all_probe_info, file = "~/Dropbox/AndersenLab/LabFolders/Katie/scripts_kse/all_probe_info.Rda")
readr::write_tsv(all_probe_info, "~/Dropbox/AndersenLab/LabFolders/Katie/git/LM-shiny/data/all_probe_info.tsv")

test <- all_eQTL %>%
    dplyr::left_join(all_probe_info %>% dplyr::select(probe:wbgene) %>% dplyr::distinct(), by = c("trait" = "probe"))

test2 <- all_probe_info %>% 
    dplyr::select(probe:wbgene) %>% 
    dplyr::distinct() %>%
    dplyr::group_by(probe) %>%
    dplyr::mutate(all_genes = paste(unique(gene), collapse = ", "),
                  all_genes2 = paste(unique(wbgene), collapse = ", "),
                  all_genes3 = paste(all_genes, all_genes2, sep = ", ")) %>%
    dplyr::select(probe, gene = all_genes3) %>%
    dplyr::mutate(gene = gsub(", NA", "", gene))

paste(unique(allRIAILsregressed$condition), collapse = " ")

# gzip - less data size
# for i in carmustine chlorothanil docetaxel etoposide fluoxetine.125 fluoxetine.250 irinotecan methotrexate.625 methotrexate.3125 thiabendazole.625 thiabendazole.125 albendazole amsacrine bortezomib chlorpyrifos dactinomycin fenbendazole.15 fenbendazole.30 mebendazole topotecan mianserin monepantel arsenicdibasic arsenictrioxide bleomycin cadmium copper deiquat FUdR mechlorethamine nickel paraquat puromycin silver vincristine zinc cisplatin.250 cisplatin.500 lysate.175 OP50 DA837 JUb68 HT115
# do
# gzip $i-GWER.annotated.tsv
# done


# how do I pair down the data significantly?
# maybe I have a peak df with extra info but only peaks and another df with only lods for all markers?
# It would help if I didn't need all the traits...

df4_peaks <- df4 %>%
    na.omit()

df4_lods <- df4 %>%
    dplyr::select(chr, pos, trait, lod)

# 20201016 - keep only 30 traits + chromosomal + GWER for both set 1 and set 2 (only get rid of flour traits, not cv/iqr)
allpeaks <- NULL
for(d in unique(allRIAILsregressed$condition)) {
    
    # read in set1 data first
    set1 <- data.table::fread(glue::glue("data/set1/{d}-GWER.chromosomal.annotated.tsv"))  %>%
        dplyr::filter(!grepl("red|green|yellow", trait)) %>%
        dplyr::select(-chromosomal) %>%
        dplyr::rename(set = `1`)
    
    # read in set2 data first
    set2 <- data.table::fread(glue::glue("data/set2/set2-{d}-GWER.chromosomal.annotated.tsv")) %>%
        dplyr::filter(!grepl("red|green|yellow", trait)) %>%
        dplyr::select(-chromosomal) %>%
        dplyr::rename(set = `2`)
    
    # bind data
    allsets <- rbind(set1, set2)
    # data.table::fwrite(allsets, "data/drug_data/{d}-GWER.chromosomal.annotated.tsv")
    fst::write_fst(allsets, glue::glue("data/drug_data/{d}-GWER.chromosomal.annotated.fst"), 100)
    # test2 <- fst::read_fst(glue::glue("data/drug_data/{d}-GWER.chromosomal.annotated.fst"))
    
    # keep only peaks
    peaks <- allsets %>%
        na.omit()
    
    allpeaks <- rbind(allpeaks, peaks)
    
}
fst::write_fst(allpeaks, glue::glue("data/drug_data/alldrugpeaks.fst"), 100)

# also filter allRIAILsregressed
allRIAILsregressed <- allRIAILsregressed %>%
    dplyr::filter(!grepl("red|green|yellow", trait)) %>%
    dplyr::select(assay, condition, control, strain, trait, phenotype)
fst::write_fst(allRIAILsregressed, "data/allRIAILsregressed.fst", 100)
