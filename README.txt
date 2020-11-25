# LM-shiny

This R shiny web application was created (2020) to easily browse linkage mapping analyses from the Andersen Lab drug conditions (phenotypes collected 2014).

***Link to Shiny app: [here]***

## Explanation of Functionality
To begin, user can select one of the provided conditions and RIAIL set (see below for details). A summary of each tab can be found as follows:

1. The initial tab (*QTL Analysis: Condition*) will plot all QTL for all traits for the selected condition and set. 
2. The second tab (*QTL Analysis: Trait*) requires the user to select a trait of interest from the drop down menu on top (`cv.EXT` is the default). This tab shows three plots follwed by a searchable dataframe for all trait-specific QTL
  - The first plot is a `plotly` interactive bar graph showing the normalized phenotype of all RIAILs in the experiment with N2 (orange) and CB4856 (blue) highlighted. Hovering your mouse over a particular bar will display which strain it is.
  - The second plot is the result of the linkage mapping analysis (LOD plot). Significant QTL are represented by a red triangle and the blue rectangles are the 95% confidence intervals
  - The third plot shows the phenotype-by-genotype splits of the RIAILs for each significant QTL (or displays "No QTL" if applicable). Strains with the N2 allele at the QTL are colored orange and strains with the CB4856 allele at the QTL are colored blue.
3. The third tab (*eQTL Overlap*) requires the user to select one significant QTL from the condition-trait using the drop down menu along the top (not applicable if no significant QTL). The resulting output is another `plotly` interactive graph showing the positions of all significant eQTL that overlaps with the drug-response QTL interval. The dotted blue line represents the peak drug-response QTL marker. eQTL are colored based on the distance between the eQTL and the gene whose expression is being measured (grey/cis < 1 Mb; yellow/distant > 1 Mb; red/diff_chr = different chromosome). eQTL point size also increases corresponding to the percent of phenotypic variance of the RIAILs explained by this locus. Hovering your mouse over an eQTL will display the probe ID, gene ID (if applicable), and genomic position of the probe.
4. The final tab (*Candidate Genes*) also requires the user to select a significant QTL. The application uses WS273 to analyze the number of genes and eQTL inside the QTL confidence interval and uses CeNDR (v.XXX) to identify genetic variants in the CB4856 strain that can help prioritize candidate genes. A datatable of all genes of interest is found below with an html link to WormBase for each gene.

Figures and code used to generate all figures (along with absolute file paths to the original data on the Andersen Lab Dropbox) can be downloaded using the "**Get Code**" button along the top below the input values.

## Available data:
- **46 conditions**: (carmustine, chlorothanil, daunorubicin, docetaxel, etoposide, fluoxetine.125, fluoxetine.250, irinotecan, methotrexate.625, methotrexate.3125, thiabendazole.625, thiabendazole.125, tunicamycin, abamectin, albendazole, amsacrine, bortezomib, chlorpyrifos, dactinomycin, fenbendazole.15, fenbendazole.30, mebendazole, topotecan, mianserin, monepantel, arsenicdibasic, arsenictrioxide, bleomycin, cadmium, copper, deiquat, FUdR, mechlorethamine, nickel, paraquat, puromycin, silver, vincristine, zinc, cisplatin.250, cisplatin.500, lysate.175, OP50, DA837, JUb68, HT115)
- **30 traits**: (cv.EXT, cv.TOF, f.ad, f.L1, f.L2L3, f.L4, iqr.EXT, iqr.TOF, mean.EXT, mean.norm.EXT, mean.TOF, median.EXT, median.norm.EXT, median.TOF, n, norm.n, q10.EXT, q10.norm.EXT, q10.TOF, q25.EXT, q25.norm.EXT, q25.TOF, q75.EXT, q75.norm.EXT, q75.TOF, q90.EXT, q90.norm.EXT, q90.TOF, var.EXT, var.TOF)
  - *Fluorescent traits removed for lack of use and in order to save space*
- **RIAIL sets**: set1 (Rockman and Kruglyak 2009) & set2 (Andersen *et al.* 2015)
<img src="riail_info.png"/>

<style type="text/css">
           body {          
           max-width:100%;
           padding:0;
           }
</style>
