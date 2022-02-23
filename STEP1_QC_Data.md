---
title: "Quality control for PDAC single-cell transcriptomic data"
author: "Mohammadali (Sam) Khalilitousi, Yongjin P. Park"
date: "2022-02-23"
font_size: 12pt
output:
    html_document:
        keep_md: yes
---



# Data preprocessing

* FASTQ files for pancreatic ductal adenocarcinoma (PDAC) were obtained from the GSA: `https://ngdc.cncb.ac.cn/gsa/browse/CRA001160`

* Two papers: [Chen _et al._ (2021)](https://pubmed.ncbi.nlm.nih.gov/33819739/) and [Peng _et al._ (2019)](https://pubmed.ncbi.nlm.nih.gov/31273297/)

* FASTQ files were transformed into spliced and unspliced count matrices using kb-python: `https://www.kallistobus.tools/`




```r
result.dir <- "Results/QC"
.mkdir(result.dir)
sample.info.file <- result.dir %&% "/Tab_Samples.csv"

if.needed(sample.info.file,
{
    .samples <-
        readxl::read_xlsx("Data/CRA001160.xlsx", 1) %>%
        select(-`Accession`) %>% 
        left_join(readxl::read_xlsx("Data/CRA001160.xlsx", 3)) %>%
        select(`ID`, `Accession`, `Sample name`, `Public description`, `Sex`)
    fwrite(.samples, sample.info.file)
})
```



## 1. Merge all the files into one file set


```r
result.dir <- "Results/QC"
.mkdir(result.dir)
merged.hdr <- result.dir %&% "/Merged"
raw.data.dirs <- list.files("Data/processed", pattern = "CRR[0-9]+", full.names = TRUE)
```

We merged 36 file sets into one file set, distinguishing the spliced and unspliced features.


```r
merged <- fileset.list(merged.hdr)

if(!all(file.exists(unlist(merged)))){

    .hdrs <- raw.data.dirs %&% "/counts_filtered"

    merged <- rcpp_mmutil_merge_file_sets(
        r_headers = .hdrs,
        r_batches = sapply(raw.data.dirs, basename),
        r_mtx = .hdrs %&% ".mtx.gz",
        r_row = .hdrs %&% ".rows.gz",
        r_col = .hdrs %&% ".cols.gz",
        output = merged.hdr)
}
```

Create a mapping between ENSEMBL gene ID and HGNC symbols.


```r
features <- fread(merged$row, col.names = "feature", header=FALSE)

parse <- function(s){
    transpose(lapply(strsplit(s, "[_.]"), function(x) c(x[1], x[length(x)])))
}

features[, c("ensembl_gene_id", "status") := parse(`feature`)]

feature.info.file <- result.dir %&% "/feature_info.RDS"

if.needed(feature.info.file,
{

    ensembl <- biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL',
                                host='uswest.ensembl.org',
                                path='/biomart/martservice',
                                dataset='hsapiens_gene_ensembl')

    ensembl.hs <- biomaRt::useDataset('hsapiens_gene_ensembl', mart=ensembl)

    .attr <- c('ensembl_gene_id',
               'hgnc_symbol',
               'chromosome_name',
               'transcription_start_site',
               'transcript_start',
               'transcript_end',
               'description',
               'percentage_gene_gc_content')

    .temp <- biomaRt::getBM(attributes=.attr,
                            filters=c('ensembl_gene_id'),
                            values=unique(features$ensembl_gene_id),
                            mart=ensembl.hs,
                            useCache = FALSE)

    .temp <- as.data.table(.temp)

    feature.info <-
        rbind(.temp[transcript_start < transcript_end,
                    .(transcript_start = min(transcript_start),
                      transcript_end = max(transcript_end),
                      percentage_gene_gc_content = mean(percentage_gene_gc_content),
                      chromosome_name = first(chromosome_name),
                      description = first(description)),
                    by = .(ensembl_gene_id, hgnc_symbol)],
              .temp[transcript_start > transcript_end,
                    .(transcript_start = max(transcript_start),
                      transcript_end = min(transcript_end),
                      percentage_gene_gc_content = mean(percentage_gene_gc_content),
                      chromosome_name = first(chromosome_name),
                      description = first(description)),
                    by = .(ensembl_gene_id, hgnc_symbol)])

    saveRDS(feature.info, file=feature.info.file)
})
feature.info <- readRDS(feature.info.file)
```

[feature information](Results/QC/feature_info.RDS)

## 2. Basic Q/C by row-wise and column-wise scores of the combined single-cell RNA-seq data matrix



```r
.mt.genes <-
    left_join(features, feature.info) %>%
    filter(chromosome_name == "MT") %>%
    select(feature) %>%
    unlist
```

```
## Joining, by = "ensembl_gene_id"
```

```r
.mt.only <- fileset.list(result.dir %&% "/MT")

if.needed(.mt.only,
{
    rcpp_mmutil_copy_selected_rows(merged$mtx,
                                   merged$row,
                                   merged$col,
                                   .mt.genes,
                                   result.dir %&% "/MT");
})
```




```r
.mkdir("Table/QC/")
.file <- "Table/QC/cell_scores.csv.gz"
.file.2 <- "Table/QC/feature_scores.csv.gz"

if.needed(c(.file, .file.2),
{

    all.scores <-
        rcpp_mmutil_compute_scores(merged$mtx,
                                   merged$row,
                                   merged$col)

    mt.scores <-
        rcpp_mmutil_compute_scores(.mt.only$mtx,
                                   .mt.only$row,
                                   .mt.only$col)

    score.dt <-
        setDT(all.scores$col) %>%
        left_join(setDT(mt.scores$col),
                  by = "name",
                  suffix = c("", ".MT")) %>%
        as.data.table
    score.dt[is.na(`sum.MT`), `sum.MT` := 0]
    score.dt[, mito := `sum.MT` / `sum` * 100]
    fwrite(score.dt, .file)

    row.score.dt <- setDT(all.scores$row)
    fwrite(row.score.dt, .file.2)
})
score.dt <- fread(.file)
row.score.dt <- fread(.file.2)
```



```r
.hist.nnz <- score.dt[, .(.N), by = .(0.1 * ceiling(10 * log10(nnz)))]
names(.hist.nnz) <- c("nnz.log10", "count")
.hist.mito <- score.dt[, .(.N), by = .(ceiling(mito/2) * 2)]
names(.hist.mito) <- c("mito", "count")

.dt <- score.dt %>%
    mutate(nnz.log10.grid = 0.1 * ceiling(10 * log10(nnz)),
           mito.grid = ceiling(mito))

.hist.joint <- .dt[, .(.N), by = .(nnz.log10.grid, mito.grid)]
```



#### Select high-quality cells that contain at least $>$ 300 non-zero genes expressed with the mitochondrial gene activity $<$ 50 %


![](STEP1_QC_Data_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

[PDF](Figures/QC/Fig_qc_nnz_mito_cutoff.pdf)



#### Of 249,457 cells, we retained 227,331 (91.13 %), discarding 22,126 cells (8.87 %) by applying these Q/C steps


```r
cell.qc.data <- fileset.list(result.dir %&% "/cell_qc")
qc.cells <- score.dt[mito < MITO.CUTOFF & nnz > NNZ.CUTOFF, .(name)] %>%
    unlist

if.needed(cell.qc.data,
{
    rcpp_mmutil_copy_selected_columns(merged$mtx,
                                      merged$row,
                                      merged$col,
                                      qc.cells,
                                      output = result.dir %&% "/cell_qc");
})
```

### Cancer-related genes are highly expressed in the scRNA-seq data


```r
gwas.data.file <- "gwas_catalog_tidy.tsv.gz"
```




```r
.file <- "Table/QC/feature_scores_qc.csv.gz"
if.needed(.file,
{
    all.scores <-
        rcpp_mmutil_compute_scores(cell.qc.data$mtx,
                                   cell.qc.data$row,
                                   cell.qc.data$col)
    .dt <- setDT(all.scores$row)
    fwrite(.dt, .file)
})
feature.scores <- fread(.file) %>%
    rename(`feature` = `name`) %>%
    left_join(features) %>%
    left_join(feature.info) %>%
    as.data.table
```

```
## Joining, by = "feature"
## Joining, by = "ensembl_gene_id"
```



```r
cancer.genes <-
    gwas.data[str_detect(`DISEASE/TRAIT`, "[Cc]ancer"), .(gene_symbol)] %>%
    unique %>%
    unlist

feature.scores[, cancer := "non-CANCER"]
feature.scores[`hgnc_symbol` %in% cancer.genes, cancer := "CANCER"]
feature.scores[str_detect(`description`, "oncogene"), cancer := "CANCER"]               
feature.scores[str_detect(`description`, "tum[ou]+r"), cancer := "CANCER"]
feature.scores[str_detect(`description`, "cancer"), cancer := "CANCER"]
```






![](STEP1_QC_Data_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

[PDF](Figures/QC/Fig_cancer_vs_noncancer.pdf)



```r
us.ratio <-
    feature.scores %>%
    filter(`nnz` >= GENE_NNZ_CUTOFF, `cv` > GENE_CV_CUTOFF) %>% 
    as.data.table %>% 
    dcast(ensembl_gene_id ~ status, value.var = "sum", fun.aggregate = sum, fill = 0) %>%
    filter(spliced > 0, unspliced > 0) %>%
    as.data.table() %>%
    (function(x) {
        .dt <- feature.scores %>%
            mutate(length = abs(transcript_end - transcript_start)) %>% 
            select(ensembl_gene_id, cancer, percentage_gene_gc_content, length)
        left_join(x, unique(.dt))
    }) %>%
    na.omit
```

```
## Joining, by = "ensembl_gene_id"
```





![](STEP1_QC_Data_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

[PDF](Figures/QC/Fig_cancer_vs_noncancer_ratio.pdf)



#### Select high-quality genes that are expressed in at least 500 cells with the coefficient of variation across cells  $>$ 1



```r
qc.data <- fileset.list(result.dir %&% "/final_qc")
qc.features <- features[ensembl_gene_id %in% us.ratio$ensembl_gene_id, .(feature)] %>%
    unique %>% 
    unlist

if.needed(qc.data,
{
    rcpp_mmutil_copy_selected_rows(cell.qc.data$mtx,
                                   cell.qc.data$row,
                                   cell.qc.data$col,
                                   qc.features,
                                   output = result.dir %&% "/final_qc");
})
```



MTX data (22,836 features; 227,331 cells; 329,824,833 non-zero elements; 6.35%):

* [mtx](Results/QC/final_qc.mtx.gz)
* [row](Results/QC/final_qc.rows.gz)
* [col](Results/QC/final_qc.cols.gz)
* [idx](Results/QC/final_qc.mtx.gz.index)
