options(stringsAsFactors = FALSE)

library(tidyverse)
library(data.table)
library(ggrepel)

`%&%` <- function(a, b) paste0(a, b)

sigmoid <- function(x) 1/(1 + exp(-x))

logit <- function(x) log(x) - log(1 - x)

meta.z <- function(z) mean(z) * sqrt(length(z))

num.sci <- function(x) format(x, digits = 2, scientific = TRUE)

num.round <- function(x, d=2) round(x, digits = d) %>% as.character()

num.int <- function(x) format(as.integer(x), big.mark = ',')

.n <- function(...) length(unique(...))

.gsub.rm <- function(x, pat) gsub(x, pattern = pat, replacement = '')

.unlist <- function(...) unlist(..., use.names = FALSE)

.select <- function(...) dplyr::select(...) %>% .unlist()

.mkdir <- function(...) dir.create(..., recursive=TRUE, showWarnings = FALSE)

.list.files <- function(...) list.files(..., full.names = TRUE)

log.msg <- function(...) {
    ss = as.character(date())
    cat(sprintf('[%s] ', ss), sprintf(...), '\n', file = stderr(), sep = '')
}

if.needed <- function(.file, .code) {
    if(!all(file.exists(unlist(.file)))) { .code }
    stopifnot(all(file.exists(unlist(.file))))
}

.gg.save <- function(filename, cat.link=TRUE, ...) {
    if(file.exists(filename)) {
        log.msg('File already exits: %s', filename)
    } else {
        ggplot2::ggsave(..., filename = filename,
                        limitsize = FALSE,
                        units = 'in',
                        dpi = 300,
                        useDingbats = FALSE)
    }
    if(cat.link){
        cat("\n\n[PDF](" %&% filename %&% ")\n\n")
    }
}

.gg.plot <- function(...) {
    ggplot2::ggplot(...) +
        ggplot2::theme_classic() +
        ggplot2::theme(plot.background = element_blank(),
                       plot.margin = unit(c(0,.5,0,.5), 'lines'),
                       panel.background = element_rect(size = 0, fill = 'gray95'),
                       strip.background = element_blank(),
                       legend.background = element_blank(),
                       legend.text = element_text(size = 6),
                       legend.title = element_text(size = 6),
                       axis.title = element_text(size = 8),
                       legend.key.width = unit(.5, 'lines'),
                       legend.key.height = unit(.3, 'lines'),
                       legend.key.size = unit(1, 'lines'),
                       axis.line = element_line(color = 'gray20', size = .5),
                       axis.text = element_text(size = 6))
}

load.data <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}

row.order <- function(mat) {
    require(cba)
    require(proxy)

    if(nrow(mat) < 3) {
        return(1:nrow(mat))
    }

    D = proxy::dist(mat, method <- function(a,b) 1 - cor(a,b, method = 'spearman'))
    D[!is.finite(D)] = 0
    h.out = hclust(D)
    o.out = cba::order.optimal(D, h.out$merge)
    return(o.out$order)
}

col.order <- function(pair.tab, .ro, ret.tab = FALSE) {

    M = pair.tab %>%
        dplyr::select(row, col, weight) %>%
        mutate(row = factor(row, .ro)) %>%
        tidyr::spread(key = col, value = weight, fill = 0)

    co = order(apply(M[, -1], 2, which.max), decreasing = TRUE)
    .co = colnames(M)[-1][co]
    if(ret.tab) {
        ret = pair.tab %>%
            mutate(row = factor(row, .ro)) %>% 
            mutate(col = factor(col, .co))
    } else {
        ret = .co
    }
    return(ret)
}

order.pair <- function(pair.tab, ret.tab=FALSE) {

    require(tidyr)
    require(dplyr)
    
    .tab = pair.tab %>% dplyr::select(row, col, weight)

    M = .tab %>% tidyr::spread(key = col, value = weight, fill = 0)
    rr = M[, 1] %>% unlist(use.names = FALSE)
    cc = colnames(M)[-1] %>% unlist(use.names = FALSE)

    ## log.msg('Built the Mat: %d x %d', nrow(M), ncol(M))
    ro = row.order(M %>% dplyr::select(-row) %>% as.matrix())

    ## log.msg('Sort the rows: %d', length(ro))
    co = order(apply(M[ro, -1], 2, which.max), decreasing = TRUE)

    ## co = row.order(t(M %>% dplyr::select(-row) %>% as.matrix()))
    ## log.msg('Sort the columns: %d', length(co))

    if(ret.tab){
        ret = pair.tab %>%
            mutate(row = factor(row, rr[ro])) %>%
            mutate(col = factor(col, cc[co]))
    } else {
        ret = list(rows = rr[ro], cols = cc[co], M = M)
    }

    return(ret)
}

#' run GOSEQ analysis
#' @param .de.ensg differentially expressed ENSEMBL ID
#' @param .ensg total ENSEMBL ID
run.goseq <- function(.de.ensg, .ensg) {

    .de.ensg <- unique(unlist(.de.ensg))
    .ensg <- unique(unlist(.ensg))

    gene.vector <- as.integer(.ensg %in% .de.ensg)
    names(gene.vector) <- .ensg

    pwf <- goseq::nullp(gene.vector, "hg19", "ensGene", plot.fit = FALSE)

    goseq::goseq(pwf, "hg19", "ensGene", use_genes_without_cat = FALSE) %>%
        as.data.table
}
