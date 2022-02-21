require(torch)
require(mmutilR)
require(data.table)
require(tidyverse)

VelocityData <- dataset(
    name = "single-cell velocity data",
    initialize = function(.hdr, MY.DEV = torch_device("cpu")) {
        self$sc.data <- sc.data <- self$.add.ext(.hdr)
        self$DEV <- MY.DEV
        self$info <- rcpp_mmutil_info(self$sc.data$mtx)
        self$genes <-
            fread(sc.data$row, header=FALSE, col.names="gene") %>%
            (function(x){
                x[, loc := 1:nrow(x)]
                x[, c("hgnc", "tag") := tstrsplit(gene, split="[_]")]
                return(x)
            }) %>%
            dcast(hgnc ~ tag, value.var = "loc",
                  fun.aggregate = function(x) x[1]) %>%
            na.omit
        self$genes[, i:= 1:nrow(self$genes)]
        self$mem.idx <- rcpp_mmutil_read_index(sc.data$idx)
        self$read.full()
    },
    read.full = function() {
        .mtx <- self$sc.data$mtx
        .idx <- self$mem.idx
        .loc <- 1:self$info$max.col
        .out <- rcpp_mmutil_read_columns_sparse(.mtx, .idx, .loc)
        .ind <- rbind(.out$col, .out$row)    # column is row & vice versa
        .sz <- c(.out$max.col, .out$max.row) # size
        log.msg("Read the full sparse matrix")
        xx <- torch_sparse_coo_tensor(.ind, .out$val, .sz,
                                      dtype = torch_float16())
        xx <- xx$to_dense()
        log.msg("Porting to the device")
        self$x.list <- lapply(1:nrow(xx), function(r) {
            cat(r,"\r",file=stderr());flush(stderr())
            xx[r, , drop=FALSE]$to_sparse()$to(dtype = torch_float16(),
                                               device=self$DEV)
        })
        log.msg("Distributed to the device")
    },
    .getitem = function(.loc) {
        xx <- torch_cat(self$x.list[.loc], dim=1)$to_dense()$to(dtype=torch_float(), device=self$DEV)
        list(spliced = xx[, self$genes$s], unspliced = xx[, self$genes$u])
    },
    .length = function() {
        self$info$max.col
    },
    dim = function() {
        nrow(self$genes)
    },
    .dim = function() {
        self$dim()
    },
    .add.ext = function(.hdr){
        .names <- c("mtx","row","col","idx")
        str_c(.hdr, c(".mtx.gz", ".rows.gz", ".cols.gz", ".mtx.gz.index")) %>%
            as.list() %>%
            (function(x) { names(x) <- .names; x })
    },
    .files.exist = function(...){
        lapply(..., file.exists) %>%
            unlist() %>%
            all()
    })

#' @param .mean mean vector
#' @param .lnvar log-variance (independent)
normal.stoch <- function(.mean, .lnvar) {
    .eps <- torch_randn_like(.lnvar)
    .sig <- torch_exp(.lnvar / 2.)
    .mean + .eps * .sig
}

kl.loss <- function(.mean, .lnvar) {
    -0.5 * torch_sum(1. + .lnvar - torch_pow(.mean, 2.) - torch_exp(.lnvar), dim = -1);
}

#' FC ReLU layer generator
#' @param in_ input dimension
#' @param out_ output dimension
#' @return a generator
#'
build.fc.relu <- function(in_, out_) {
    nn_sequential(nn_linear(in_, out_),
                  nn_relu())
}

#' Build a stack of layers
#' @param in_ input dimension (first layer)
#' @param layers dimensions for the subsequent layers
#' @param generator a shared generator function
#' @param name a header for the names
#'
build.stack <- function(in_, layers, generator = nn_linear, name = "", stdizer = NULL) {
    ret <-
        nn_module(
            classname = "stack.layers",
            initialize = function(in_,
                                  .layers,
                                  .generator,
                                  .name) {
                d.prev <- in_
                for(l in 1:length(.layers)) {
                    d.this <- .layers[l]
                    .name.l <- paste0(name, ".a", (l - 1))
                    self$add_module(.name.l, module = .generator(d.prev, d.this))
                    if(!is.null(stdizer)){
                        .name.l <- paste0(name, ".d", (l - 1))
                        self$add_module(.name.l, module = stdizer(d.this))
                    }
                    d.prev <- d.this
                }
            },
            forward = function(input) {
                for(module in private$modules_) {
                    input <- module(input)
                }
                input
            })
    return(ret(in_, layers, generator, name))
}

