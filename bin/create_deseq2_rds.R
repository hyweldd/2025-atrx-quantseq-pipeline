#!/usr/bin/env Rscript
# Written by Hywel Dunn-Davies, 2025-02-18


# Setup ----

RPROJECT_DIR <- "Rproject_mbp"
renv::load(project = RPROJECT_DIR)

library(optparse)

source("functions.R")


# Option parsing ----

option_list <- list(
    make_option(c("-i", "--count_file"    ), type="character", default=NULL    , metavar="path"   , help="Count file matrix where rows are genes and columns are samples."),
    make_option(c("-c", "--control"       ), type="character", default=''      , metavar="string" , help="Control sample prefix."),
    make_option(c("-t", "--treatment"     ), type="character", default=''      , metavar="string" , help="Treatment sample prefix."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'    , metavar="path"   , help="Output directory.")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$count_file)){
    print_help(opt_parser)
    stop("Please provide a counts file.", call.=FALSE)
}

if (opt$control == ''){
    print_help(opt_parser)
    stop("Please provide a treatment.", call.=FALSE)

}

if (opt$control == ''){
    print_help(opt_parser)
    stop("Please provide a control.", call.=FALSE)
}

assert_that(opt$control != opt$treatment)


# Main code ----

output_fn <- opt_to_output_fn(opt)
output_fp <- file.path(opt$outdir, output_fn)

if (file.exists(output_fp)){
    stop("Output file already exists.", call.=FALSE)
}

deseq_results <- opt_to_deseq_results(opt)

saveRDS(deseq_results$dds, output_fp)

write_csv(deseq_results$combined_results_table, file = paste0(output_fp, ".deseq2_results.csv"))
