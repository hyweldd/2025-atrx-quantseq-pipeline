# Functions written by hywel, 18/02/2025

library(assertthat)
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(purrr)
library(DESeq2)


# Functions ----

load_count_tables <- function(opt, round_counts = TRUE) {
  
  if (is.null(opt$count_file)){
    print_help(opt_parser)
    stop("Please provide a counts file.", call.=FALSE)
  }
  
  assert_that(opt$control != opt$treatment)
  
  count_table <- read_tsv(file=opt$count_file, col_names=TRUE) %>% 
    mutate(gene_uuid=paste(gene_id, gene_name, sep=":::")) %>% 
    select(-gene_id, -gene_name)
  
  control_table <- count_table %>%
    select(gene_uuid, starts_with(opt$control))
  
  assert_that(ncol(control_table)>1)
  
  print(control_table %>% names)
  
  print(opt$treatment)
  
  treatment_table <- count_table %>%
    select(gene_uuid, starts_with(opt$treatment))
  
  assert_that(ncol(treatment_table)>1)
  
  print(treatment_table %>% names)
  
  count_tables <- list(trt = treatment_table, control = control_table)
  
  if(round_counts == TRUE) {
    count_tables <- lapply(count_tables, function(x) mutate(x, across(-gene_uuid, round))) 
  }
  
  count_tables
  
}

create_countdata <- function(trt_ctrl_tables) {
  
  assert_that(nrow(trt_ctrl_tables$control)==nrow(trt_ctrl_tables$trt))
  
  countdata <- trt_ctrl_tables$control %>%
    inner_join(trt_ctrl_tables$trt) %>%
    na.omit
  
  assert_that(nrow(countdata)==nrow(trt_ctrl_tables$control))
  
  countdata_mtx <- countdata %>%
    column_to_rownames("gene_uuid") %>%
    as.matrix
  
}

create_pairwise_design_matrix <- function(trt_ctrl_tables) {
  
  control_dm <- tibble(
    sample = names(trt_ctrl_tables$control %>% select(-gene_uuid)),
    condition = "control"
  ) 
  
  trt_dm <- tibble(
    sample = names(trt_ctrl_tables$trt %>% select(-gene_uuid)),
    condition = "treatment"
  ) 
 
  dm_df <- bind_rows(control_dm, trt_dm)
  
  dm_df$condition <- dm_df$condition %>% as.factor %>% relevel("control")
  
  dm_df
   
}

create_pairwise_dds <- function(coldata, countdata) {
  
  ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                   colData = coldata,
                                   design = ~ condition)
  
}

filter_dds <- function(dds, smallestGroupSize = 3) {
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
}

run_deseq2 <- function(dds) {
  dds <- DESeq(dds)
}

create_combined_results_table <- function(dds) {
  
  raw_counts <- counts(dds, normalized=FALSE) %>% as.data.frame %>% 
    rownames_to_column("gene_uuid") %>% 
    rename_with(~ paste0("raw_count_", .x, recycle0 = TRUE), -gene_uuid)
    
  norm_counts <- counts(dds, normalized=TRUE) %>% as.data.frame %>% 
    rownames_to_column("gene_uuid") %>% 
    rename_with(~ paste0("normalized_count_", .x, recycle0 = TRUE), -gene_uuid)
  
  res_df <- results(dds) %>% as.data.frame %>% 
    rownames_to_column("gene_uuid")
  
  crt_df <- res_df %>% 
    inner_join(norm_counts) %>%
    inner_join(raw_counts)
  
  assert_that(nrow(res_df)==nrow(crt_df))
  
  crt_df %>% 
    arrange(padj) %>%
    mutate(
      ensembl_gene_id=gsub(':::.*$', '', gene_uuid), 
      gene_symbol=gsub('^.*:::', '', gene_uuid)
    ) %>% 
    select(ensembl_gene_id, gene_symbol, everything()) %>% 
    select(-gene_uuid) %>% 
    return

}

opt_to_deseq_results <- function(opt) {
  
  trt_ctrl_tables <- load_count_tables(opt, round_counts = TRUE)
  
  countdata <- trt_ctrl_tables %>% create_countdata
  coldata <- trt_ctrl_tables %>% create_pairwise_design_matrix
  
  dds <- create_pairwise_dds(coldata = coldata, countdata = countdata)
  
  smallest_group_size <- min(
    ncol(trt_ctrl_tables$trt)-1,
    ncol(trt_ctrl_tables$control)-1
  )
  
  dds <- filter_dds(dds = dds, smallestGroupSize = smallest_group_size)
  
  dds <- run_deseq2(dds)
  
  crt_df <- create_combined_results_table(dds = dds)
  
  list(
    dds = dds,
    combined_results_table = crt_df
  ) %>% return
  
}

opt_to_output_fn <- function(opt) {
  
  paste0("trt_", opt$treatment, "-vs-ctrl_", opt$control, ".rds") %>% 
    return
  
}

get_deseq2_results_tables <- function() {
  
  ISG_444_TSS_bed <- read_table("ISG_444_TSS.bed.txt", col_names = FALSE) %>% 
    select(X4)
  
  names(ISG_444_TSS_bed) <- c("gene_symbol")
  
  deseq2_results_table_dir <- "deseq2_results_tables"

  deseq2_results_tables <- list.files(deseq2_results_table_dir, full.names = FALSE) %>%
    set_names %>%
    as.list %>%
    map(~read_csv(file.path(deseq2_results_table_dir, .))) %>%
    map(function(x) { select(x, gene_symbol, log2FoldChange, pvalue, padj) }) %>% 
    map(function(x) { ISG_444_TSS_bed %>% left_join(x) %>% mutate(log2FoldChange=ifelse(is.na(log2FoldChange),0,log2FoldChange)) })

}

create_interferon_genes_heatmap_simple <- function(deseq2_results_tables = get_deseq2_results_tables()) {
  
  background_df <- deseq2_results_tables[["trt_ATRXKI_CONTROL-vs-ctrl_ATRXWT_CONTROL.rds.deseq2_results.csv"]] %>% 
    select(gene_symbol, log2FoldChange) %>% 
    arrange(desc(log2FoldChange)) %>% 
    dplyr::rename(cc=log2FoldChange)
  
  gg_df <- deseq2_results_tables[["trt_ATRXKI_cGASKO-vs-ctrl_ATRXWT_cGASKO.rds.deseq2_results.csv"]] %>%
    select(gene_symbol, log2FoldChange) %>%
    arrange(desc(log2FoldChange)) %>%
    dplyr::rename(gg=log2FoldChange)

  data_df <- background_df %>%
    left_join(gg_df) %>%
    tidyr::pivot_longer(cols = cc:gg, names_to = "nm", values_to = "val") %>% 
    mutate(nm = as.factor(nm)) %>% 
    left_join(background_df)
  
  plt <- data_df %>% head(50) %>%
    ggplot(
      aes(
        x = c("val"),
        y = reorder(gene_symbol, cc),
        fill = val
      )
    ) + geom_tile() + facet_grid(~nm)
  
}


simplify_results_df <- function(results_df, filter_padj=NULL, filter_pval=NULL) {
  
  if(!is.null(filter_padj)) {
    results_df <- results_df %>% 
      mutate(
        log2FoldChange = ifelse(padj<filter_padj,log2FoldChange,NA)
      ) 
  }
  
  if(!is.null(filter_pval)) {
    results_df <- results_df %>% 
      mutate(
        log2FoldChange = ifelse(pvalue<filter_pval,log2FoldChange,NA)
      ) 
  }
  
  results_df <- results_df %>% 
    select(gene_symbol, log2FoldChange) %>% 
    arrange(desc(log2FoldChange))
}

get_experiment_name_mappings <- function(deseq2_results_tables) {
  
  deseq2_results_tables %>% 
    names() %>% 
    set_names() %>% 
    map(function(x) {
        str_replace_all(x, "(.rds|.deseq2_results.csv|trt_|ctrl_)", "") %>% 
        str_replace_all("[_-]", " ") %>% 
        str_replace_all(' vs', '\nvs')
      }) %>% 
    return
  
}

get_simple_results_df <- function(table_key, deseq2_results_tables, filter_padj, filter_pval) {
  experiment_name_mapping <- (get_experiment_name_mappings(deseq2_results_tables))[[table_key]]
  simple_results_df <- deseq2_results_tables[[table_key]] %>% 
    simplify_results_df(filter_padj = filter_padj, filter_pval = filter_pval)
  names(simple_results_df) <- c("gene_symbol", experiment_name_mapping)
  simple_results_df
    #dplyr::rename(function(x) return(experiment_name_mapping)=log2FoldChange) %>% return
}


create_interferon_genes_heatmap_source_table <- function(deseq2_results_tables = get_deseq2_results_tables(), filter_padj=NULL, filter_pval=NULL, sig_genes=NULL) {
  
  gsrdf <- partial(get_simple_results_df, 
    deseq2_results_tables = deseq2_results_tables, 
    filter_padj = filter_padj,
    filter_pval = filter_pval
  )
  
  experiments_to_plot <- c(
    "trt_ATRXKI_CONTROL-vs-ctrl_ATRXWT_CONTROL.rds.deseq2_results.csv",
    "trt_ATRXKI_cGASKO-vs-ctrl_ATRXWT_cGASKO.rds.deseq2_results.csv",
    "trt_ATRXKI_STING1KO-vs-ctrl_ATRXWT_STING1KO.rds.deseq2_results.csv",
    "trt_ATRXKI_IRF3KO-vs-ctrl_ATRXWT_IRF3KO.rds.deseq2_results.csv",
    "trt_ATRXKI_ANIFRO-vs-ctrl_ATRXKI_WATER.rds.deseq2_results.csv",
    "trt_ATRXWT_cGASKO-vs-ctrl_ATRXWT_CONTROL.rds.deseq2_results.csv",
    "trt_ATRXWT_STING1KO-vs-ctrl_ATRXWT_CONTROL.rds.deseq2_results.csv",
    "trt_ATRXWT_IRF3KO-vs-ctrl_ATRXWT_CONTROL.rds.deseq2_results.csv",
    "trt_ATRXWT_ANIFRO-vs-ctrl_ATRXWT_WATER.rds.deseq2_results.csv"
  )
  
  dfs <- experiments_to_plot %>%
    set_names %>% 
    map(gsrdf)
  
  experiment_name_mappings <- get_experiment_name_mappings(deseq2_results_tables = deseq2_results_tables)
  factor_levels <- experiments_to_plot %>% set_names() %>% map_chr(~experiment_name_mappings[[.]])
  
  background_df <- deseq2_results_tables[["trt_ATRXKI_CONTROL-vs-ctrl_ATRXWT_CONTROL.rds.deseq2_results.csv"]] %>%
    simplify_results_df %>%
    dplyr::rename('control_log2FC'=log2FoldChange)
 
  
  data_df <- dfs %>% 
    reduce(left_join)
  
  if(!is.null(sig_genes)) {
    data_df <- data_df %>% filter(gene_symbol %in% sig_genes)
  }
  
  data_df <- data_df %>% 
    tidyr::pivot_longer(cols = -gene_symbol, names_to = "nm", values_to = "log2FC") %>% 
    mutate(nm = as.factor(nm) %>% fct_relevel(factor_levels)) %>% 
    left_join(background_df)
  
}

get_sig_genes_pval <- function(deseq2_results_tables = get_deseq2_results_tables()) {
  
  controls_df <- deseq2_results_tables$`trt_ATRXKI_CONTROL-vs-ctrl_ATRXWT_CONTROL.rds.deseq2_results.csv`
  controls_df <- controls_df %>% 
    filter(!is.na(pvalue)) %>% 
    filter(pvalue<0.05)

  controls_df %>% select(gene_symbol) %>% 
    unlist %>% 
    return
  
}

get_sig_genes_padj <- function(deseq2_results_tables = get_deseq2_results_tables()) {
  
  controls_df <- deseq2_results_tables$`trt_ATRXKI_CONTROL-vs-ctrl_ATRXWT_CONTROL.rds.deseq2_results.csv`
  controls_df <- controls_df %>% 
    filter(!is.na(padj)) %>% 
    filter(padj<0.1)

  controls_df %>% select(gene_symbol) %>% 
    unlist %>% 
    return
  
}

create_interferon_genes_heatmap <- function(data_df, deseq2_results_tables = get_deseq2_results_tables()) { 
  
  data_df %>% 
    ggplot(aes(x=nm, y=reorder(gene_symbol, control_log2FC), fill = log2FC)) + 
      geom_tile() + 
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank()
      ) + 
      scale_fill_gradient2(low = "red", 
                          high = "green", 
                          mid = "white",
                          # Make the gradient2 "switch" at a midpoint of 0
                          midpoint = 0)
  
  
}

create_interferon_genes_pheatmap <- function(data_df, plot_title="Heatmap") {

  ph_data_df <- data_df %>% 
    select(-control_log2FC) %>% 
    pivot_wider(names_from = "nm", values_from = "log2FC") %>% 
    column_to_rownames("gene_symbol")
  
  heatmap_plot <- pheatmap::pheatmap(ph_data_df %>% select(1:5), 
                        cluster_rows = T, cluster_cols = F, 
                        clustering_distance_cols = 'euclidean',
                        clustering_distance_rows = 'euclidean',
                        clustering_method = 'ward.D',
                        annotation_names_row = T, 
                        fontsize_row = 10,
                        fontsize_col = 7,
                        angle_col = 45,
                        legend_breaks = c(-6, -4, -2, 0, 2, 4, 6, 8), 
                        legend_labels = c("-6", "-4", "-2", "0", "2", "4", "6", "log2FC"), 
                        show_colnames = T, show_rownames = T, 
                        main = plot_title)
}

create_clustered_heatmaps <- function() {
  
  pheatmap_padj <- create_interferon_genes_heatmap_source_table(
    deseq2_results_tables = get_deseq2_results_tables(),
    sig_genes = get_sig_genes_padj(),
    filter_padj = 0.1
  ) %>% create_interferon_genes_pheatmap(plot_title = "Clustered heatmap (adjusted p-value < 0.1)")
  
  pheatmap_pval <- create_interferon_genes_heatmap_source_table(
    deseq2_results_tables = get_deseq2_results_tables(),
    sig_genes = get_sig_genes_pval(),
    filter_pval = 0.05
  ) %>% create_interferon_genes_pheatmap(plot_title = "Clustered heatmap (p-value < 0.05)")
  
  list(
    padj=pheatmap_padj,
    pval=pheatmap_pval
  ) %>% return
  
}
