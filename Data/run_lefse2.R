pacman::p_load(tidyverse,
               microbiomeMarker)

devtools::source_url("https://github.com/HuaZou/MicrobiomeAnalysis/blob/main/R/utilities_lefse.R?raw=TRUE")
devtools::source_url("https://github.com/HuaZou/MicrobiomeAnalysis/blob/main/R/utilities.R?raw=TRUE")
devtools::source_url("https://github.com/yiluheihei/microbiomeMarker/blob/master/R/lefse-utilities.R?raw=TRUE")

run_lefse2 <-function(ps,
                      group,
                      subgroup = NULL,
                      taxa_rank = "all",
                      transform = c("identity",
                                    "log10", "log10p"),
                      norm = "CPM",
                      norm_para = list(),
                      kw_cutoff = 0.05,
                      lda_cutoff = 2,
                      bootstrap_n = 30,
                      bootstrap_fraction = 2 / 3,
                      wilcoxon_cutoff = 0.05,
                      multigrp_strat = FALSE,
                      strict = c("0",
                                 "1", "2"),
                      sample_min = 10,
                      only_same_subgrp = FALSE,
                      curv = FALSE)

{
  transform <- match.arg(transform, c("identity", "log10", 
                                      "log10p"))
  strict <- match.arg(strict, c("0", "1", "2"))
  strict <- as.numeric(strict)
  summarized <- check_tax_summarize(ps)
  if (summarized && norm != "CPM") {
    stop("`norm` must be a 'CPM' or 'none' while `ps` has been summarized", 
         call. = FALSE)
  }
  ps <- preprocess_ps(ps)
  ps <- transform_abundances(ps, transform = transform)
  norm_para <- c(norm_para, method = norm, object = list(ps))
  ps_normed <- do.call(normalize, norm_para)
  sample_meta <- sample_data(ps_normed)
  grp_info <- lefse_format_grp(sample_meta, group, subgroup = subgroup)
  grp <- grp_info$group
  subgrp <- grp_info$subgroup
  grp_hie <- grp_info$group_hie
  check_taxa_rank(ps, taxa_rank)
  if (taxa_rank == "all") {
    ps_summarized <- summarize_taxa(ps_normed)
  }
  else if (taxa_rank == "none") {
    ps_summarized <- extract_rank(ps_normed, taxa_rank)
  }
  else {
    ps_summarized <- aggregate_taxa(ps_normed, taxa_rank) %>% 
      extract_rank(taxa_rank)
  }
  otus <- abundances(ps_summarized, norm = FALSE)
  otus_test <- as.data.frame(t(otus), stringsAsFactors = FALSE)
  feature <- tax_table(ps_summarized)@.Data[, 1]
  names(otus_test) <- feature
  kw_p <- purrr::map_dbl(otus_test, ~kruskal.test(.x, grp)$p.value)
  na_ind <- is.na(kw_p)
  if (sum(na_ind) >= 1) {
    otus_test <- otus_test[!na_ind]
    kw_p <- kw_p[!na_ind]
  }
  sig_ind <- kw_p <= kw_cutoff
  sig_otus <- otus_test[, sig_ind]
  features_nms <- names(sig_otus)
  wilcoxon_p <- purrr::map2_lgl(sig_otus, features_nms, ~test_rep_wilcoxon(subgroup, 
                                                                           grp_hie, .x, .y, wilcoxon_cutoff = wilcoxon_cutoff, 
                                                                           multicls_strat = multigrp_strat, strict = strict, sample_min = sample_min, 
                                                                           only_same_subcls = only_same_subgrp, curv = curv))
  sig_otus <- sig_otus[, wilcoxon_p]
  otus_enriched_group <- get_feature_enrich_group(grp, sig_otus)
  ldas <- bootstap_lda(sig_otus, boot_n = bootstrap_n, class = grp, 
                       sample_fract = bootstrap_fraction)
  lefse_res <- data.frame(feature = names(sig_otus), enrich_group = otus_enriched_group$group, 
                          ef_lda = ldas, pvalue = kw_p[sig_ind][wilcoxon_p], stringsAsFactors = FALSE)
  lefse_sig <- filter(lefse_res, .data$ef_lda >= lda_cutoff) %>% 
    arrange(.data$enrich_group, desc(.data$ef_lda))
  lefse_out <- return_marker(lefse_sig, lefse_res)
  lefse_out$padj <- lefse_out$pvalue
  row.names(lefse_out) <- paste0("marker", seq_len(nrow(lefse_out)))
  tax <- matrix(feature) %>% tax_table()
  row.names(tax) <- row.names(otus)
  mm <- microbiomeMarker(marker_table = lefse_out, norm_method = get_norm_method(norm), 
                         diff_method = "lefse", otu_table = otu_table(otus, taxa_are_rows = TRUE), 
                         sam_data = sample_data(ps_normed), tax_table = tax)
  mm
}
