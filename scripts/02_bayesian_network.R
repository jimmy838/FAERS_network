# Bootstrap-averaged Bayesian network estimation and numerical summaries.
source("config.R", encoding = "UTF-8")
source(file.path("R", "load_dependencies.R"), encoding = "UTF-8")
setwd(project_root)

drug <- analysis_drug
network_file <- file.path(output_dir, paste0("network_data_FAERS_", drug, ".rds"))
if (!file.exists(network_file)) {
  stop("Run 01_cooccurrence_network.R first: ", network_file)
}

network_data <- data.table::as.data.table(readRDS(network_file))
binary_long <- unique(network_data[, .(primaryid, pt)])
binary_wide <- tidyr::pivot_wider(
  binary_long[, value := 1],
  names_from = pt,
  values_from = value,
  values_fill = 0,
  values_fn = length
)
patient_ids <- binary_wide$primaryid
binary_matrix <- as.matrix(as.data.frame(binary_wide[, -1]))
rownames(binary_matrix) <- patient_ids
storage.mode(binary_matrix) <- "numeric"
binary_matrix[binary_matrix > 0] <- 1

original_names <- colnames(binary_matrix)
clean_names <- make.unique(make.names(original_names))
name_map <- data.frame(original = original_names, clean = clean_names)
bn_data <- as.data.frame(binary_matrix)
colnames(bn_data) <- clean_names
bn_data[] <- lapply(bn_data, factor, levels = c(0, 1))

death_clean <- name_map$clean[tolower(name_map$original) == "death"]
if (!length(death_clean)) stop("The required death node was not found.")
blacklist <- data.frame(
  from = death_clean,
  to = setdiff(clean_names, death_clean),
  stringsAsFactors = FALSE
)

set.seed(random_seed)
batch_size <- 50L
n_batches <- ceiling(n_bootstrap / batch_size)
bootstrap_results <- vector("list", n_batches)
progress <- txtProgressBar(min = 0, max = n_batches, style = 3, char = "=")
for (i in seq_len(n_batches)) {
  current_R <- min(batch_size, n_bootstrap - (i - 1L) * batch_size)
  bootstrap_results[[i]] <- bnlearn::boot.strength(
    data = bn_data,
    R = current_R,
    algorithm = "hc",
    algorithm.args = list(score = "bic", blacklist = blacklist)
  )
  setTxtProgressBar(progress, i)
}
close(progress)

bn_boot <- do.call(rbind, bootstrap_results) |>
  dplyr::group_by(from, to) |>
  dplyr::summarise(
    strength = mean(strength),
    direction = mean(direction),
    .groups = "drop"
  ) |>
  as.data.frame()
bn_boot$from <- as.character(bn_boot$from)
bn_boot$to <- as.character(bn_boot$to)
class(bn_boot) <- c("bn.strength", "data.frame")
attr(bn_boot, "nodes") <- names(bn_data)
attr(bn_boot, "method") <- "bootstrap"
attr(bn_boot, "threshold") <- 0.5
attr(bn_boot, "algorithm") <- list(algo = "hc")

saveRDS(bn_boot, file.path(output_dir, paste0("bn_bootstrap_", drug, ".rds")))
averaged <- bnlearn::averaged.network(bn_boot)
arcs_clean <- as.data.frame(bnlearn::arcs(averaged))
names(arcs_clean) <- c("from_clean", "to_clean")
strength_edges <- bn_boot[bn_boot$strength >= attr(averaged, "threshold"), ]
arcs_clean <- merge(
  arcs_clean,
  strength_edges[, c("from", "to", "strength", "direction")],
  by.x = c("from_clean", "to_clean"),
  by.y = c("from", "to"),
  all.x = TRUE
)
arcs_clean$from <- name_map$original[match(arcs_clean$from_clean, name_map$clean)]
arcs_clean$to <- name_map$original[match(arcs_clean$to_clean, name_map$clean)]
if (anyNA(arcs_clean$from) || anyNA(arcs_clean$to)) {
  stop("At least one Bayesian-network node could not be restored to its original name.")
}

write.csv(arcs_clean, file.path(output_dir, paste0("bn_edgelist_", drug, ".csv")), row.names = FALSE)
saveRDS(name_map, file.path(output_dir, paste0("bn_name_map_", drug, ".rds")))
saveRDS(
  as.matrix(arcs_clean[, c("from", "to")]),
  file.path(output_dir, paste0("net_bayesian_averaged_", drug, "_edge_list.rds"))
)

bn_graph <- igraph::graph_from_data_frame(
  arcs_clean[, c("from", "to", "strength")],
  directed = TRUE
)
if ("death" %in% igraph::V(bn_graph)$name) {
  stopifnot(length(igraph::neighbors(bn_graph, "death", mode = "out")) == 0L)
}
saveRDS(bn_graph, file.path(output_dir, paste0("bn_igraph_", drug, ".rds")))

out_degree <- igraph::degree(bn_graph, mode = "out")
in_degree <- igraph::degree(bn_graph, mode = "in")
centrality <- data.frame(
  node = igraph::V(bn_graph)$name,
  out_degree = unname(out_degree),
  in_degree = unname(in_degree),
  total_degree = unname(out_degree + in_degree),
  betweenness = round(igraph::betweenness(bn_graph, directed = TRUE, normalized = TRUE), 4),
  stringsAsFactors = FALSE
)
centrality$node_type <- ifelse(centrality$in_degree == 0, "source",
                              ifelse(centrality$out_degree == 0, "sink", "intermediate"))
centrality <- centrality[order(-centrality$total_degree, centrality$node), ]
write.csv(centrality, file.path(output_dir, paste0("bn_centrality_", drug, ".csv")), row.names = FALSE)

adjacent_rows <- lapply(head(centrality$node, 10), function(node) {
  adjacent <- names(igraph::neighbors(bn_graph, node, mode = "all"))
  data.frame(
    node = node,
    total_degree = length(adjacent),
    adjacent_nodes = paste(sort(adjacent), collapse = "; "),
    stringsAsFactors = FALSE
  )
})
adjacent <- do.call(rbind, adjacent_rows)
bn_statistics <- data.frame(
  metric = c("reports", "nodes", "edges", "density", "mean_edge_strength",
             "maximum_edge_strength", "bootstrap_iterations", "threshold"),
  value = c(
    uniqueN(network_data$primaryid),
    vcount(bn_graph),
    ecount(bn_graph),
    if (vcount(bn_graph) > 1L) ecount(bn_graph) / (vcount(bn_graph) * (vcount(bn_graph) - 1L)) else NA_real_,
    if (ecount(bn_graph)) mean(igraph::E(bn_graph)$strength) else NA_real_,
    if (ecount(bn_graph)) max(igraph::E(bn_graph)$strength) else NA_real_,
    n_bootstrap,
    attr(averaged, "threshold")
  )
)
write.csv(bn_statistics, file.path(output_dir, paste0("bn_network_statistics_", drug, ".csv")),
          row.names = FALSE)
summary_book <- openxlsx::createWorkbook()
openxlsx::addWorksheet(summary_book, "Centrality")
openxlsx::writeData(summary_book, "Centrality", centrality)
openxlsx::addWorksheet(summary_book, "Edges")
openxlsx::writeData(summary_book, "Edges", arcs_clean)
openxlsx::addWorksheet(summary_book, "Adjacent_nodes")
openxlsx::writeData(summary_book, "Adjacent_nodes", adjacent)
openxlsx::addWorksheet(summary_book, "Network_info")
openxlsx::writeData(summary_book, "Network_info", bn_statistics)
openxlsx::saveWorkbook(
  summary_book,
  file.path(output_dir, paste0("bn_summary_", drug, ".xlsx")),
  overwrite = TRUE
)
