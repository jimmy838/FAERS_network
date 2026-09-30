# Estimate numerical co-occurrence networks from prepared FAERS input.

source("config.R", encoding = "UTF-8")
source(file.path("R", "load_dependencies.R"), encoding = "UTF-8")
source(file.path("R", "network_functions.R"), encoding = "UTF-8")
setwd(project_root)

input_file <- file.path(input_dir, paste0("network_input_", analysis_drug, ".rds"))
input <- readRDS(input_file)
required_input <- c("results_pt", "reac", "drug_pids", "procedure_pts")
if (!all(required_input %in% names(input))) {
  stop("The input RDS must contain: ", paste(required_input, collapse = ", "))
}

results_pt <- data.table::as.data.table(input$results_pt)
reac <- data.table::as.data.table(input$reac)
drug_pids <- input$drug_pids
procedure_pts <- input$procedure_pts
drug <- tolower(analysis_drug)

signal_pts <- results_pt[
  IC_lower > 0 & (ROR_lower > 1 | (PRR_median >= 2 & chi_square >= 4) | EB05 > 2),
  event
]
case_threshold <- switch(
  drug,
  tofacitinib = ceiling(length(drug_pids) / 200),
  upadacitinib = ceiling(length(drug_pids) / 300),
  stop("Unsupported study drug: ", drug)
)

event_counts <- reac[primaryid %in% drug_pids & pt %in% signal_pts,
                     .(report_count = uniqueN(primaryid)), by = pt]
network_pts <- event_counts[report_count >= case_threshold, pt]
network_pts <- network_pts[!tolower(network_pts) %in% tolower(procedure_pts)]
network_data <- unique(reac[primaryid %in% drug_pids & pt %in% network_pts,
                            .(primaryid, pt)])
if (nrow(network_data) == 0L) stop("No PTs passed the network filters.")

saveRDS(network_data, file.path(output_dir, paste0("network_data_FAERS_", drug, ".rds")))
pt_counts <- network_data[, .(report_count = uniqueN(primaryid)), by = pt]
saveRDS(pt_counts, file.path(output_dir, paste0("pt_report_counts_", drug, ".rds")))

wide <- tidyr::pivot_wider(
  network_data[, value := 1], names_from = pt, values_from = value,
  values_fill = 0, values_fn = length
)
report_ids <- wide$primaryid
binary <- as.matrix(as.data.frame(wide[, -1]))
rownames(binary) <- report_ids
storage.mode(binary) <- "numeric"
binary[binary > 0] <- 1

set.seed(random_seed)
ising <- IsingFit::IsingFit(binary, progressbar = TRUE, plot = FALSE)$weiadj
ising[ising < 0] <- 0
phi <- calculate_phi_network(binary)
ppmi <- calculate_ppmi_network_bootstrap(binary, n_bootstrap = n_bootstrap)

graphs <- list(
  Ising = igraph::graph_from_adjacency_matrix(ising, mode = "undirected", weighted = TRUE),
  Phi = igraph::graph_from_adjacency_matrix(phi, mode = "undirected", weighted = TRUE),
  PPMI = igraph::graph_from_adjacency_matrix(ppmi, mode = "undirected", weighted = TRUE)
)

annotate_clusters <- function(graph) {
  if (igraph::vcount(graph) == 0L) return(graph)
  if (igraph::ecount(graph) == 0L) {
    igraph::V(graph)$cluster <- seq_len(igraph::vcount(graph))
  } else {
    communities <- igraph::cluster_fast_greedy(graph)
    igraph::V(graph)$cluster <- igraph::membership(communities)
  }
  igraph::V(graph)$size <- igraph::degree(graph)
  graph
}

for (name in names(graphs)) {
  graph <- graphs[[name]]
  graph <- igraph::delete_edges(graph, igraph::E(graph)[igraph::E(graph)$weight <= 0])
  graph <- annotate_clusters(graph)
  graphs[[name]] <- graph
  igraph::V(graph)$report_count <- pt_counts$report_count[
    match(igraph::V(graph)$name, pt_counts$pt)
  ]
  saveRDS(graph, file.path(output_dir, paste0(tolower(name), "_network_", drug, ".rds")))
  write.csv(
    igraph::as_adjacency_matrix(graph, sparse = FALSE, attr = "weight"),
    file.path(output_dir, paste0(name, "_", drug, "_adjacency.csv")),
    row.names = TRUE
  )
  write.csv(
    graph_edge_table(graph, name),
    file.path(output_dir, paste0(tolower(name), "_edges_", drug, ".csv")),
    row.names = FALSE
  )
}

node_stats <- data.frame(
  node = colnames(binary),
  report_count = pt_counts$report_count[match(colnames(binary), pt_counts$pt)],
  Ising_degree = igraph::degree(graphs$Ising)[colnames(binary)],
  Phi_degree = igraph::degree(graphs$Phi)[colnames(binary)],
  PPMI_degree = igraph::degree(graphs$PPMI)[colnames(binary)],
  row.names = NULL
)
write.csv(node_stats, file.path(output_dir, paste0("node_statistics_", drug, ".csv")),
          row.names = FALSE)

cluster_table <- data.frame(
  node = colnames(binary),
  Ising_cluster = igraph::V(graphs$Ising)$cluster,
  Phi_cluster = igraph::V(graphs$Phi)$cluster,
  PPMI_cluster = igraph::V(graphs$PPMI)$cluster,
  stringsAsFactors = FALSE
)
write.csv(cluster_table, file.path(output_dir, paste0("cluster_membership_", drug, ".csv")),
          row.names = FALSE)

network_stats <- rbind(
  network_summary(graphs$Ising, "Ising"),
  network_summary(graphs$Phi, "Phi"),
  network_summary(graphs$PPMI, "PPMI")
)
write.csv(network_stats, file.path(output_dir, paste0("network_statistics_", drug, ".csv")),
          row.names = FALSE)

membership <- lapply(graphs, function(graph) {
  labels <- igraph::V(graph)$cluster
  names(labels) <- igraph::V(graph)$name
  labels
})
network_similarity <- data.frame(
  network_pair = c("Ising-Phi", "Ising-PPMI", "Phi-PPMI"),
  edge_jaccard = c(
    edge_jaccard(graphs$Ising, graphs$Phi),
    edge_jaccard(graphs$Ising, graphs$PPMI),
    edge_jaccard(graphs$Phi, graphs$PPMI)
  ),
  mean_cluster_purity = c(
    symmetric_purity(membership$Ising, membership$Phi),
    symmetric_purity(membership$Ising, membership$PPMI),
    symmetric_purity(membership$Phi, membership$PPMI)
  )
)
write.csv(network_similarity, file.path(output_dir, paste0("network_similarity_", drug, ".csv")),
          row.names = FALSE)

summary_book <- openxlsx::createWorkbook()
openxlsx::addWorksheet(summary_book, "Network_statistics")
openxlsx::writeData(summary_book, "Network_statistics", network_stats)
openxlsx::addWorksheet(summary_book, "Node_statistics")
openxlsx::writeData(summary_book, "Node_statistics", node_stats)
openxlsx::addWorksheet(summary_book, "Cluster_membership")
openxlsx::writeData(summary_book, "Cluster_membership", cluster_table)
openxlsx::addWorksheet(summary_book, "Network_similarity")
openxlsx::writeData(summary_book, "Network_similarity", network_similarity)
for (name in names(graphs)) {
  openxlsx::addWorksheet(summary_book, paste0(name, "_edges"))
  openxlsx::writeData(summary_book, paste0(name, "_edges"), graph_edge_table(graphs[[name]], name))
}
openxlsx::saveWorkbook(
  summary_book,
  file.path(output_dir, paste0("network_summary_", drug, ".xlsx")),
  overwrite = TRUE
)
