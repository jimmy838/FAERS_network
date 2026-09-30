# Numerical functions for association networks.

calculate_phi_network <- function(binary_matrix, significance_threshold = 0.01) {
  n <- nrow(binary_matrix)
  p <- ncol(binary_matrix)
  values <- matrix(0, p, p, dimnames = list(colnames(binary_matrix), colnames(binary_matrix)))
  p_values <- matrix(1, p, p)
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      x <- binary_matrix[, i]
      y <- binary_matrix[, j]
      n11 <- sum(x == 1 & y == 1)
      n10 <- sum(x == 1 & y == 0)
      n01 <- sum(x == 0 & y == 1)
      n00 <- sum(x == 0 & y == 0)
      denominator <- sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
      if (!is.finite(denominator) || denominator == 0) next
      phi <- (n11 * n00 - n10 * n01) / denominator
      values[i, j] <- values[j, i] <- max(phi, 0)
      p_values[i, j] <- p_values[j, i] <- pchisq(n * phi^2, df = 1, lower.tail = FALSE)
    }
  }
  cutoff <- significance_threshold / choose(p, 2L)
  values[p_values > cutoff] <- 0
  values
}

calculate_ppmi_network_bootstrap <- function(binary_matrix, k = 1,
                                               n_bootstrap = 1000,
                                               significance_threshold = 0.01) {
  n <- nrow(binary_matrix)
  p <- ncol(binary_matrix)
  marginal <- colSums(binary_matrix)
  observed <- crossprod(binary_matrix)
  smoothed_n <- n + k * p
  smoothed_marginal <- marginal + k
  pmi <- log2((observed + k) * smoothed_n / outer(smoothed_marginal, smoothed_marginal))
  diag(pmi) <- 0
  observed <- pmax(pmi, 0)
  null_values <- array(0, dim = c(p, p, n_bootstrap))
  for (b in seq_len(n_bootstrap)) {
    shuffled <- apply(binary_matrix, 2, sample)
    marginal_b <- colSums(shuffled) + k
    pmi_b <- log2((crossprod(shuffled) + k) * smoothed_n / outer(marginal_b, marginal_b))
    diag(pmi_b) <- 0
    null_values[, , b] <- pmax(pmi_b, 0)
  }
  p_values <- matrix(1, p, p)
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      if (observed[i, j] > 0) {
        p_values[i, j] <- p_values[j, i] <- mean(null_values[i, j, ] >= observed[i, j])
      }
    }
  }
  cutoff <- significance_threshold / choose(p, 2L)
  observed[p_values > cutoff] <- 0
  dimnames(observed) <- list(colnames(binary_matrix), colnames(binary_matrix))
  observed
}

graph_edge_table <- function(graph, network_name) {
  if (igraph::ecount(graph) == 0L) {
    return(data.frame(node_1 = character(), node_2 = character(), weight = numeric(), network = character()))
  }
  edges <- igraph::as_data_frame(graph, what = "edges")
  names(edges)[1:2] <- c("node_1", "node_2")
  edges$weight <- as.numeric(edges$weight)
  edges$network <- network_name
  edges[order(-edges$weight, edges$node_1, edges$node_2), c("node_1", "node_2", "weight", "network")]
}

largest_component <- function(graph) {
  if (igraph::vcount(graph) == 0L) return(NULL)
  component_ids <- igraph::components(graph)$membership
  largest <- which.max(igraph::components(graph)$csize)
  igraph::induced_subgraph(graph, which(component_ids == largest))
}

ring_lattice <- function(nodes, edges) {
  if (nodes < 2L || edges < 1L) return(igraph::make_empty_graph(nodes))
  pairs <- utils::combn(seq_len(nodes), 2L)
  distance <- abs(pairs[1, ] - pairs[2, ])
  distance <- pmin(distance, nodes - distance)
  keep <- head(order(distance, pairs[1, ], pairs[2, ]), min(edges, ncol(pairs)))
  igraph::graph_from_edgelist(t(pairs[, keep, drop = FALSE]), directed = FALSE)
}

small_worldness <- function(graph, simulations = 30L) {
  observed <- largest_component(graph)
  if (is.null(observed) || igraph::vcount(observed) < 3L || igraph::ecount(observed) < 2L) {
    return(NA_real_)
  }
  path_length <- suppressWarnings(igraph::mean_distance(observed, directed = FALSE,
                                                        unconnected = FALSE, weights = NA))
  clustering <- suppressWarnings(igraph::transitivity(observed, type = "globalundirected"))
  lattice <- ring_lattice(igraph::vcount(observed), igraph::ecount(observed))
  lattice_clustering <- suppressWarnings(igraph::transitivity(lattice, type = "globalundirected"))
  if (!is.finite(path_length) || !is.finite(clustering) || !is.finite(lattice_clustering) ||
      path_length <= 0 || lattice_clustering <= 0) return(NA_real_)
  random_lengths <- replicate(simulations, {
    random_graph <- igraph::sample_gnm(igraph::vcount(observed), igraph::ecount(observed),
                                       directed = FALSE, loops = FALSE)
    random_component <- largest_component(random_graph)
    if (is.null(random_component)) {
      NA_real_
    } else {
      suppressWarnings(igraph::mean_distance(random_component, directed = FALSE,
                                              unconnected = FALSE, weights = NA))
    }
  })
  random_lengths <- random_lengths[is.finite(random_lengths)]
  if (!length(random_lengths)) return(NA_real_)
  mean(random_lengths) / path_length - clustering / lattice_clustering
}

network_summary <- function(graph, network_name) {
  nodes <- igraph::vcount(graph)
  edges <- igraph::ecount(graph)
  degrees <- igraph::degree(graph)
  node_names <- igraph::V(graph)$name
  if (is.null(node_names)) node_names <- as.character(seq_len(nodes))
  cluster_labels <- igraph::V(graph)$cluster
  if (is.null(cluster_labels)) cluster_labels <- seq_len(nodes)
  data.frame(
    network = network_name,
    total_nodes = nodes,
    connected_nodes = sum(degrees > 0),
    isolated_nodes = sum(degrees == 0),
    total_edges = edges,
    density = if (nodes > 1L) edges / choose(nodes, 2L) else NA_real_,
    clusters = length(unique(cluster_labels)),
    modularity = if (edges > 0L) {
      membership <- igraph::V(graph)$cluster
      igraph::modularity(graph, membership)
    } else NA_real_,
    small_worldness = small_worldness(graph),
    central_node = if (nodes) node_names[which.max(degrees)] else NA_character_,
    maximum_degree = if (nodes) max(degrees) else NA_real_
  )
}

edge_jaccard <- function(graph_a, graph_b) {
  edges_a <- graph_edge_table(graph_a, "a")
  edges_b <- graph_edge_table(graph_b, "b")
  key_a <- paste(pmin(edges_a$node_1, edges_a$node_2), pmax(edges_a$node_1, edges_a$node_2), sep = "||")
  key_b <- paste(pmin(edges_b$node_1, edges_b$node_2), pmax(edges_b$node_1, edges_b$node_2), sep = "||")
  union_size <- length(union(key_a, key_b))
  if (!union_size) return(NA_real_)
  length(intersect(key_a, key_b)) / union_size
}

cluster_purity <- function(cluster_a, cluster_b) {
  common <- intersect(names(cluster_a), names(cluster_b))
  if (!length(common)) return(NA_real_)
  a <- cluster_a[common]
  b <- cluster_b[common]
  sum(vapply(unique(a), function(id) max(table(b[names(a)[a == id]])), numeric(1))) / length(common)
}

symmetric_purity <- function(cluster_a, cluster_b) {
  mean(c(cluster_purity(cluster_a, cluster_b), cluster_purity(cluster_b, cluster_a)))
}
