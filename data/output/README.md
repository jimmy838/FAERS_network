# Network output files

The analysis scripts write generated files to `results/network/`.

## Co-occurrence networks

For each study drug, `01_cooccurrence_network.R` produces:

- `ising_network_<drug>.rds`, `phi_network_<drug>.rds`, and `ppmi_network_<drug>.rds`: weighted `igraph` objects with node names, edge weights, report counts, cluster labels, and node sizes.
- `<Network>_<drug>_adjacency.csv`: adjacency matrices.
- `<network>_edges_<drug>.csv`: edge lists with node names and weights.
- `node_statistics_<drug>.csv`: report counts and degree values.
- `cluster_membership_<drug>.csv`: cluster assignment for each node.
- `network_statistics_<drug>.csv`: node, edge, density, modularity, small-worldness, and central-node statistics.
- `network_similarity_<drug>.csv`: edge Jaccard similarity and mean cluster purity.
- `network_summary_<drug>.xlsx`: numerical summaries and edge tables in one workbook.

The `.rds` graph files are the main inputs for interactive visualization. CSV and XLSX files are intended for numerical review and reporting.

## Bayesian networks

`02_bayesian_network.R` produces:

- `bn_igraph_<drug>.rds`: the directed Bayesian network as an `igraph` object.
- `bn_bootstrap_<drug>.rds`: bootstrap edge strength and direction estimates.
- `bn_edgelist_<drug>.csv`: averaged-network edges with strength and direction.
- `bn_centrality_<drug>.csv`: in-degree, out-degree, total degree, betweenness, and node type.
- `bn_network_statistics_<drug>.csv`: network size, density, edge strength, and bootstrap settings.
- `bn_summary_<drug>.xlsx`: centrality, edges, adjacent nodes, and network information.

## Shiny visualization

Upload one of the following `.rds` graph files to the Shiny application:

```text
ising_network_tofacitinib.rds
phi_network_tofacitinib.rds
ppmi_network_tofacitinib.rds
ising_network_upadacitinib.rds
phi_network_upadacitinib.rds
ppmi_network_upadacitinib.rds
```

The application is available at [Posit Connect Cloud](https://01a00e5b-12e3-4cd6-3574-6cada8d91575.share.connect.posit.cloud/). Use the `.rds` graph objects for visualization; use CSV and XLSX files for numerical inspection.

Generated outputs are excluded from the public code package by default. Run the scripts locally to recreate them.
