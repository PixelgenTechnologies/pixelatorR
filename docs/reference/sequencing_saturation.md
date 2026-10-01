# Calculate Sequencing Saturation

This function calculates the sequencing saturation of a sample or a
graph component, which can be applied to unique reads, nodes, or edges.

## Usage

``` r
sequencing_saturation(graph_elements, graph_reads)
```

## Arguments

- graph_elements:

  The number of graph elements (e.g., nodes, edges).

- graph_reads:

  The number of reads or supporting elements.

## Value

The sequencing saturation of the graph expressed as a percentage.

## Details

The sequencing saturation is calculated using the formula: \$\$S = 100
\times \left(1 - \frac{E}{R}\right)\$\$ where:

- \\S\\ is the saturation (as a percentage),

- \\E\\ is the number of unique reads or graph elements (e.g., nodes,
  edges), and

- \\R\\ is the number of reads or supporting elements in the graph.

This function can be used for calculating the saturation of:

- Reads: Number of deduplicated reads vs. reads

- Nodes: Number of nodes vs. reads

- Edges: Number of edges vs. reads

- Other graph elements: Adjust the input accordingly

## Examples

``` r
# For a graph with 300 unique reads, 100 nodes, and 200 edges,
# sequenced at 400 total reads

# Read sequencing saturation
sequencing_saturation(300, 400)
#> [1] 25

# Node sequencing saturation
sequencing_saturation(100, 400)
#> [1] 75

# Edge sequencing saturation
sequencing_saturation(200, 400)
#> [1] 50
```
