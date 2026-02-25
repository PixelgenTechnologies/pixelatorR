# Calculate edge weights for pMDS

**\[experimental\]**

## Usage

``` r
cos_distance_weights(g, pivots)

prob_distance_weights(g, k = 5, min_weight = 0.001)
```

## Arguments

- g:

  A `tbl_graph` object

- pivots:

  Number of pivots

- k:

  Number of steps

- min_weight:

  Minimum allowed weight when computing transition probabilities.

## Cosine weights

The main idea is that nodes in close spatial proximity should have
similar distance profiles to all other nodes in the graph. First, we
calculate the pairwise distances between all nodes and a set of randomly
selected pivot nodes. After a double centering transformation, we can
compute the cosine distance \\d(u, v)\\ between the distance profiles
for nodes \\u\\ and \\v\\ connected by an edge in the graph.

The cosine distance is a measure of similarity between two vectors and
ranges from -1 to 1, where 1 indicates maximal similarity and -1
indicates maximal dissimilarity. Since we want to assign higher weights
to nodes that are dissimilar, we convert the distances as follows:

\$\$w\_{u,v} = 2 - d(u, v)\$\$

## Transition probability weights

Weights are calculated for edges in the graph from the bidirectional
transition probability of a random walker. For an undirected edge
\\e\_{i,j}\\, we define the bidirectional transition probability
\\P(e\_{u,v}\\ as the probability of transitioning from \\u\\ to \\v\\
and back to \\u\\ in \\k\\ steps. Note that it doesn't matter what node
the random walker starts from since it has to go in both directions. In
this computation, we also allow self-loops in the graph to slow the
transition of the random walker.

Then, we assume that there is a relationship between the bidirectional
transition probability and an exponential function which decays with
distance:

\$\$P(e\_{u,v}) \propto e^{-d(u, v)}\$\$

by solving for \\d(u, v)\\, we can define an edge weight that is
proportional to distance:

\$\$w\_{u,v} = -\log(P(e\_{u,v}))\$\$
