List a few graph properties.

```bash
## With temporal graphs:
$ ./a.out -t loop ~/data/email-Eu-core-temporal.txt
$ ./a.out -t loop ~/data/CollegeMsg.txt
$ ...

# Using graph /home/subhajit/data/email-Eu-core-temporal.txt ...
# - temporal-edges: 332335 (including repeated ones)
# ---
# order: 751 size: 6952 {}
# order: 751 size: 7089 {} (loopDeadEnds)
# order: 751 size: 7089 {} (transposeWithDegree)
# - sample:  1 (of 10)
# - order:   751 (vertices)
# - size:    7089 (edges)
# - span:    1005 (vertex-id limit)
# - density: 1.256913e-02 (fully connectedness fraction)
# - degree-min:   1 (minimum out-degree)
# - degree-max:   89 (maximum out-degree)
# - degree-avg:   9.439414 (average out-degree)
# - indegree-min: 0 (minimum in-degree)
# - indegree-max: 86 (maximum in-degree)
# - indegree-avg: 9.439414 (average in-degree)
# - dead-ends:      0 (vertices with no outgoing edges)
# - dead-ends-rec:  0 (vertices pointing to dead ends)
# - beginnings:     23 (vertices with no incoming edges)
# - beginnings-rec: 26 (vertices pointed by beginnings)
# - self-loops:     137 (vertices pointing to self)
# - components:    168 (strongly connected)
# - component-min: 1 (minimum component size)
# - component-max: 582 (maximum component size)
# - component-avg: 4.470238 (average component size)
# - inidenticals:    4 (vertices with identical in-edges)
# - inidentical-min: 2 (minimum in-identical size)
# - inidentical-max: 2 (maximum in-identical size)
# - inidentical-avg: 2.000000 (average in-identical size)
# - outidenticals:    13 (vertices with identical out-edges)
# - outidentical-min: 2 (minimum out-identical size)
# - outidentical-max: 2 (maximum out-identical size)
# - outidentical-avg: 2.000000 (average out-identical size)
# - chains:    4 (vertices connected in a line)
# - chain-min: 1 (minimum chain size)
# - chain-max: 2 (maximum chain size)
# - chain-avg: 1.250000 (average chain size)
# - blockgraph-order:   168 (vertices)
# - blockgraph-size:    170 (edges)
# - blockgraph-span:    168 (vertex-id limit)
# - blockgraph-density: 6.023243e-03 (fully connectedness fraction)
# - blockgraph-degree-min:   0 (minimum out-degree)
# - blockgraph-degree-max:   134 (maximum out-degree)
# - blockgraph-degree-avg:   1.011905 (average out-degree)
# - blockgraph-indegree-min: 0 (minimum in-degree)
# - blockgraph-indegree-max: 22 (maximum in-degree)
# - blockgraph-indegree-avg: 1.011905 (average in-degree)
# - blockgraph-dead-ends:      139 (vertices with no outgoing edges)
# - blockgraph-dead-ends-rec:  168 (vertices pointing to dead ends)
# - blockgraph-beginnings:     24 (vertices with no incoming edges)
# - blockgraph-beginnings-rec: 168 (vertices pointed by beginnings)
# - blockgraph-self-loops:     0 (vertices pointing to self)
# - blockgraph-depth-min: 1 (minimum path length from beginnings)
# - blockgraph-depth-max: 5 (maximum path length from beginnings)
# - blockgraph-depth-avg: 3.791667 (average path length from beginnings)
# ---
# order: 802 size: 10532 {}
# order: 802 size: 10677 {} (loopDeadEnds)
# order: 802 size: 10677 {} (transposeWithDegree)
# - sample:  2 (of 10)
# - order:   802 (vertices)
# - size:    10677 (edges)
# - span:    1005 (vertex-id limit)
# - density: 1.659971e-02 (fully connectedness fraction)
# ...
```
