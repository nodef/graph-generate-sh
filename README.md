Generate a batch of deletions/insertions of edges upon a fixed graph.

```bash
## With fixed graphs:
## -s: batch size (50-50 mix of deletions and insertions)
$ ./a.out -s 10000 ~/data/web-Stanford.mtx
$ ./a.out -s 10000 ~/data/web-Google.mtx
$ ...

# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 {}
# order: 281903 size: 2312669 {} (selfLoopDeadEnds)
# order: 281903 size: 2312669 {} (transposeWithDegree)
# - 115217 214128
# - 96560 153236
# - 244152 237905
# - 129538 226933
# ...
# + 274818 63930
# + 270085 147496
# + 150545 703
# + 133155 269970
# ...
#
# Loading graph /home/subhajit/data/web-Google.mtx ...
# ...
```
