Perform certain operations upon a fixed graph.

```bash
## REWRITE
## -------

# Rewrite MTX as plain edges file (no weights).
$ ./a.out rewrite ~/data/web-Google.mtx --output ~/data/web-Google.edges

# Rewrite MTX as plain edges file after making it symmetric.
$ ./a.out rewrite ~/data/web-Google.mtx --transform symmetricize --output ~/data/web-Google_symmetric.edges


# Rewrite MTX after making it symmetric and seting edge weights to 1.
$ ./a.out rewrite ~/data/web-Google.mtx --transform symmetricize,set-weights --output ~/data/web-Google_symmetric.mtx

# Generate MTX after adding self-loops to dead ends.
$ ./a.out rewrite ~/data/web-Google.mtx --transform loop-deadends --output ~/data/web-Google_nodeadends.mtx
```

<br>

```bash
## DELTA
## -----

# Generate a graph delta (batch) of 10000 edges (samples) with an equal mix of insertions and deletions.
$ ./a.out delta ~/data/web-Google.mtx --samples 10000

# Generate a graph delta and write to file web-Google(+/-).delta.
$ ./a.out delta ~/data/web-Google.mtx --samples 10000 --output web-Google.delta

# Generate 5 graph deltas and write to file web-Google(+/-).delta.
$ ./a.out delta ~/data/web-Google.mtx --samples 10000 --count 5 --output web-Google.delta
```
