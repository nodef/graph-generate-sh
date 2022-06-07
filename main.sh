src="graph-generate"
out="/home/resources/Documents/subhajit/$src.log"
ulimit -s unlimited
printf "" > "$out"

# Download program
rm -rf $src
git clone https://github.com/puzzlef/$src
cd $src

# Run
g++ -std=c++17 -O3 main.cxx

log() {
  stdbuf --output=L "$@" 2>&1 | tee -a "$out"
}

# generate-delta <graph> <batch_size> <repeat>
generate-delta() {
  G="$1"; B="$2"; R="$3"; S=$(( 2*B ))
  cp ~/data/"${G}.mtx" ~/out/"${G}.mtx"
  log ./a.out rewrite  ~/out/"${G}.mtx" -t "unsymmetricize,set-weights" -o ~/out/"${G}_unsymmetric.mtx"
  log ./a.out rewrite  ~/out/"${G}.mtx"                 -o ~/out/"${G}.edges"
  log ./a.out delta    ~/out/"${G}.mtx" -s "$S" -c "$R" -o ~/out/"${G}_$B.delta"
}

# generate-apply <graph> <batch_size> <repeat>
generate-apply() {
  G="$1"; B="$2"; R="$3"
  for ((i=0; i<$R; i++)); do
    log ./a.out rewrite ~/out/"${G}.mtx" -t "--$HOME/out/${G}_$B-$i.delta,set-weights" -o ~/out/"${G}_$B-$i.mtx"
    log ./a.out rewrite ~/out/"${G}.mtx" -t	"++$HOME/out/${G}_$B+$i.delta,set-weights" -o ~/out/"${G}_$B+$i.mtx"
    log ./a.out rewrite ~/out/"${G}_$B-$i.mtx" -o ~/out/"${G}_$B-$i.edges"
    log ./a.out rewrite ~/out/"${G}_$B+$i.mtx" -o ~/out/"${G}_$B+$i.edges"
  done
}

# main <graph>
generate-batches() {
  log echo "# generate-batches $1"
  generate-delta "$1" 500   5
  generate-apply "$1" 500   5
  generate-delta "$1" 1000  5
  generate-apply "$1" 1000  5
  generate-delta "$1" 2000  5
  generate-apply "$1" 2000  5
  generate-delta "$1" 5000  5
  generate-apply "$1" 5000  5
  generate-delta "$1" 10000 5
  generate-apply "$1" 10000 5
}

# find-gini-coeff <graph>
find-gini-coeff() {
  G="$1"
  log echo "# find-gini-coeff $1"
  cp ~/data/"${G}.mtx" ~/out/"${G}.mtx"
  log nvgraph pagerank ~/out/"${G}.mtx" -o ~/out/"${G}.yaml"
  rm ~/out/"${G}.mtx"
  log ./a.out rewrite ~/data/"${G}.mtx" -t "loop-deadends" -o ~/out/"${G}_loop.mtx"
  log nvgraph pagerank ~/out/"${G}_loop.mtx" -o ~/out/"${G}_loop.yaml"
  rm ~/out/"${G}_loop.mtx"
  log ./a.out rewrite ~/data/"${G}.mtx" -t "loop-all" -o ~/out/"${G}_loopall.mtx"
  log nvgraph pagerank ~/out/"${G}_loopall.mtx" -o ~/out/"${G}_loopall.yaml"
  rm ~/out/"${G}_loopall.mtx"
}

main() {
  find-gini-coeff "$@"
}

main web-Stanford
main web-BerkStan
main web-Google
main web-NotreDame
main soc-Slashdot0811
main soc-Slashdot0902
main soc-Epinions1
main coAuthorsDBLP
main coAuthorsCiteseer
main soc-LiveJournal1
main coPapersCiteseer
main coPapersDBLP
main indochina-2004
main italy_osm
main great-britain_osm
main germany_osm
main asia_osm

main arabic-2005
main uk-2005
main it-2004
# main soc-Epinions1
# main soc-LiveJournal1
main wiki-Talk
main cit-Patents
# main coPapersDBLP
main amazon-2008
# main italy_osm
main Linux_call_graph
