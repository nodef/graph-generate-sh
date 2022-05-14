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

generate-delta() {
  log cp ~/data/"$1.mtx" ~/out/"$1.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t unsymmetricize,set-weights -o ~/out/"$1_unsymmetric.mtx"
  log ./a.out rewrite ~/out/"$1.mtx"                -o ~/out/"$1.edges"
  log ./a.out delta   ~/out/"$1.mtx" -s 10000 -c 5  -o ~/out/"$1.delta"
}

generate-apply-mtx() {
  log ./a.out rewrite ~/out/"$1.mtx" -t --$HOME/out/$1-0.delta,set-weights -o ~/out/"$1-0.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t	++$HOME/out/$1+0.delta,set-weights -o ~/out/"$1+0.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t	--$HOME/out/$1-1.delta,set-weights -o ~/out/"$1-1.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t	++$HOME/out/$1+1.delta,set-weights -o ~/out/"$1+1.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t	--$HOME/out/$1-2.delta,set-weights -o ~/out/"$1-2.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t	++$HOME/out/$1+2.delta,set-weights -o ~/out/"$1+2.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t	--$HOME/out/$1-3.delta,set-weights -o ~/out/"$1-3.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t	++$HOME/out/$1+3.delta,set-weights -o ~/out/"$1+3.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t	--$HOME/out/$1-4.delta,set-weights -o ~/out/"$1-4.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t	++$HOME/out/$1+4.delta,set-weights -o ~/out/"$1+4.mtx"
}

generate-apply-edges() {
  log ./a.out rewrite ~/out/"$1-0.mtx" -o ~/out/"$1-0.edges"
  log ./a.out rewrite ~/out/"$1+0.mtx" -o ~/out/"$1+0.edges"
  log ./a.out rewrite ~/out/"$1-1.mtx" -o ~/out/"$1-1.edges"
  log ./a.out rewrite ~/out/"$1+1.mtx" -o ~/out/"$1+1.edges"
  log ./a.out rewrite ~/out/"$1-2.mtx" -o ~/out/"$1-2.edges"
  log ./a.out rewrite ~/out/"$1+2.mtx" -o ~/out/"$1+2.edges"
  log ./a.out rewrite ~/out/"$1-3.mtx" -o ~/out/"$1-3.edges"
  log ./a.out rewrite ~/out/"$1+3.mtx" -o ~/out/"$1+3.edges"
  log ./a.out rewrite ~/out/"$1-4.mtx" -o ~/out/"$1-4.edges"
  log ./a.out rewrite ~/out/"$1+4.mtx" -o ~/out/"$1+4.edges"
}

main() {
  generate-apply-edges "$@"
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
