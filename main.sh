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

main() {
  log cp ~/data/"$1.mtx" ~/out/"$1.mtx"
  log ./a.out rewrite ~/out/"$1.mtx" -t symmetricize,set-weights -o ~/out/"$1_symmetric.mtx"
  log ./a.out rewrite ~/out/"$1_symmetric.mtx"                   -o ~/out/"$1_symmetric.edges"
  log ./a.out delta   ~/out/"$1_symmetric.mtx" -s 10000 -c 5     -o ~/out/"$1_symmetric.delta"
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

# stdbuf --output=L ./a.out rewrite ~/data/web-Stanford.mtx        -t symmetricize,set-weights -o ~/data/web-Stanford_symmetric.mtx       2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/web-BerkStan.mtx        -t symmetricize,set-weights -o ~/data/web-BerkStan_symmetric.mtx       2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/web-Google.mtx          -t symmetricize,set-weights -o ~/data/web-Google_symmetric.mtx         2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/web-NotreDame.mtx       -t symmetricize,set-weights -o ~/data/web-NotreDame_symmetric.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/soc-Slashdot0811.mtx    -t symmetricize,set-weights -o ~/data/soc-Slashdot0811_symmetric.mtx   2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/soc-Slashdot0902.mtx    -t symmetricize,set-weights -o ~/data/soc-Slashdot0902_symmetric.mtx   2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/soc-Epinions1.mtx       -t symmetricize,set-weights -o ~/data/soc-Epinions1_symmetric.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/coAuthorsDBLP.mtx       -t symmetricize,set-weights -o ~/data/coAuthorsDBLP_symmetric.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/coAuthorsCiteseer.mtx   -t symmetricize,set-weights -o ~/data/coAuthorsCiteseer_symmetric.mtx  2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/soc-LiveJournal1.mtx    -t symmetricize,set-weights -o ~/data/soc-LiveJournal1_symmetric.mtx   2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/coPapersCiteseer.mtx    -t symmetricize,set-weights -o ~/data/coPapersCiteseer_symmetric.mtx   2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/coPapersDBLP.mtx        -t symmetricize,set-weights -o ~/data/coPapersDBLP_symmetric.mtx       2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/indochina-2004.mtx      -t symmetricize,set-weights -o ~/data/indochina-2004_symmetric.mtx     2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/italy_osm.mtx           -t symmetricize,set-weights -o ~/data/italy_osm_symmetric.mtx          2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/great-britain_osm.mtx   -t symmetricize,set-weights -o ~/data/great-britain_osm_symmetric.mtx  2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/germany_osm.mtx         -t symmetricize,set-weights -o ~/data/germany_osm_symmetric.mtx        2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/asia_osm.mtx            -t symmetricize,set-weights -o ~/data/asia_osm_symmetric.mtx           2>&1 | tee -a "$out"

# stdbuf --output=L ./a.out rewrite ~/data/arabic-2005.mtx       -t symmetricize,set-weights -o ~/data/arabic-2005_symmetric.mtx        2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/uk-2005.mtx           -t symmetricize,set-weights -o ~/data/uk-2005_symmetric.mtx            2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/it-2004.mtx           -t symmetricize,set-weights -o ~/data/web-Stanford_symmetric.mtx       2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/soc-Epinions1.mtx     -t symmetricize,set-weights -o ~/data/soc-Epinions1_symmetric.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/soc-LiveJournal1.mtx  -t symmetricize,set-weights -o ~/data/soc-LiveJournal1_symmetric.mtx   2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/wiki-Talk.mtx         -t symmetricize,set-weights -o ~/data/wiki-Talk_symmetric.mtx          2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/cit-Patents.mtx       -t symmetricize,set-weights -o ~/data/cit-Patents_symmetric.mtx        2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/coPapersDBLP.mtx      -t symmetricize,set-weights -o ~/data/coPapersDBLP_symmetric.mtx       2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/amazon-2008.mtx       -t symmetricize,set-weights -o ~/data/amazon-2008_symmetric.mtx        2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/italy_osm.mtx         -t symmetricize,set-weights -o ~/data/italy_osm_symmetric.mtx          2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out rewrite ~/data/Linux_call_graph.mtx  -t symmetricize,set-weights -o ~/data/Linux_call_graph_symmetric.mtx   2>&1 | tee -a "$out"

# stdbuf --output=L ./a.out -s 10000 ~/data/web-Stanford.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/web-Google.mtx        2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/indochina-2004.mtx    2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/uk-2002.mtx           2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/webbase-2001.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/twitter7.mtx          2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/sk-2005.mtx           2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/arabic-2005.mtx       2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/uk-2005.mtx           2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/it-2004.mtx           2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/soc-Epinions1.mtx     2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/soc-LiveJournal1.mtx  2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/wiki-Talk.mtx         2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/cit-Patents.mtx       2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/coPapersDBLP.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/amazon-2008.mtx       2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/italy_osm.mtx         2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/Linux_call_graph.mtx  2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/web-Stanford.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/web-Google.mtx        2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/indochina-2004.mtx    2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/uk-2002.mtx           2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/webbase-2001.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/twitter7.mtx          2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/sk-2005.mtx           2>&1 | tee -a "$out"
