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
stdbuf --output=L ./a.out -s 10000 ~/data/arabic-2005.mtx       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/uk-2005.mtx           2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/it-2004.mtx           2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/soc-Epinions1.mtx     2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/soc-LiveJournal1.mtx  2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/wiki-Talk.mtx         2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/cit-Patents.mtx       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/coPapersDBLP.mtx      2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/amazon-2008.mtx       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/italy_osm.mtx         2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -s 10000 ~/data/Linux_call_graph.mtx  2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/web-Stanford.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/web-Google.mtx        2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/indochina-2004.mtx    2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/uk-2002.mtx           2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/webbase-2001.mtx      2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/twitter7.mtx          2>&1 | tee -a "$out"
# stdbuf --output=L ./a.out -s 10000 ~/data/sk-2005.mtx           2>&1 | tee -a "$out"
