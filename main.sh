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
stdbuf --output=L ./a.out -t loop ~/data/email-Eu-core-temporal.txt 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/CollegeMsg.txt             2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-mathoverflow.txt        2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-askubuntu.txt           2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-askubuntu-a2q.txt       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-askubuntu-c2a.txt       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-askubuntu-c2q.txt       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-superuser.txt           2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-superuser-a2q.txt       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-superuser-c2a.txt       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-superuser-c2q.txt       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/wiki-talk-temporal.txt     2>&1 | tee -a "$out"
stdbuf --output=L ./a.out -t loop ~/data/sx-stackoverflow.txt       2>&1 | tee -a "$out"
