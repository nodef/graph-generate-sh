#pragma once
#include <string>
#include <algorithm>

using std::string;
using std::remove;




// COUNT-LINES
// -----------
// For counting temporal edges.

int countLines(const char *x) {
  int a = 1;
  for (; *x; x++) {
    if (*x == '\r' || *x == '\n') a++;
    else if (*x == '\r' && *(x+1) == '\n') x++;
  }
  return a;
}

int countLines(const string& x) {
  return countLines(x.c_str());
}




// REMOVE-ALL
// ----------

void removeAll(string& x, char c) {
  auto ie = remove(x.begin(), x.end(), c);
  x.erase(ie, x.end());
}

string removeAll(const string& x, char c) {
  string a = x; removeAll(a, c);
  return a;
}
