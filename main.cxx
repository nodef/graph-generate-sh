#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include "inc/main.hxx"
#include "options.hxx"

using namespace std;




#pragma region MAIN
/**
 * Main function.
 * @param argc argument count
 * @param argv argument values
 * @returns zero on success, non-zero on failure
 */
int main(int argc, char **argv) {
  Options o = readOptions(argc, argv);
  printf("%s\n", helpMessage());
  return 0;
}
#pragma endregion
