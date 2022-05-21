#include "magres.h"

using namespace MagRes;

int main(int argc, const char **argv) {
  try {
    MagresFile magres;
    magres.parse_from_file(argv[1]);
    std::cout << magres;
  }
  catch (exception_t& exc) {
    std::cerr << "Parsing failed: " << exc << '\n';
    return 1;
  }
  catch (notmagres_exception_t&) {
    std::cerr << "Not a new format magres file\n";
    return 2;
  }
  return 0;
}
