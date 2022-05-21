//#include "Action.h"
#include <iostream>

struct Extras_Proxy_ {
  Extras_Proxy_() {
    std::cout << "Registering extras\n";
  }
};

//declare
static Extras_Proxy_ extras_proxy_;

