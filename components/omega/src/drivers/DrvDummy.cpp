// OCN dummy driver

#include <iostream>

void dummy(int argc, char **argv);

int main(int argc, char **argv) {

   std::cout << std::endl << "Starting driver..." << std::endl;

   dummy(argc, argv);

   std::cout << "Stopped driver." << std::endl;

   return 0;
}
