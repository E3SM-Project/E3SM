// OCN dummy driver

#include <iostream>

using namespace std;

void dummy(int argc, char **argv);

int main(int argc, char **argv) {

   cout << "Starting driver..." << endl;

   dummy(argc, argv);

   cout << "Stopped driver." << endl;

   return 0;
}
