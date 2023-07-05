#include <iostream>

using namespace std;

double dummy(int argc, char **argv);

///////////////////////////////////////////////////////////////////////////////////////
// THE MAIN PROGRAM STARTS HERE
///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

   int retval = 0;

   cout << "Starting tests..." << endl;

   dummy(argc, argv);

   cout << "Tests are finished." << endl;

   return retval;
}
