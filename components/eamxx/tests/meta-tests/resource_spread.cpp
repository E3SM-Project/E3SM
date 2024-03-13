#include <chrono>
#include <thread>

int main(int, char**)
{
  // Nothing needs to be done here except sleeping to give a chance
  // for tests to run concurrently.
  const auto seconds_to_sleep = 5;
  std::this_thread::sleep_for(std::chrono::seconds(seconds_to_sleep));
  return 0;
}
