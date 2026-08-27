#include <iostream>
#include <string_view>

#include <dexpr/supported_functions.hpp>

namespace {

void print_functions() {
    std::cout << "Supported functions\n\n";

    for (const auto& function : dexpr::supported) {
        std::cout << "  " << function << '\n';
    }
}

void print_help() {
    std::cout <<
R"(Usage:
    dexpr functions
    dexpr help
)";
}

} // namespace

int main(int argc, char* argv[]) {
    if (argc == 1) {
        print_help();
        return 0;
    }

    std::string_view command{argv[1]};

    if (command == "functions") {
        print_functions();
        return 0;
    }

    if (command == "help") {
        print_help();
        return 0;
    }

    std::cerr << "Unknown command: " << command << '\n';
    return 1;
}
