#include <iostream>
#include <string_view>

#include <dexpr/ast.hpp>
#include <dexpr/lexer.hpp>
#include <dexpr/parser.hpp>
#include <dexpr/supported_functions.hpp>

namespace {

void print_functions() {
    // Builtins only; a model's own vocabulary is a superset of this.
    std::cout << "Builtin functions\n\n";

    for (const auto& entry : dexpr::builtin_functions()) {
        // Each spec spans a few lines, so separate entries with a blank one.
        std::cout << "  " << entry.second << "\n\n";
    }
}

// Checks grammar plus the builtins. A component's own functions are unknown
// here; it checks those by calling validate_calls() with its own registry.
int check_expression(std::string_view input) {
    dexpr::ast::ExprPtr expr;
    try {
        dexpr::parser::Parser parser{dexpr::Lexer{std::string{input}}};
        expr = parser.parse();
    } catch (const std::exception& e) {
        std::cerr << "does not parse:\n" << e.what() << '\n';
        return 1;
    }

    std::cout << "parsed as: " << dexpr::ast::to_string(*expr) << '\n';

    try {
        dexpr::validate_calls(*expr, dexpr::builtin_functions());
    } catch (const dexpr::ValidationError& e) {
        std::cerr << e.what();
        return 1;
    }

    std::cout << "ok\n";
    return 0;
}

// The builtins against their own examples; same check a component runs.
int check_registry() {
    try {
        dexpr::validate_registry(dexpr::builtin_functions());
    } catch (const dexpr::ValidationError& e) {
        std::cerr << e.what();
        return 1;
    }
    std::cout << "builtin functions: "
              << dexpr::builtin_functions().names().size() << " ok\n";
    return 0;
}

void print_help() {
    std::cout <<
R"(Usage:
    dexpr functions          list the builtin functions
    dexpr check <expr>       parse <expr> and check its calls
    dexpr check-registry     check every builtin function against its example
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

    if (command == "check") {
        if (argc < 3) {
            std::cerr << "check needs an expression\n";
            return 1;
        }
        return check_expression(argv[2]);
    }

    if (command == "check-registry") {
        return check_registry();
    }

    if (command == "help") {
        print_help();
        return 0;
    }

    std::cerr << "Unknown command: " << command << '\n';
    return 1;
}
