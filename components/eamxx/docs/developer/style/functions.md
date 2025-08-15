# Functions and Methods

## Naming

- Please name functions (methods) to be ***descriptive*** and ***verb**-ish* so
  that the name describes what the function *does*.
      - For example, `install_flux_capacitor()` instead of `f_capacitor_method()`
        or `ifc()`.
      - To note, if your function cannot be named verb-ishly, that is probably a
        sign that there is a fundamental issue with your function.
- In general, functions should use `snake_case` and not contain capitalization,
  unless capitalizing that quantity is a standard convention (e.g.,
  `find_Doc_Brown()`).

## Usage

- In situations where multiple class-member variables are passed as function
  arguments, favor passing the object, rather than the individual variables,
  for the sake of shorter lines and function signatures.
- As a rule of thumb, developers should prefer passing arguments by reference,
  rather than by value, especially in cases for which the argument is ***large***
  or structurally ***complex***.
      - E.g., prefer `func(Type &x) { ... }` versus `func(Type x) { ... }`.
      - This holds true much of the time because passing by value creates a copy
        of the argument for use by the function, and this increases memory
        pressure and a resultant performance penalty.
- To be more in-depth on the topic, we echo the guidance from the
  [***C++ Core Guidelines***](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#S-introduction)
  and would encourage the curious to read more deeply on the topic.

      > - For "in" parameters, pass cheaply-copied types by value and others by
      > reference to `const`.
      > - For "in-out" parameters, pass by reference to non-`const`.
      > - For "out" parameters, prefer return values to output parameters.

      - The authors go on to explain that "cheap to copy" includes variables
        holding 2 or 3 words--for example a double, pointer, or reference.
      - To illustrate, the below function adheres to the provided guidance.

      ```c++
      SmallStruct func(LargeNComplex &notSmall, double forty_two,
                       const int *laser, const BigStruct &lotsaData) {

        // SmallStruct object is only assigned-to and then returned, so it is
        // defined as the return type
        SmallStruct ans;
        // forty_two, a scalar double, is considered small and so passed by value
        ans.val = forty_two;
        // lotsaData, a BigStruct object, is only accessed and so is
        // passed by const reference
        ans.other_val = lotsaData.small_piece;

        // we pass laser by value and as const because pointers are cheap to copy,
        // and the value is not changed
        for (int i = 0; i < forty_two; ++i) {
          ans.val_vector.push_back(laser[i] + notSmall.epsilon[i]);
        }

        // notSmall is large and also modified in the function, so we pass by
        // reference
        // it is also accessed, above, so it must be an argument and not the
        // return value
        notSmall.ans_struct = ans;
        return ans;
      }
      ```
