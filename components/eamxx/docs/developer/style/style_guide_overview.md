# EAMxx C++ Style Guide

EAMxx enforces some standards on ***style*** and ***format***.
For the purpose of this guide, we draw a distinction between these two related
topics and loosely define them as:

- **Style**
      - Style is concerned with how the code ***reads*** on a lexical level as
        well as with the structural and technical decisions that determine how
        the code ***runs***.
      - For example:
          - Descriptive, as opposed to terse, variable names.
          - Employing small, single-purpose functions.
          - Usage or non-usage of advanced C++ features.
          - The usage and content of comments.
- **Format**
      - Format is the domain of how the code ***looks*** or how the code is
        organized.
      - The good news for formatting is that we enforce a strict, LLVM-based
        formatting standard, and we employ
        [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html)
        to automatically conduct and enforce this standard.

More detail regarding each of these topics is contained in the following
sections.
