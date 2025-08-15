# Types, Classes, Structures

## Naming

- Please name types (classes, structures/`struct`, enumerations/`enum`) as
  ***descriptive nouns***,
  that describe what the type *is* or what its *purpose* is.
      - As an example, `AtmosphericManipulator` and not `Manipulator` or `AtMnp`.
      - If you find yourself unable to name your type this way, that may
        indicate that it should be broken into independent classes, variables,
        or functions.
            - See FIXME: below for an example.
- Types, classes, or structs should make use of `UpperCamelCase`.
      - Variables containing an object or an instance of a type, class, or struct,
        should follow the
        [naming convention for variables](variables.md#naming)
        and be `lower_snake_case`.

## Usage

- There should be a logical grouping among the variables and classes contained
  in a type, such that the name captures that association.
      - To demonstrate, consider the following `Fruit` class.

          ```c++
          class Fruit {
            enum FruitName {Apple, Banana, Orange};

            FruitName m_fruit_name;

            void is_juiceable(FruitName fruit_) { ... };

            enum Color {Red, Green, Blue};

            Color fruit_color;

            bool is_favorite_color(Color color_) { ... };
          };
          ```

      - It would be better to break this into a separate `Fruit` and `Color`
        class because:
            1. `Color` is not inherently associated with fruits--it could also
               belong to, for instance, a `Vegetable` class.
            2. There is no strong reason the `is_favorite_color()` needs to be a
               part of `Fruit`, since it would work the same way if it knew nothing
               about fruits.
- Related to the previous guideline, types/classes/structs should ideally
  encapsulate ***concepts*** to improve both readability and usability.
      - This implies that even a small struct has a reason to exist when it
        serves the goal of shrinking or simplifying the outer class.
      - An example of this in EAMxx is the
        [`TimeInterval`](https://github.com/E3SM-Project/E3SM/blob/75b5b0a0c9078e18736860b2445a8975d7de750d/components/eamxx/src/share/util/eamxx_time_stamp.hpp#L114)
        struct that only contains 4 class variables and 3 relatively simple
        methods; yet, it serves to keep the
        [`DataInterpolation`](https://github.com/E3SM-Project/E3SM/blob/75b5b0a0c9078e18736860b2445a8975d7de750d/components/eamxx/src/share/util/eamxx_data_interpolation.hpp#L14)
        class smaller and better-organized.

### Organization

For the sake of consistency to aid readability and findability, we recommend
classes be organized (ordered) according to the following rules.

- Public interfaces should appear at the top of a class declaration.
      - These are the typically what other developers or users will be looking
        for when reading or interacting with your code.
- Place methods first and class variables last, according to the same logic as
  above.
- Group methods of the same "kind" or "purpose."
      - E.g., group *getters* together and separate them from state-changing
        methods that are also grouped together.
      - Avoid interleaving functions of different kind.
- Seek to group class-variable declarations together in a meaningful way.
      - This could be according to the actual C++ ***type***, but even better
        would be done according to how the variables will be used or what
        other parts of the code they will interact with.
      - To illustrate, in the snippet below, a class storing 2 grid-types and
        a `bool` should prefer `class A` to `class B`

          ```c++
          class A {
            MyGrid source;
            MyGrid target;
            bool forward;
          };
          class B {
            MyGrid source;
            bool forward;
            MyGrid target;
          };
          ```

We acknowledge that these rules do not form a complete logical schema, and so,
we defer to the judgement of the developer and their willingness to consider
the poor souls that will one day read their code. :slightly_smiling_face:
