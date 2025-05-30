# Types

## Naming

- Please name types (classes, structures) as ***descriptive nouns***,
  that describe what the type *is* or what its *purpose* is.
    - As an example, `AtmosphericManipulator` and not `Manipulator` or `AtMnp`.
    - If you find yourself unable to name your type this way, that may
      indicate that it should be broken into independent classes, variables,
      or functions.
- Types should make use of `CamelCase` and begin capitalized
  (e.g., not `camelCase`).
- There should be a logical grouping among the variables and classes contained
  in a type, such that the name captures that association.
    - To demonstrate, consider the following `Fruit` class.

        ```c++
        class Fruit {
          enum FruitName {apple, banana, orange};
        
          FruitName m_fruit;
        
          void is_juiceable(FruitName fruit_) { ... };
        
          enum Color {red, green, blue};

          Color fruit_color;
        
          bool classify_good_color(Color color_) { ... };
        }
        ```

    - It would be better to break this into a separate `Fruit` and `Color`
      class because:
        1. `Color` is not inherently associated with fruits--it could also
           belong to, for instance, a `Vegetable` class.
        2. There is no strong reason the `classify_good_color()` needs to be a
           part of `Fruit`, since it would work the same way if it knew nothing
           about fruits.
