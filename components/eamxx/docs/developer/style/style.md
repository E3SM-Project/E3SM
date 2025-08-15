# EAMxx Code Style Standards

EAMxx does not currently impose strict styling standards, other than those of
the autoformatter.
However, if we can follow consistent style and usage guidelines, we can make
EAMxx developers' lives easier and turnaround times quicker.
In the age of modern text editors with autocomplete and running our code on
machines with enormous storage capacity, there is no real need to save screen
or disk space with minimal naming strategies or other meaning-obscuring style
choices.
So, please adhere to the spirit of these guidelines, and your fellow developers
will thank you for it!

## General Guidelines

- Please give ***descriptive***, rather than ***terse***, names to your
  variables, functions, classes, etc.
- Avoid idiosyncratic naming or styling conventions.
      - If the utility of a style decision would not be immediately apparent to
        a reasonable developer, please avoid its usage.
- In the hierarchy of coding concerns, correctness and speed take the top
  tiers, but below that level, please favor readability and clarity over
  space savings, line-breaking, heavy use of arcane language features, etc.
      - That said, if an arcane language feature ***is*** the correct tool to
        be used, please do so.
          - And perhaps add a comment to aid the non-Jedi-Master developers who
            read the code next. :slightly_smiling_face:
- With regard to comments, the general wisdom is that the correct amount of
  comments is the amount required to make the code clear and easy to understand.
      - To quote
        [Jeff Atwood](https://blog.codinghorror.com/code-tells-you-how-comments-tell-you-why/)
        from his essay about code comments:
      > Only at the point where the code *cannot* be made easier to understand
      should you begin to add comments.
      - However, since that does not offer much in the way of prescriptive
        guidance, here are some guidelines gathered from around the internet.
          - A general rule of thumb, though not always true, is that comments
            should explain ***why*** something is being done in the code,
            rather than ***what*** the code is doing.[^butwhattabout]
          - Comments should dispel confusion and not cause it.
          - Comments should not duplicate the code.[^dupe]
          - Provide links or references when using code from elsewhere or
            whenever it is otherwise appropriate.
          - A bad comment is worse than no comment.
          - Add comments to explain non-standard usage or unidiomatic code.

[^butwhattabout]: An obvious exception to this is explaining complex or opaque
parts of the code that cannot be made simpler--for instance, a clever arithmetic
trick in a complicated interpolation scheme.
[^dupe]: For example, this type of comment does not add any information and
is unnecessary.

    ```c++
    // perform initialization
    this->initialize();
    ```
