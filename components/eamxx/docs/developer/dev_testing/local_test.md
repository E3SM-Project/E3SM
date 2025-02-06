# Testing Locally

==Add this both to "Full Model" and "Standalone"==

E.g., "Configuring your local machine"

## E.g., Testing on Non-supported Machines via Local Configuration Files

- `${HOME}/.cime` directory:

!!! Note
    Don't run `ne1024` on your laptop!


## Order-of-operations for `test-all-scream`

#### Placeholder
``` mermaid
graph LR
  A[Start] --> B{Error?};
  B -->|Yes| C[Hmm...];
  C --> D[Debug];
  D --> B;
  B ---->|No| E[Yay!];
```
