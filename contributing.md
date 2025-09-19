# Contributing to OptSurvCutR

First off, thank you for considering contributing to `OptSurvCutR`. It's people like you that make the R open-source community such a fantastic place.

We welcome any and all contributions, from bug reports and feature suggestions to code contributions.

## How Can I Contribute?

### Reporting Bugs

If you find a bug, please open an issue on our [GitHub Issues page](https://github.com/paytonyau/OptSurvCutR/issues).

When reporting a bug, please include:
* Your operating system and R session information (`sessionInfo()`).
* A reproducible example of the bug. The `reprex` package is a great way to create one.
* Any error messages you received.

### Suggesting Enhancements

If you have an idea for a new feature or an improvement to an existing one, please open an issue on our [GitHub Issues page](https://github.com/paytonyau/OptSurvCutR/issues). Please provide a clear and detailed description of the feature you are proposing and why it would be valuable.

### Code Contributions via Pull Requests

We welcome code contributions! If you would like to contribute code, please follow these steps:

1.  Fork the repository.
2.  Create a new branch for your feature (`git checkout -b feature/AmazingFeature`).
3.  Make your changes. Please follow the existing code style.
4.  If you are adding new functionality, please add tests for it in the `tests/testthat/` directory.
5.  Ensure that `devtools::check()` runs without any errors, warnings, or notes.
6.  Open a Pull Request.

## Code Style

Please try to follow the general coding style used throughout the package. We generally follow the [tidyverse style guide](https://style.tidyverse.org/).