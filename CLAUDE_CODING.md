This file contains guidance for Claude when writing R and C++ code and other code for R packages.

## Setup

Before working on a task, read all files in `R/` and `src/` that are relevant
to the task, plus NAMESPACE and DESCRIPTION. When uncertain about scope, read
broadly rather than narrowly --- it is better to understand the surrounding code
than to make changes in isolation. Also read any vignettes that touch the same
functionality. The package uses roxygen2 documentation.

## Code style

First, prefer boring code over clever code for readability and maintainability.
This does not mean prefer for() loops over other approaches --- since
vectorization is both much faster and a core piece of R language use. In fact,
often avoiding a for() loop can make code more clear than a for() loop because
of all of the overhead necessary to setup the objects to modify, etc.

Second, comment the code explaining **why** the particular sections and lines of
code are there more than **what** they are doing.

## File organization

Group functions by conceptual purpose, one file per coherent unit. Do not let
files grow past roughly 300 lines. When a new function is conceptually distinct
from the existing contents of a file (different mechanism, different external
dependency, different layer of abstraction), create a new file rather than
appending.

## Testing

Write well commented unit tests that will go into tests/testthat **before**
refactoring. Tests should represent both that the code runs but most
importantly they should represent the statistical principles underlying and
justifying the code. If I am writing code to square numbers then I should have
tests of squaring numbers --- these are more important tests than tests that
when I provide a numeric input I get a numeric output. Tests should pass before
we can judge the code to be correct.

It is not appropriate to remove failing tests rather than fixing source code or
asking for clarification from me. You can skip tests when you can't quickly
resolve failures, but you will need to remember these as future tasks.

Always prioritize readable, maintainable tests over comprehensive coverage.

## Build discipline

Run `devtools::document()` after adding or modifying roxygen2 documentation.
Run `devtools::check()` before considering a task complete. When adding new
exported functions, bump the patch version in DESCRIPTION (e.g., 0.0.3.0 →
0.0.3.1).

## Workflow with me

Pause for my review at these checkpoints:
1. After writing tests and before writing implementation.
2. After writing implementation and before running `devtools::check()`.
3. Whenever a design decision arises that the plan does not already resolve.

Do not proceed past a checkpoint without my input.

