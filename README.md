# ASCE Sediment Transport Equations in Python

A comprehensive, open-source digital library of the equations presented in the American Society of Civil Engineers (ASCE) Manual on Sedimentation Engineering.

## Project Goal

The primary objective of this project is to translate the extensive collection of sediment transport equations from the ASCE Manual into a clear, usable, and open-source library of Python functions. This aims to make these critical formulas more accessible for researchers, engineers, and students.

## Methodology

The initial codebase was generated through an AI-assisted process. An AI model was used to parse the PDF of the ASCE manual, identify mathematical equations, and translate them into Python code. While this approach allows for rapid development, it is susceptible to errors in interpretation, transcription, and context.

## ⚠️ Important Disclaimer: Project in Alpha Stage

This project is currently in an **unverified, alpha stage**. The Python functions have been automatically generated and **have not yet been systematically verified** against the source material.

**It is almost certain that errors exist in the code.** These may include:
* Incorrect mathematical operations
* Mismatched variables or constants
* Missing context or applicability limits

> **DO NOT USE ANY EQUATION FROM THIS REPOSITORY FOR PROFESSIONAL OR ACADEMIC WORK WITHOUT FIRST PERSONALLY VERIFYING ITS ACCURACY AGAINST THE OFFICIAL ASCE MANUAL.** The maintainers of this repository assume no liability for any consequences arising from the use of unverified code.

## How to Contribute

Community involvement is essential to achieving the goal of a fully verified and reliable library. We strongly encourage users to contribute by correcting errors.

If you find an equation that is incorrect, please help us fix it by following these steps:
1.  **Fork** the repository.
2.  Create a new branch for your changes (e.g., `fix/equation-5-12`).
3.  Correct the equation in the relevant Python file. Please add comments referencing the page and equation number from the manual.
4.  **Submit a Pull Request** to the `main` branch.
5.  In your pull request description, please clearly explain the error you found and cite the page number in the manual that contains the correct formula.

Your contributions will help make this a valuable resource for the entire hydraulic engineering community.

## Roadmap
* [ ] **Phase 1:** Initial AI-driven generation of all equations. (*In Progress*)
* [ ] **Phase 2:** Community-driven verification and correction of all equations. (*In Progress*)
* [ ] **Phase 3:** Development of unit tests for verified equations to ensure long-term stability.
* [ ] **Phase 4:** Refactoring the code into a more structured Python package.