# Contributing to onda2d-leapfrog

## About This Project

This repository contains a **college project** from a computational physics course. It was created for educational purposes to demonstrate the implementation of the leapfrog method for 2D wave propagation.

## Project Status: Educational/Archival

⚠️ **This is not an active development project.** It exists primarily as:
- A learning resource for students studying numerical methods
- An archive of coursework
- A reference implementation of the leapfrog finite difference method

## Can I Still Contribute?

While we're not actively developing new features, you're welcome to:

### ✅ Encouraged
- **Fork the repository** for your own experiments
- **Use it as a learning resource** for your own studies
- **Submit bug reports** if you find critical errors in the implementation
- **Share improvements** you've made in your fork (we might merge them!)
- **Ask questions** via issues about how the code works
- **Translate comments** from Portuguese to English
- **Add visualization scripts** (Python, Gnuplot, etc.)

### ⚠️ Limited Acceptance
We may accept pull requests for:
- Bug fixes that correct the numerical method
- Documentation improvements
- Code comments/clarity enhancements
- Additional output formats
- Build system improvements (Makefile, CMake)

### ❌ Unlikely to Accept
- Major refactoring or modernization
- Adding complex features
- Changing the core algorithm
- Making it "production-ready"

## Why These Limitations?

This project serves as a **snapshot of coursework**. Dramatically changing it would reduce its value as:
1. A historical record of learning
2. A simple, understandable reference implementation
3. An example of "real student code" for other learners

If you want to build something more advanced, we encourage you to fork and create your own evolution of this project!

## How to Contribute

If you'd like to submit a contribution:

1. **Fork the repository**
2. **Create a feature branch**
   ```bash
   git checkout -b fix/your-fix-name
   ```
3. **Make your changes**
   - Keep changes minimal and focused
   - Maintain the existing code style
   - Add comments if adding new functionality
4. **Test your changes**
   ```bash
   gcc leapfrog2d.c -o leapfrog2d -lm
   ./leapfrog2d
   ```
5. **Submit a pull request**
   - Describe what you changed and why
   - Reference any related issues

## Code of Conduct

Be kind and respectful. Remember:
- This was a **learning project** - the code isn't perfect, and that's okay!
- Questions are always welcome
- Constructive criticism is appreciated
- We're all here to learn

## Questions?

Open an issue! We're happy to discuss:
- How the code works
- The numerical method
- Physics of wave propagation
- Improving your own fork

## Alternative Projects

If you're looking for more advanced wave simulation tools, check out:
- **OpenFOAM**: Professional CFD software
- **Gerris**: Free CFD solver
- **Basilisk**: Adaptive mesh refinement for fluids
- **SWE-1D/2D**: Shallow water equation solvers

## License

All contributions will be under the same MIT License as the main project.

---

Thank you for your interest in this educational project! 🌊
