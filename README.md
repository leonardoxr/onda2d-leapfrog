# 2D Wave Propagation Simulator - Leapfrog Method

> **📚 College Project Notice**: This was developed as a college project for studying numerical methods in physics simulations. While functional and educational, it's not production-grade code - so don't take it too seriously! 😊

## Overview

A numerical simulation of 2D wave propagation using the **Leapfrog finite difference method**. This program simulates shallow water waves (like tsunamis or ocean waves) propagating over variable depth terrain.

The simulator solves the 2D wave equation with terrain-dependent propagation, modeling how waves behave when traveling over different depths - similar to how tsunamis amplify as they approach shallow coastal waters.

## What It Does

- **Simulates 2D wave propagation** over a 71x71 grid
- Uses the **leapfrog time-stepping scheme** for numerical stability
- Models **variable terrain depth** (λ) affecting wave speed
- Generates initial Gaussian wave profile ("bell curve" shaped wave)
- Outputs simulation data for visualization

## Mathematical Background

The code implements a discretized version of the 2D wave equation:

```
u_new(x,y) = 2·u_now(x,y) - u_old(x,y) + Δt²·∇²u
```

Where:
- `u(x,y,t)` = wave height at position (x,y) and time t
- `λ(x,y)` = terrain depth function (affects wave speed)
- `rx, ry` = Δt/Δx and Δt/Δy (time/space ratios)

The leapfrog method uses three time levels (past, present, future) for second-order accuracy in time.

## Features

- **Initial Conditions**:
  - Gaussian wave profile centered at origin
  - Inverted Gaussian terrain (simulates a depression/hole in the seafloor)

- **Boundary Conditions**:
  - Custom handling for all edges and corners
  - Reflective boundaries (waves reflect at domain edges)

- **Grid**: 71 × 71 spatial points
- **Time Steps**: 300 iterations
- **Output**: ASCII data file (`onda.txt`) suitable for plotting

## Compilation and Usage

### Requirements
- GCC compiler (or any C compiler)
- Basic Unix/Linux environment (or Windows with MinGW)

### Compile

```bash
gcc leapfrog2d.c -o leapfrog2d -lm
```

The `-lm` flag links the math library (required for `exp()`, `pow()` functions).

### Run

```bash
./leapfrog2d
```

This generates `onda.txt` containing the simulation results.

### Output Format

The output file `onda.txt` contains diagonal slice data:
```
x y u(x,y)
0 0 0.234
1 1 0.456
...
```

Only points where `x == y` are printed (diagonal of the grid), showing wave evolution along this line.

## Visualization

While the code doesn't include built-in visualization, you can plot the results using:

### Using Gnuplot
```bash
gnuplot
> plot 'onda.txt' with lines
```

### Using Python
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('onda.txt')
plt.plot(data[:, 0], data[:, 2])
plt.xlabel('Position')
plt.ylabel('Wave Height')
plt.title('2D Wave Propagation (Diagonal Slice)')
plt.show()
```

## Code Structure

```
leapfrog2d.c
├── main()                 # Main simulation loop
├── gauss()                # Generates Gaussian initial wave
├── H()                    # Terrain depth function
├── atualizar_onda()       # Updates wave state (main algorithm)
└── delta_u()              # Computes spatial derivatives
```

### Key Parameters

Located at top of `leapfrog2d.c`:
- `nx_, ny_`: Grid dimensions (71 × 71)
- `tmax`: Number of time steps (300)
- `rx, ry`: Courant numbers (0.25) - controls numerical stability

**Note**: For numerical stability, the Courant condition requires `rx, ry < 0.5`.

## Physics Explained Simply

Imagine dropping a pebble in a pond:
1. Creates ripples (our initial Gaussian wave)
2. Ripples spread outward (wave propagation)
3. Shallow areas slow the wave (terrain effect via λ)
4. Ripples reflect off pond edges (boundary conditions)

This code simulates exactly that, but in 2D and with more control over terrain shape!

## Known Limitations

⚠️ **Remember: This is a college project!**

- Only saves diagonal slice (not full 2D grid at each timestep)
- Hardcoded parameters (no command-line options)
- Portuguese comments in source code
- No real-time visualization
- Fixed grid size (requires recompilation to change)
- Basic file I/O (could overflow for large grids)

## Educational Value

This project demonstrates:
- ✅ Finite difference methods for PDEs
- ✅ Leapfrog time integration
- ✅ Numerical stability considerations
- ✅ Boundary condition handling
- ✅ C programming for scientific computing
- ✅ Array manipulation and memory management in C

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Since this is an educational/archival project, we're not actively seeking contributions. However, feel free to:
- Fork it for your own experiments
- Use it as a learning resource
- Report any critical bugs (though remember - college project! 😄)

## Acknowledgments

Developed as part of computational physics coursework, exploring numerical methods for solving partial differential equations.

---

**Disclaimer**: This code was written for educational purposes. If you're looking for production-quality wave simulation, check out professional CFD packages like OpenFOAM, Gerris, or commercial options. But if you want to understand the basics of numerical wave propagation - you're in the right place!
