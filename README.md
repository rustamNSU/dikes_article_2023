# MAGMA DIKES

## Naming conventions
- **camelCase** for variables
- **UpperCase** for class and structures
- **under_scores** for function, e.g. solve_lubrication(..)

## Custom structures
### Mesh
- **xc** - coordinates of the center of the mesh elements
- **leftBoundary** - coordinates of the mesh elements left boundary
- **rightBoundary** - coordinates of the mesh elements right boundary
- **dx** - mesh elements width
- **n** - number of mesh elements

### Reservoir
- **E** - Young's modulus
- **nu** - Poisson's coefficient
- **toughness** - kIC
- **leakoff** - Carter's leakoff coefficient
- **Ep**
- **Kp**
- **Cp**
- **mu** - viscosity of fluid (magma)
- **mup**

## Filepaths and folders
Name of file path or folder should be relative (for compatibility in different platforms and devices). For this purpose use next commands:
- pwd
- fullfile
- genpath

## Comments in matlab code
- The end lines in nested loops can have identifying comments