
# Stable Neo-Hookean Flesh Simulation
[Stable Neo-Hookean](https://graphics.pixar.com/library/StableElasticity/paper.pdf) is a hyperelastic energy that remain stable under large deformations. This is the FEM implementation with Stable Neo-Hookean energy in [libigl](https://libigl.github.io)

## Finite Element Method
 In order to enable the physics simulation, all deformation laws have to be discretized. Tetrahedral meshes are common discrete volumetric geometry representations. All the quantity inside the volume can be represent use a shape function $\phi_i(x)$ where $x \in R^3$ is a position inside the space volume.

## Deformations
In order to determine forces due to deformation, we need to know how the nearby points have moved relative to one another. This information is captured by the **Deformation Gradient** $F$.
## The Strain Energy
$$
\psi(F) = \frac{\mu}{2}(I_C - 3) + \frac{\lambda}{2}(J - 1 - \frac{\mu}{\lambda})^2
$$
where $J=\text{det}(F)$, $I_C = \text{tr}(F^TF)$, $\mu$ and $\lambda$ are Lamé constants.


## Dependnecy 
```
sudo apt-get install git
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install libx11-dev
sudo apt-get install mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev
sudo apt-get install libxrandr-dev
sudo apt-get install libxi-dev
sudo apt-get install libxmu-dev
sudo apt-get install libblas-dev
sudo apt install libxinerama-dev libxcursor-dev
```

## How to run?
Build:
```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```
Run with Stable Neo-Hookean
```
./stable_neohooken stable
```
Run without Stable Neo-Hookean
```
./stable_neohooken nonstable
```

## Notice
The framework of this project that I used is based on [assignment](https://github.com/dilevin/CSC417-a3-finite-elements-3d) from Prof. David I.W. Levin.