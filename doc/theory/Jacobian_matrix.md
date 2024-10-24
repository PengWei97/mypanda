- **创建时间**: 2024年10月24日
- **创建者**: Peng Wei, GPT4
- **文件目的**: 本文件旨在总结和记录 Jacobian matrix 的相关知识。
---

In **continuum mechanics** and the **finite element method (FEM)**, the **Jacobian matrix** plays a crucial role in solving nonlinear systems of equations and transforming variables between different coordinate systems. The Jacobian matrix can refer to slightly different things depending on the context, but generally, it captures how a set of variables changes with respect to another set of variables, allowing transformations and linearizations of nonlinear problems.

Here’s a breakdown of the role of the Jacobian matrix in these two contexts:

---

## 1. **Jacobian Matrix in Continuum Mechanics**

In **continuum mechanics**, the Jacobian matrix often describes the relationship between the reference (undeformed) and current (deformed) configurations of a material body. Specifically, it captures the **local deformation gradient** in problems involving large deformations. 

### Deformation Gradient
For a body undergoing deformation, the **deformation gradient tensor** $\mathbf{F}$ is a key quantity, which relates the differential change in position in the deformed configuration $\mathbf{x}$ to the differential change in the reference configuration $\mathbf{X}$. Mathematically:

$
\mathbf{F} = \frac{\partial \mathbf{x}}{\partial \mathbf{X}}
$

Here, $\mathbf{x}$ is the position in the deformed configuration, and $\mathbf{X}$ is the position in the reference configuration.

- The Jacobian determinant $J$ in continuum mechanics is often denoted as:

$
J = \text{det}(\mathbf{F})
$

This determinant represents the ratio of the volume change due to the deformation. If $J = 1$, the volume remains unchanged, while $J > 1$ indicates expansion, and $J < 1$ indicates compression.

### Importance in Continuum Mechanics:
- The Jacobian matrix ($\mathbf{F}$) is central to the formulation of **strain tensors** (such as Green-Lagrange strain) and the computation of **stress tensors** in nonlinear elasticity.
- It enables the transformation of integrals between reference and deformed configurations, which is key in computing physical quantities like strain energy, stress, or forces during large deformations.

---

## 2. **Jacobian Matrix in the Finite Element Method (FEM)**

In the **finite element method**, the Jacobian matrix is used to relate the **local (element) coordinates** to the **global (system) coordinates** during the process of numerical integration (such as Gauss quadrature) and element formulation.

### Mapping Between Coordinate Systems
When solving finite element problems, an element is typically described in **local coordinates** (natural coordinates), but the equations must be evaluated in the **global coordinates** of the entire system. The Jacobian matrix in FEM allows for this transformation. It maps the differential changes between the natural coordinates $(\xi, \eta, \zeta)$ and the global Cartesian coordinates $(x, y, z)$.

For a given element, the position of a point in global coordinates $\mathbf{x}$ is expressed in terms of the shape functions $N_i(\xi, \eta, \zeta)$ and the nodal coordinates $\mathbf{x}_i$:

$
\mathbf{x} = \sum_{i=1}^n N_i(\xi, \eta, \zeta) \mathbf{x}_i
$

The **Jacobian matrix** $\mathbf{J}$ in FEM relates the differentials of the global coordinates to the local coordinates as:

$
\mathbf{J} = \frac{\partial \mathbf{x}}{\partial \boldsymbol{\xi}} = \begin{bmatrix}
\frac{\partial x}{\partial \xi} & \frac{\partial x}{\partial \eta} & \frac{\partial x}{\partial \zeta} \\
\frac{\partial y}{\partial \xi} & \frac{\partial y}{\partial \eta} & \frac{\partial y}{\partial \zeta} \\
\frac{\partial z}{\partial \xi} & \frac{\partial z}{\partial \eta} & \frac{\partial z}{\partial \zeta}
\end{bmatrix}
$

where $\xi, \eta, \zeta$ are the local (natural) coordinates of the element.

### Importance in FEM:
- **Transformation of variables**: The Jacobian matrix is necessary for transforming integrals (such as strain energy, internal forces, or stiffness) from the local coordinates (natural coordinates) to global coordinates.
- **Numerical integration**: During Gauss quadrature, the Jacobian matrix is used to compute the volume (or area) element in the global coordinate system. The integration over the element is carried out in natural coordinates, but the Jacobian accounts for how the element is mapped to the global system.
  
- **Deformation in non-linear FEM**: In problems involving large deformations, the Jacobian matrix is key to updating the element geometry and computing the deformed state.

---

### Summary

- **In continuum mechanics**, the **Jacobian matrix** (or the deformation gradient $\mathbf{F}$) describes how the material deforms from the reference configuration to the current configuration. It is crucial for calculating quantities like strains, stresses, and volume changes.
  
- **In the finite element method**, the **Jacobian matrix** is used to transform between the local element coordinates (natural coordinates) and global coordinates. It plays a central role in numerical integration, mesh deformation, and solution of nonlinear problems.

Both in continuum mechanics and FEM, the Jacobian matrix is a core mathematical tool that allows handling nonlinearities, geometric transformations, and integration in complex, multidimensional problems.