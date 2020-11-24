# Finite volume methods for conservation laws
This documentation describes the implementation of first-order finite volume methods on unstructured triangular meshes to solve the shallow water, heat, and nonlinear heat equations.

## General finite volume methods
Before we solve specific examples, we should review finite volume methods. The methods here are written in a general way so that we can apply them to unstructured triangular meshes. The eventual form can be used for any shape of mesh, including regular structured meshes or unstructured mixed meshes.

Consider a scalar conservation law of the form

$$ \frac{\partial u}{\partial t} + \nabla \cdot \vec{q}(u, \nabla u) = f$$
