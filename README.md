# Finite volume methods for conservation laws
This documentation describes the implementation of first-order finite volume methods on unstructured triangular meshes to solve the shallow water, heat, and nonlinear heat equations.

## General finite volume methods
Before we solve specific examples, we should review finite volume methods. The methods here are written in a general way so that we can apply them to unstructured triangular meshes. The eventual form can be used for any shape of mesh, including regular structured meshes or unstructured mixed meshes.

Consider a scalar conservation law for some quantity \\(u\\) of the form

$$ \frac{\partial u}{\partial t} + \nabla \cdot \vec{q}(u, \nabla u) = f$$

Now we integrate over a control area \\(\Omega\\) (this could easily be extended to a volume in 3D space),


$$ \\iint_\Omega \frac{\partial u}{\partial t} dA + \\iint_\Omega\nabla \cdot \vec{q}(u, \nabla u)dA = \\iint_\Omega f dA$$

Apply the divergence theorem to the second term, and pull the time derivative out of the first integral:

$$ \frac{\partial}{\partial t} \\iint_\Omega u dA + \int\_{\partial\Omega}\vec{q}(u, \nabla u)\cdot\vec{n}dl = \\iint_\Omega f dA$$
