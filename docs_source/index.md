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

Now, let \\(\bar{u}\\) be the average value of \\(u\\) over \\(\Omega\\), \\(F\\) be the average of \\(f\\), and \\(|\Omega|\\) be the area. Then we have

$$ |\Omega|\frac{\partial \bar{u}}{\partial t} + \int\_{\partial\Omega}\vec{q}(u, \nabla u)\cdot\vec{n}dl = |\Omega| F$$

This yields an important equation for the average value of \\(u\\),

$$ \frac{\partial \bar{u}}{\partial t} + \frac{1}{|\Omega|}\int\_{\partial\Omega}\vec{q}(u, \nabla u)\cdot\vec{n}dl = F$$

This equation says that the average value of \\(u\\) only changes due to the net flux flowing out through the boundary (first term) or the source term \\(F\\). This result should not be surprising or new, but it is an important step along the way to the semi-discrete form.

To discretize in space, suppose that \\(\Omega\\) is a polygon with \\(N\\) edges (extend to 3D by considering a polyhedron with \\(N\\) faces). Then we can break up the contour integral to be the sum of contributions from each edge. On each edge, we replace the line integral by the average value \\(\vec{q}_k\\) multiplied by the edge length \\(l_k\\). Therefore, we have

$$ \frac{\partial \bar{u}}{\partial t} + \frac{1}{|\Omega|} \sum_{k=1}^N \vec{q}_k \cdot \vec{n}_k l_k = F$$

This is called the semi-discrete form, since the spatial component has been discretized but the temporal component has not. Note that the method is not yet complete since we have not specified how to calculate the edge values \\(\vec{q}_k\\).
