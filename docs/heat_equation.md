# Heat equation
The heat equation is usually written

$$\frac{\partial u}{\partial t} = \gamma \nabla^2 u,$$

but it can also be written as a conservation law,
$$\frac{\partial u}{\partial t} + \nabla \cdot \vec{q} = f$$,
with
$$\vec{q} = -\gamma\nabla u$$.

Therefore, we get the same semi-discrete form,

$$ \frac{\partial \bar{u}}{\partial t} + \frac{1}{|\Omega|} \sum_{k=1}^N \vec{q}_k \cdot \vec{n}_k l_k = F.$$

The primary difficulty in solving the heat equation is that the flux depends directly on the gradient of \\(u\\). The numerical flux function we use is simply the exact flux,

$$q_k = -\gamma \nabla u\big\rvert_{\Gamma_k}.$$

The difficulty lies in calculating $\nabla u$. We have a few different approaches of varying complexity.

## Green-Gauss method
The Green-Gauss method uses the Green-Gauss/divergence theorem to calculate the gradient. We follow the same approach as we did to derive the general finite volume semi-discrete form, and find

$$\nabla u = \frac{1}{|\Omega|} \sum_{k=1}^N u_k \vec{n}_k$$,
where the edge value \\(\vec{u}_k\\) is calculated as the average of the neighbouring elements,
$$\vec{u}_k = \frac{1}{2}(\vec{u}_1 + \vec{u}_2)$$.

This method is simple to implement, but it has been shown to be *inconsistent* on unstructured meshes. That is, the gradient calculation converges to an incorrect value as the mesh is refined.

## Green-Gauss hybrid method
We can fix the issues with the Green-Gauss method by being smarter in calculating the edge value \\(u_k\\). Let \\(\{\Gamma_i\\}_{i=1}^m\\) be the set of elements that share a node with edge \\(\Gamma_k\\). Then, for an arbitrary element, we can expand

$$u_i = u_k + \frac{\partial u_k}{\partial x}\Delta x + \frac{\partial u_k}{\partial y}\Delta y.$$

Combining the \\(m\\) equations, we find a matrix system

$$\begin{bmatrix} 1 && \Delta x_1 && \Delta y_1 \\\ 1 && \Delta x_2 && \Delta y_2 \\\ && \vdots && \\\ 1 && \Delta x_m && \Delta y_m \end{bmatrix} \begin{bmatrix} u_k \\\ \frac{\partial u_k}{\partial x} \\\ \frac{\partial u_k}{\partial y} \end{bmatrix} = \begin{bmatrix} u_1 \\\ u_2 \\\ vdots \\\ u_m \end{bmatrix}$$
