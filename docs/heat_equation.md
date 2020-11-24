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

The difficulty lies in calculating $\nabla u$.
