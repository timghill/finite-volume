# Shallow water equations
Our finite volume methods were originally developed for the shallow water equations. In a few senses this was an odd choice:

 * Our eventual model has a scalar state variable.
 * Our model is diffusive, so we do not need to add numerical diffusion.

Nevertheless, they are a good test for our methods and provide a good case study.

## Two-dimensional shallow water equations
The shallow water solver solves the two-dimensional system of equations

$$\frac{\partial}{\partial t} \begin{bmatrix} h \\\ hu \\\ hv \end{bmatrix} + \frac{\partial}{\partial x} \begin{bmatrix} hu \\\ hu^2 + \frac{1}{2}gh^2 \\\ huv \end{bmatrix} + \frac{\partial}{\partial y} \begin{bmatrix} hv \\\ huv \\\ hv^2 + \frac{1}{2}gh^2 \end{bmatrix} = 0$$

We note that the conserved quantities in these equations are not the physical variables \\((h, u, v)\\), but are the depth and momenta \\(\vec{u} = (h, m\_x, m\_y)\\), where \\(m\-x = hu\\) and \\(m\_y = hv\\). In terms of these conserved quantities, the equations become

$$\frac{\partial}{\partial t} \begin{bmatrix} h \\\ m\_x \\\ m\_y \end{bmatrix} + \frac{\partial}{\partial x} \begin{bmatrix} m\_x \\\ \frac{m\_x^2}{h} + \frac{1}{2}gh^2 \\\ \frac{m\_x m\_y}{h} \end{bmatrix} + \frac{\partial}{\partial y} \begin{bmatrix} m\_y \\\ \frac{m\_x m\_y}{h} \\\ \frac{m\_y^2}{h} + \frac{1}{2}gh^2 \end{bmatrix} = 0$$

We also note the symmetry between \\(x\\) and \\(y\\). We recognize the \\(\frac{\partial}{\partial x}\\) term as representing the flux in the \\(x\\)-direction, and the \\(\frac{\partial}{\partial y}\\) term as representing the flux in the \\(y\\)-direction (for example by integrating over a square). To apply finite volume methods, we therefore need to write the flux generally for normal and tangent momentum components.

Let \\(m_n\\) by the normal component of momentum, and \\(m_\tau\\) be the tangent component with respect to a a cell boundary. Then, we have that the flux is

$$ \vec{f}(\vec{u}) = \begin{bmatrix} f_h \\\ f_n \\\ f_\tau \end{bmatrix} = \begin{bmatrix} m_n \\\ \frac{m_x^n}{h} + \frac{1}{2}gh^2 \\\ \frac{m_n m_\tau}{h} \end{bmatrix}$$

Note that the first component says that convergence of the normal component of momentum causes an increase in water depth. This is exactly what we want! Now, we can write the system of equations as

$$\frac{\partial \vec{u}}{\partial t} + \frac{\partial}{\partial n} \vec{f}(\vec{u}) + \frac{\partial}{\partial \tau} \vec{f}(\vec{u}) = 0$$

where $\vec{n} and \vec{\tau}$ are the normal and tangent vectors.

## Finite volume methods
As usual, we integrate over a control volume \\(\Omega\\). However, we note that only the normal component of flux will act to exchange mass across an edge. We also note that we will replace the flux \\(\vec{f}\\) with the *numerical* flux \\(\vec{f}^\*(\vec{u}_1, \vec{u}_2)\\), which depends on the state vector on both sides of the boundary. Therefore, we end up with

$$\frac{\partial \vec{u}_i}{\partial t} + \frac{1}{|\Omega_i|} \sum_j \vec{f}^\*(\vec{u}_i, \vec{u}_j)l_j = 0$$

where \\(\vec{u}_i\\) is the *average* of \\(\vec{u}\\) on element \\(\Omega_i\\), \\(|\Omega_i|\\) is the area of the element, and \\(l_j\\) is the length of edge \\(\Gamma_j\\).

### Numerical flux function
We use a simple first-order numerical flux function, where we take the average of values on neighbouring elements. We add the minimum diffusion required for stability. This turns out to be

$$\vec{f}^\*(\vec{u}_1, \vec{u}_2) = \frac{1}{2}[\vec{f}(\vec{u}_1) + \vec{f}(\vec{u}_2)] - \frac{\alpha}{2}(\vec{u}_2 - \vec{u}_1)$$
where \\(\alpha = \frac{|m_n|}{h} + \sqrt{h}\\) (which can be calculated from the eigenvalues of the flux Jacobian).

## Boundary conditions
We apply ideal wall boundary conditions. That is, we implement ghost cells where we reverse the momentum components

$$\vec{u}_j = \begin{bmatrix} h \\\ -m\_{n,i} \\\ -m\_{\tau,i} \end{bmatrix}$$

This ensures no mass crosses the boundary.
