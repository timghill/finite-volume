# Shallow water equations
The shallow water solver solves the two-dimensional system of equations

$$\frac{\partial}{\partial t} \begin{bmatrix} h \\\ hu \\\ hv \end{bmatrix} + \frac{\partial}{\partial x} \begin{bmatrix} hu \\\ hu^2 + \frac{1}{2}gh^2 \\\ huv \end{bmatrix} + \frac{\partial}{\partial y} \begin{bmatrix} hv \\\ huv \\\ hv^2 + \frac{1}{2}gh^2 \end{bmatrix} = 0$$

We note that the conserved quantities in these equations are not the physical variables \\((h, u, v)\\), but are the depth and momenta \\(\vec{u} = (h, m\_x, m\_y)\\), where \\(m\-x = hu\\) and \\(m\_y = hv\\). In terms of these conserved quantities, the equations become

$$\frac{\partial}{\partial t} \begin{bmatrix} h \\\ m\_x \\\ m\_y \end{bmatrix} + \frac{\partial}{\partial x} \begin{bmatrix} m\_x \\\ \frac{m\_x^2}{h} + \frac{1}{2}gh^2 \\\ \frac{m\_x m\_y}{h} \end{bmatrix} + \frac{\partial}{\partial y} \begin{bmatrix} m\_y \\\ \frac{m\_x m\_y}{h} \\\ \frac{m\_y^2}{h} + \frac{1}{2}gh^2 \end{bmatrix} = 0$$

We also note the symmetry between \\(x\\) and \\(y\\). We recognize the \\(\frac{\partial}{\partial x}\\) term as representing the flux in the \\(x\\)-direction, and the \\(\frac{\partial}{\partial y}\\) term as representing the flux in the \\(y\\)-direction (for example by integrating over a square). To apply finite volume methods, we therefore need to write the flux generally for normal and tangent momentum components.

Let \\(m_n\\) by the normal component of momentum, and \\(m_t\\) be the tangent component with respect to a a cell boundary. Then, we have that the flux is

$$ \vec{f}(\vec{u}) = \begin{bmatrix} f_h \\\ f_n \\\ f_\Tau \end{bmatrix} = \begin{bmatrix} m_n \\\ \frac{m_x^n}{h} + \frac{1}{2}gh^2 \\\ \frac{m_n m_t}{h} \end{bmatrix}$$

Note that the first component says that convergence of the normal component of momentum causes an increase in water depth. This is exactly what we want! Now, we can write the system of equations as

$$\frac{\partial \vec{u}}{\partial t} + \frac{\partial}{\partial n} \vec{f}(\vec{u}) + \frac{\partial}{\partial \Tau} \vec{f}(\vec{u}) = 0$$.
