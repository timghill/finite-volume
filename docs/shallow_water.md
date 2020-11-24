# Shallow water equations
The shallow water solver solves the two-dimensional system of equations

\begin{align}
x &= \nabla^2 \psi \\\
y &= - \nabla^2 z
\end{align}

$$\frac{\partial}{\partial t} \begin{bmatrix} h \\\ hu \\\ hv \end{bmatrix}$$

$$\frac{\partial}{\partial x} \begin{bmatrix} hu \\\ hu^2 + \frac{1}{2}gh^2 \\\ huv \end{bmatrix}$$
$$\frac{\partial}{\partial y} \begin{bmatrix} hv \\\ huv \\\ hv^2 + \frac{1}[2}gh^2 \end{bmatrix}$$
$$ = 0$$

We note that the conserved quantities in these equations are not the physical variables \\((h, u, v)\\), but are the depth and momenta \\((h, m_x, m_y)\\), where \\(m_x = hu\\) and \\(m_y = hv\\).
