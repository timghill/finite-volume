# Shallow water equations
The shallow water solver solves the two-dimensional system of equations

\begin{align}
x &= \nabla^2 \psi \\\
y &= - \nabla^2 z
\end{align}

$$\frac{\partial}{\partial t} \begin{bmatrix} h \\\ hu \\\ hv \end{bmatrix} + \frac{\partial}{\partial x} \begin{bmatrix} hu \\\ hu^2 + \frac{1}{2}gh^2 \\\ huv \end{bmatrix} + \frac{\partial}{\partial y} \begin{bmatrix} hv \\\ huv \\\ hv^2 + \frac{1}[2}gh^2 \end{bmatrix} = 0.$$
