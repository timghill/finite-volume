# Shallow water equations
The shallow water solver solves the two-dimensional system of equations

{% \begin{align*}
\frac{\partial}{\partial t} \begin{bmatrix} h \\ hu \\ hv \end{bmatrix} + \frac{\partial}{\partial x} \begin{bmatrix} hu \\ hu^2 + \frac{1}{2}gh^2 \\ huv \end{bmatrix} + \frac{\partial}{\partial y} \begin{bmatrix} hv \\ huv \\ hv^2 + \frac{1}{2} gh^2 \end{bmatrix}.
\end{align*} %}

\begin{align*}
x &= \nabla y \\
y &= -\nabla x \\
end{align*}
