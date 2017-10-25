# timoshenko-fem

This project is a small example of codeing the finite element method. The coded
example is to determine the modal frequencies of a Timoshenko cantilever.

<!--
The kinematics of a Timoshenko beam are:

$$ u_1(x,y,z) = y\theta(x), \qquad u_3(x,y,z) = w(x). $$

With these kinematric, the non-zero strains are:

$$ 
    S_1 &= y\frac{\partial\theta}{\partial x},
    \qquad
    S_6 &= \theta + \frac{\partial w}{\partial x}.
$$

With these non-zero strains, the constitutive equation of the linear material 
become:

$$ T_1 = ES_1 \qquad T_6 = \kappa G S_6, $$
    
where $E$ is Young's modulus, $\kappa$ is a correction factor that reduces the 
shear stress in the finite element model, and $G$ is the shear modulus. The
shear modulus is a function of Young's modulus and Poisson's ratio:

$$ G = \frac{E}{2(1+\nu)}. $$
   
To develop the finite element equations, the strain and kinetic energy of the 
beam are derived and Hamilton's principle is applied. The strain energy is:

$$ U = \frac{1}{2}\int_V T_1S_1 + T_6S_6 \; dV $$

where I is the moment of inertia of the beam and A is the cross sectional area:
 
$$ I = \frac{t_zt_y^3}{12}, \qquad A = t_zt_y.$$

The shape functions are given by:

$$
    \begin{bmatrix} w \\ \theta_z \end{bmatrix}
$$

The element and mass matrices are:
Hamilton's principle must hold system described. The kinetic energy and strain energy 
$$
\begin{align}
    K_e = \int_{-1}^{1}EI\left(\frac{1}{a}\frac{dN}{d\xi})^T
        \left(\frac{1}{a}\frac{dN}{d\xi})
$$
-->

## References

Matlab Codes for Finite Element Analysis; A. J. M. Ferreira; 2009; Springer
