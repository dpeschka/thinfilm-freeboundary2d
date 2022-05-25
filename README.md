# Time-dependent thin-film free boundary problems with dynamic contact angle

MATLAB and deal.ii code complementing the data in the JCP paper "Model hierarchies and higher-order discretisation of time-dependent thin-film free boundary problems with dynamic contact angle". The codes in this repository solve a similar problem as the one stored in the repository

https://github.com/dpeschka/thinfilm-freeboundary

but in higher dimensions using higher-order methods in space and in time.

The general problem solved is to find a function $h(t):\omega\to\mathbb{R}$ and a time-dependent set $\omega(t)\subset\mathbb{R}^2$ such that

$$
\partial_t h - \nabla\cdot (m(h)\nabla\pi) = 0, \qquad \pi=-\sigma\nabla^2 h + \partial_h (hf)
$$

for all $x\in\omega(t)$ with natural boundary conditions $\nu\cdot\nabla\pi=0$ on $\partial\omega(t)$ supplemented with a dynamic contact angle

$$
\partial_t h + n(\nabla h)|\nabla h|^2\zeta^2 =0, \qquad \zeta=|\nabla h|^{-1}(-\tfrac{\sigma}{2}|\nabla h|^2 + s)
$$

for all $x\in\partial\omega(t)$ with kinematic condition $\partial_t h + \dot{x}\cdot\nabla h=0$ emerging from $h(t)=0$ on $\partial\omega(t)$. This evolution is a gradient flow for the energy

$$
\mathcal{E}(q) = \int_\omega \tfrac{\sigma}{2}|\nabla h|^2 + s + h f(x,h)\,{\rm d}x,
$$

with $q=(\omega,h)$. In the paper we explain the derivation of different weak formulations for this problem and different limits of weak and strong contact line dissipation.
