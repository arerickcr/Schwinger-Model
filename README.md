# Schwinger-Model

Quantum Electrodynamics in $$(1 + 1)$$-dimensions, also known as the Schwinger model, describes a system where a 2-d fermion interacts with a $$U(1)$$ gauge field. The most general action is:

$$\begin{align}
    S&=\int d^2x\left( \frac{1}{2e^2}F_{01}^2+\frac{\theta}{2\pi}F_{01}+i\bar{\psi}\gamma^{\mu}D_{\mu}\psi - m \bar{\psi} \psi \right)\notag\\
    &=\int d^2x \left( \frac{1}{2}(\partial_{\mu}\phi)^2-\frac{e^2}{2\pi}\phi^2+m e C\cos{(2\sqrt{\pi}\phi-\theta)}\right),
\end{align}$$

where in the second line we have used bosonization and $$C$$ is a known constant. When $$m\neq 0$$, the theory is not solvable, and this is why we use the help of DMRG implemented in iTensor.

In order to do so, we discretize our space and put fermions in a staggered form, which together with Jordan-Wigner transformation, gives us the following Hamiltonian:

$$\begin{align}
    H=\sum_{n=1}^{N}\left[L(n)+\frac{\theta}{2\pi} \right]^2+x\sum_{n=1}^{N-1}\left[ \sigma^{+}(n)e^{i\theta(n)}\sigma^{-}(n+1)+\textrm{c.c.}\right]+\frac{\mu}{2}\sum_{n=1}^{N}(-1)^{n}\sigma_{3}(n),
\end{align}$$

where $$\mu=\frac{2m_{latt}}{ae^2}$$, $$m_{latt}=m-\frac{e^2a}{8}$$, $$x=\frac{1}{e^2a^2}$$ and $$a$$ is the lattice spacing. Implementing Gauss' law, we can rewrite our previous expression as:

$$\begin{align}
    H& = \sum_{n=1}^{N}\left[\sum_{k=1}^{n}\frac{1}{2}\left(\sigma_3(k)+(-1)^{k}\right) + L + \frac{\theta}{2\pi} \right]^2 +x\sum_{n=1}^{N-1}\left[ \sigma^{+}(n)\sigma^{-}(n+1)+\textrm{c.c.}\right] + \notag\\
    & \quad \quad - \varepsilon x \left[ (-i)^{N}\sigma^{+}(N)e^{i\theta}\sigma^{-}(1)\prod_{\ell=1}^{N-1}\sigma_{3}(\ell) + i^{N}\sigma^{+}(1)e^{-i\theta}\sigma^{-}(N)\prod_{\ell=1}^{N-1}\sigma_{3}(\ell)\right] +\frac{\mu}{2}\sum_{n=1}^{N}(-1)^{n}\sigma_{3}(n),
\end{align}$$

with $$\varepsilon=0$$ and $$L=0$$ for Open Boundary Conditions (OBC) and $$\varepsilon=1$$ for Periodic Boundary Conditions(PBC). 

For the PBC case, since we have an extra (bosonic) degree of freedom $$L$$, then in order to describe the Hilbert space, we consider the basis $$\ket{s_{1}}\ket{s_{2}}\dots \ket{s_{N}}\ket{\ell}$$. Where $$L\ket{\ell}=\ell\ket{\ell}$$ and $$e^{\pm i\theta}\ket{\ell}=\ket{\ell \pm 1}$$ and $$\ell\in \mathbb{Z}$$. So, the Hilbert space is infinite dimensional, that means that we have to truncate our basis up to some fixed $$L_{max}$$. Thus $$-L_{max}\leq \ell \leq L_{max}$$ and there are $$(2L_{max}+1)$$ states for each chain of spins $$\ket{s_{1}}\ket{s_{2}}\dots \ket{s_{N}}$$. 

For more details on the Schwinger model and its lattice representation, I have compiled a detailed set of notes HERE, which are based on the following papers:

- [Schwinger, '62 (Gauge Invariance and Mass)] (https://doi.org/10.1103/PhysRev.125.397)
- [Kogut and Susskind '75 (Hamiltonian formulation of Wilson’s lattice gauge theories)] (https://doi.org/10.1103/PhysRevD.11.395)
- [Banks, Susskind, and Kogut, '75 (Strong Coupling calculations of lattice gauge theories: (1 + 1)−dimensional exercises)] (https://doi.org/10.1103/PhysRevD.13.1043)
- [Coleman, Jackiw, and Susskind '75 (Charge Shielding and Quark Confinement in the Massive Schwinger Model] (https://doi.org/10.1016/0003-4916(75)90212-2)
- [Coleman '76 (More about the Massive Schwinger Model)] (https://doi.org/10.1016/0003-4916(76)90280-3)
- [Hamer, Weihong, and Oitmaa, ’97 (Series expansions for the massive Schwinger model in Hamiltonian lattice theory)] (https://doi.org/10.1103/PhysRevD.56.55)
- [Byrnes '03 (Density matrix renormalization group : a new approach to lattice gauge theory] (https://doi.org/10.1103/PhysRevD.66.013002)
- [Tong '18 (Quantum Field Theory on the line)] (https://www.damtp.cam.ac.uk/user/tong/gaugetheory/72d.pdf)
- [Dempsey, Klebanov, Pufu, and Zan, '23 (Discrete chiral symmetry and mass shift in lattice hamiltonian approach to schwinger model)] (https://doi.org/10.48550/arXiv.2206.05308)
- [Hamer and Barbert, '88 (Finite-size scaling in hamiltonian field theory)] (https://doi.org/10.1016/B978-0-444-87109-1.50021-2)
- [Campostrini, Pelissetto, and Vicari, '14 (Finite-size scaling at quantum transitions)] (https://doi.org/10.1103/PhysRevB.89.094516)
- [Calabrese and Cardy, '04 (Entanglement entropy and quantum field theory)] (https://doi.org/10.1088/1742-5468/2004/06/P06002)
