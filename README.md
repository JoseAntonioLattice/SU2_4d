# Usage

Use Linux.

Open a terminal in the directory of the project and run:
```
make install
make compile
```

To run the simulation, type:
```
make run
```
Open the input_parameters.par file in the input directory and modify the parameters of the simulation. 

# Theory
We simulate the SU(2) pure gauge field theory on a 4-dimensional lattice with volume $V = L^3\times L_t$.

$U_{\mu}(x)\in\text{SU}(2)$ is a link variable at lattice point $x$ and direction $\mu$ defined in the gauge group SU(2).
We use the standard Wilson plaquette action
$$S[U] = -\frac{\beta}{N}\sum_x\sum_{\mu<\nu} \text{Re} \, \text{Tr} \left[ 1 - U_{\mu\nu}(x)\right], $$
where $\beta = 2N/g_0$. For SU(2), $N = 2$ and $g_0$ is the gauge coupling. 

The plaquette variable $U_{\mu\nu}(x)$ is defined as the ordered product
$$U_{\mu\nu}(x) = U_{\mu}(x)U_{\nu}(x+\hat\mu)U_{\mu}^{\dagger}(x+\hat\nu)U_{\nu}^{\dagger}(x).$$


