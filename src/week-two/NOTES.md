# Markov processes

## Simple example - Discrete

Let the transition function be `Q(z, z')`, then `P`: `P[i, j] = Q(i, j)`.

A stationary distribution is given by `z_t = π(0) * P^t`.

such that the stationary distribution is `π * (I - P) = 0`, by solving the **eigenvalue problem**. 

### Eigenvalue problem

If the probability matrix does not imply that any state can be reached by another, then it is not solvable.

### Stability

A Markov chain is **irreducible** and **aperiodic** if `P(z | any) > 0`, i.e. you can always end up in another state given your initial one.

A markov chain is **stable** if `π(0) P^t -> π` for every `π(0)`. (*Note*: the transition matrix is a contraction mapping!)

## General

- Markov process forgets history: `E[f(X(t)) | H] = E[f(X(t)) | X(t-2)]`
- Every time we construct a process such that `X(t+1) ~ Q` then it is a Markov process

Let a stochastic difference equation, any equation such that `z(t+1) = g(z(t), w(t))`.

### E.g.

- *AR(1)*: `z(t+1) = r*z(t) + e(t+1)`, where `e(t+1) ~ N(0,1)`. Then note that `Q(z, A) = Φ(a - r z | a in A)`  

- *Log normal*: `log(z(t+1)) = r * log(z(t)) + e(t+1)`, where `e(t+1) ~ N(0,1)`. Then,  `Q(z, A) = Φ(log(a) - r*log(z) | a in A)`  

## Stability

The Markov operator (conditional expectation),

$$
T{f}(z) = \int{f(z') \cdot Q(z, z') dz}
$$

Then, as before, stability requires $$T^{(t)}\{\pi_0\} \xRightarrow{t \xRightarrow{} \infty} \pi$$

Furthermore,

- Sufficient mixing, for any state, every other state needs to have a positive probability. - Use the Dobrushin coefficient to measure the commonality between two distribution (do not want it to be small). In particular, a **positive Dobrushin** coefficient implies **nonexpansiveness**.
- No divergence - **Tightness**
- No state should accumulate probability mass (no mass point) - **Uniform Integrability**

# Equilibrium of RGB

Let `z(t)` be a Markov process. Then the Bellman equation,

```
V(z, k) = max{ u(c) +  β E[V(z', k')] } s.t. c + k' - (1-d)k = z f(k)
```

is a contraction mapping iff `z` is stable.

## Discretization

### Tauchen (1986) - Fix size, different mass

Consider *AR(1)* for illustration. Note that it has variance, `s^2 / (1 - r^2)`. Objective is to find, `{z_1, ..., z_N}` discretization from real line and associate probability.

We set 

```julia
z_N = m √(s^2 / (1 - r^2))
```

and, if symmetric, `z_1 = -z_N` (note that this depends on the form of the distribution), and take a perfect partition between `z_1` and `z_N`.

Then the probability transition matrix, is the probability of falling around the interval `z_i`, conditional on being in a given `z`, in this case,

```
P[i, 1] = F((z_i + d/2 - r * z_i) / s) # first row
P[i, j] = F((z_i + d/2 - r * z_i) / s) -  F((z_i - d/2 - r * z_i) / s) # i-th row
P[i, N] = 1 -  F((z_i - d/2 - r * z_i) / s) # last row
```

### Tauchen, importance sampling -  Fix mass, different size

Instead of using equispace grid point, we look at the same probability for every interval, and change the partition. Then the bin bound are,

$$
b_i = F^{-1}(i / N), \ \ i \in {0, 1, ..., N}
$$

such that the grid points are, 

$$
z_i = E[x | x in (b_{i - 1}, b_i)]
$$

### Rouwenhorst (1995)

Rouwenhorst has better accuracy then Tauchen when the process has high persistence. Construct the chain as follows,

$$P_2 = \begin{pmatrix} p & 1 - p \\ 1 - q & q \end{pmatrix}$$

and then use the recursive definition,

$$P_N = p \begin{pmatrix} P_{N-1} & 0 \\ 0 & 0 \end{pmatrix} + (1-p) \begin{pmatrix} 0 & P_{N-1} \\ 0 & 0 \end{pmatrix} + (1-q) \begin{pmatrix}0 & 0 \\  P_{N-1} & 0 \end{pmatrix} + q \begin{pmatrix}0 & 0 \\  0 & P_{N-1} \end{pmatrix} $$

then use a recursive system. Then it has been shown that, the discretised version is a binomial distribution that needs to be matched to the target distribution using the first two moments.

## Integration

Consider a partition [a, b],

### Midpoint rule

The rule yields,

$$\int_{x_i - h/2}^{x_i + h/2} f(x) dx \approx f(x) \cdot h $$

yields cumulatively.

### Trapezoid rule

Pick the endpoints, and use the linear approximation, 

$$\tilde{f}(x)$$

and compute the below surface,
 
$$\int_{x_i}^{x_{i+1}} f(x) dx \approx \frac{h}{2} \left[f(x_i) + f(x_{i+1})]\right) $$

### Simpson's rule

Approximate the integral with a quadratic function over a triplet of points.

$$\int_{x_{i-1}}^{x_{i+1}} f(x) dx \approx \frac{h}{3} \left[f(x_{i-1}) + f(x_{i+1}) + 4 \cdot f(x_i)]\right) $$

### Gaussian Quadrature

We can choose nodes more accurately (increases accuracy from `n - 1` to `2n - 1`).

Orthogonal polynomials are such that (Judd, 1998)

$$
\langle\phi_k, \phi_j \rangle = \int_{a}^b \phi_k(x) \cdot \phi_k(x) \cdot w(x) dx = 0
$$

then,

$$
\int_{a}^b f(x) \cdot w(x)   dx \approx \sum^n_{i=0} w_i f(x_i) + \frac{f^{(2n)(\xi)}}{q^2_n(2n)!}
$$

where $$x_i: \Psi(x_i) = 0,  \Psi \text{is orthonormal for } w(x)$$

### Gauss-Chebyshev

Are function orthonormal to $$(1-x^2)^{-0.5}$$ so that,

$$
\int_{-1}^{1} f(x) (1-x^2)^{-0.5} dx \approx \sum^{n}_{i=1} w_i f(x_i), \text{ where } w_i = \frac{\pi}{n}, x_i = cos\left(\frac{2i - 1}{2n} \pi\right)
$$

Note that this, allows as to map any, 

$$
\int^b_{a} f(y) dy \mapsto \int^1_{0} f(x) w(x) dx
$$

In this case the mapping is,

$$
x = \left(y \mapsto -1 + \frac{2(y-a)}{b-a} \right)
$$

# Problem set 2

## a.

Given an *AR(1)* process,

$$
z_{t+1} = \rho \ z_t + (1- \rho^2) \ \epsilon_{t+1}
$$

Then, 

$$
V[z_{t+1}] = \rho^2 \ V[z_t] +  (1- \rho^2)^2 \cdot \sigma^2\\
V_z = \rho^2\cdot V_z +   (1- \rho^2)^2 \cdot \sigma^2 \\
 (1- \rho^2) \cdot  V_z = (1- \rho^2)^2 \cdot \sigma^2 \\
 V_z =  (1- \rho^2) \cdot \sigma^2
$$
