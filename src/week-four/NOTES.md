## Aiyagari Model (Standard incomplete market)

### Decentralization

Solutions are not equivalent to planner problem, `price -> individuals -> decision -> prices -> (...)`. So solved in two step:

1. Given $p$ solve for decision rule
2. Given decision (heterogenous agents) solve for $p$

### Heterogenous agents

Consumers are different in wealth $a$, which yields aggregate demand $\lambda(a)$ which in turns yields prices $a_t \mapsto p(\cdot, \lambda(a_t))$.

### Partial equilibrium: Income fluctuation

$$
V(a, y) = max_{c, a'} u(c) + \beta \sum_{y'\in Y} \Gamma(y' \vert y) \cdot V(a', y')
$$

s.t.

$$
c + a' \leq \underbrace{(1+r)}_{ex. in partial}\cdot a + w\cdot y
$$

where $y$ (productivity) follows the markov process $\Gamma$.

We can use usual methods to solve for $c(a, y)$ and $a'(a, y)$. Note that we have two laws of motion for $(a, y)$ of which we need to find the **stationary distribution**.

We call the measure 

$$
\lambda \ : \ A \times Y \xRightarrow{} [0, 1]
$$

*stationary* if, 

$$
\lambda(\mathcal{A} \times \mathcal{Y}) = \int_{A\times Y} Q((a, y), \mathcal{A} \times \mathcal{Y}) \lambda(a, y)
$$

such that,

$$
Q((a, y), \mathcal{A} \times \mathcal{Y}) = \mathbf{1}\{ a' \in \mathcal{A} \} \cdot \sum_{y' \in \mathcal{Y}} \Gamma(y' \vert y)
$$

In practice we want to find a stationary $\lambda$. Algorithm,

1. Find $c(a, y)$ and $a'(a, y)$ - (note that the constraint can be occasionally binding if $a'(a, y) > \underline{a}$)
2. For every $\mathcal{A} \times \mathcal{Y}$ construct the transition $Q((a, y), \mathcal{A} \times \mathcal{Y})$
3. Find stationary for $\lambda$

### General Equilibrium

Look for market clearing prices. Let aggregate assets,

$$
A(r) = \int a \cdot \lambda(a, y \vert r)
$$ 

If the model is *Huggett*, market clears at $A(r) = 0$.

If the model is $Aiyagari$, firms demand funds, such that

$$
K(r): r + \delta = F_k(K, L) \\
w = F_L(K, L)
$$

hence $A(r) = K(r)$.

Use root-finding to solve for $r$ (or $w$). 

### Calculating $\lambda(a, y)$

Let $Y \sim \Gamma$ be the state space and $\gamma: \gamma^{T} = \gamma^{T}\Gamma$, found with the left eigenvector (with eigenvalue $1$).

#### Iterative method

Then the invariant cumulative asset distribution, given $\hat{y}\in Y$, with $g(a \vert y) = a'(a, y)$,

$$
\Lambda(a', y') = \sum_{y \in Y} \Gamma(y' \vert y) \cdot \Lambda\left( g^{-1}(a' \vert y), y \right)
$$

Note, that we assume $\frac{\partial}{\partial a} g(a, y) > 0$.

How? Start (1) with spreading the agents uniformly across the asset space,

$$
A_0(a_k, y_j) = \frac{a_k - \underline{a}}{\overline{a} - \underline{a}} \cdot \gamma
$$

then update distribution on grid points (2) such that,

$$
A_1(a_k, y_j) = \sum_{y_i \in Y} \Gamma(y_j \vert y_i) \cdot \Lambda_0\left( g^{-1}(a_k \vert y_i), y_i \right)
$$

We can get $g^{-1}(a_k \vert y_i)$ solving non-linearly as $a_k = a'(a, y_i)$. Note, that $\Lambda_i$ needs to be interpolated. 

#### Eigenvector method

We are looking for the density instead of the cumulative.

For a given transition matrix, the invariant distribution $\lambda' = \lambda' Q$, is associated with the left-eigenvector of $Q$ (re-normalized to sum to $1$).

The overall transition $Q \in R(NJ \times NJ)$ satisfies, 

$$
Q((a', y') \vert (a, y)) = Q_a(a' \vert a, y) \otimes \Gamma
$$

Note that $Q$ is sparse, since $Q(a' \vert a, y) = 0 \ \ \forall a' \notin [a_k, a_{k+1}]$. Hence one needs to find the interval in which $a'$ falls. Then,

$Q_a(a_k \vert a, y) = \frac{a_{k+1} - a'}{a_{k+1} - a_k}$.