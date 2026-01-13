# Adult Mosquitoes - BQ

The `BQ` model describes a discrete-time, behavioral state model for
adult mosquitoes moving on point sets. The model has two states: blood
feeding $B,$ and egg laying $Q.$

## Point Sets

We define this class of micro-simulation models on point sets
representing the locations of resources: haunts where mosquitoes rest
and where a mosquito might find and feed on a vertebrate host; and
aquatic habitats where mosquitoes could find aquatic habitats and lay
eggs.

- Let $\left\{ b \right\}$ denote the point set where blood feeding
  could occur;

- Let $\left\{ q \right\}$ denote the point set where egg laying could
  occur;

## Dispersal

Movement among point sets is modeled using matrices that describe where
mosquitoes move each time step. We assume that every surviving mosquito
ends up somewhere. Mortality while searching is referenced to the point
where the search starts.

The proportion moving from a point in one set to another, from
$x \in \left\{ x \right\}$ to $y \in \left\{ y \right\}$, is described
by a matrix $\Psi_{y\leftarrow x}$ or equivalently $\Psi_{yx}$.
Similarly, the proportion moving from a point in one set to a point in
the other is $\Psi_{xx}$. Since there are a finite number of
destinations, each column is a probability mass function (PMF).

- Let $\Psi_{b\leftarrow q}$ denote a matrix describing the location
  where mosquitoes end their flights in $\left\{ b \right\}$ starting
  from each point in $\left\{ q \right\}$.

- Let $\Psi_{b\leftarrow b}$ denote a matrix describing the location
  where mosquitoes end their flights in $\left\{ b \right\}$ starting
  from each point in $\left\{ b \right\}$.

- Let $\Psi_{q\leftarrow b}$ denote a matrix describing the location
  where mosquitoes end their flights in $\left\{ q \right\}$ starting
  from each point in $\left\{ b \right\}$.

- Let $\Psi_{q\leftarrow q}$ denote a matrix describing the location
  where mosquitoes end their flights in $\left\{ q \right\}$ starting
  from each point in $\left\{ q \right\}$.

## Variables & Terms

In the simulation models, the number of adult mosquitoes emerging from
each aquatic habitat on each day is a vector denoted $\Lambda_{t}$ that
can be either passed as a parameter or simulated using a population
dynamic model for aquatic immature population dynamics. Adult mosquitoes
are located at one of these points at each point at each point in time
(*e.g.*, one day), represented as a set of vectors: $B_{t}$, the number
seeking blood at haunts; or $Q_{t}$, the number attempting to lay eggs
at habitats.

- Let $B_{t}$ denote the number of *blood feeding* mosquitoes at
  $\left\{ b \right\}$

- Let $Q_{t}$ denote the number of *egg laying* mosquitoes at
  $\left\{ q \right\}$

- Let $\Lambda_{t}$ denote the number of recently emerged mosquitoes at
  $\left\{ q \right\}$

## Parameters

Mosquito bionomics can depend on both behavioral state and location,
including daily survival and daily blood feeding or egg laying success.
Note that survival is linked to the point where the search started, but
foraging success is linked to the point where the mosquito ends up after
moving in a time step.

**Survival:**

- Let $p_{b}$ denote the probability of surviving for a mosquito at
  $\left\{ b \right\}$ at time $t$

- Let $p_{q}$ denote the probability of surviving for a mosquito at
  $\left\{ q \right\}$ at time $t$

**Foraging Success:**

- Let $\psi_{b}$ denote the probability of surviving for a mosquito at
  $\left\{ b \right\}$ at time $t$ and let
  ${\widehat{\psi}}_{b} = 1 - \psi_{b}$

- Let $\psi_{q}$ denote the probability of surviving for a mosquito at
  $\left\{ q \right\}$ at time $t$ and let
  ${\widehat{\psi}}_{q} = 1 - \psi_{q}$

## Dynamics

In the basic feeding cycle model over one time step, mosquitoes either
attempt to blood feed or attempt to lay eggs. The result of an attempt
is either survival or death, and if the mosquito survives, success or
failure. A success moves a mosquito to the other state, and if they
fail, they must try again. Either way, the mosquito moves to a point in
the set for the resource they seek: the diagonal of
$\Psi_{x\leftarrow x}$ is the probability of staying. The dynamics are:

$$\begin{bmatrix}
B_{t} \\
Q_{t} \\

\end{bmatrix} = \begin{bmatrix}
{\Psi_{bq}p_{q}\Lambda_{t - 1}} \\
0 \\

\end{bmatrix} + \begin{bmatrix}
{\Psi_{bb} \cdot \text{diag}\left( {\widehat{\psi}}_{b} \right)} & {\Psi_{bq} \cdot \text{diag}\left( \psi_{q} \right)} \\
{\Psi_{qb} \cdot \text{diag}\left( \psi_{b} \right)} & {\Psi_{qq} \cdot \text{diag}\left( {\widehat{\psi}}_{q} \right)} \\
 & 
\end{bmatrix}\begin{bmatrix}
{p_{b}B_{t - 1}} \\
{p_{q}Q_{t - 1}} \\

\end{bmatrix}$$
