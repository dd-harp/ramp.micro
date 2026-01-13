# Immature Mosquitoes - L

## 

For purposes of integrating the effects of movement and aquatic ecology,
we assume each female lays $o$ eggs, so the total number of eggs laid
each day in each site is:

$$\eta_{t} = o\psi_{q}Q_{t}$$

We assume the number of immature mosquitoes in aquatic habitats, denoted
$L_{t}$, is subject to density dependence, which could delay maturation
or increase mortality. The maturating fraction is,
$\theta e^{- \xi L_{t}}$ where $\xi > 0$, and surviving fraction is
$p_{L}e^{- \zeta L_{t}}$ where $\zeta > 0$. The parameters are
site-specific, so that some habitats can vary in quality. The dynamics
are:

$$L_{t + 1} = \eta_{t} + p_{L}e^{- \zeta L_{t}}\left( 1 - \theta e^{- \xi L_{t}} \right)L_{t}.$$

The number of adult females emerging each day is half the mosquitoes who
both survived and matured:

$$\Lambda_{t + 1} = \frac{p_{L}\theta e^{- {(\xi + \zeta)}L_{t}}}{2}L_{t}$$

While this model is adequate for the needs of this study, the model may
not robustly capture some relevant features of mosquito ecology or
larval source management, such as when larval population structure and
delays are important.
