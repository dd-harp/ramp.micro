# **`ramp.micro`** 

# Microsimulation for Mosquito Ecology  & Pathogen Transmission

This software was developed to lower the costs of building and analyzing micro-simulation models for mosquito ecology and for mosquito-borne pathogen transmission dynamics and control.
It is a GitHub R Package that can be installed using `devtools` by running this code from an R prompt:

```
library(devtools)
devtools::install_github("dd-harp/ramp.micro")
```

After it has been installed, it can be loaded in the normal way: 

```
library(ramp.micro)
```

*** 

![**Figure 1:** In microsimulation models, mosquitoes are in behavioral states (*e.g.* egg laying, *L,* or blood feeding, *F*) and move among points where the required resources are found (*e.g.* aquatic habitats, $\left\{l \right\},$ or blood feeding sites, $\left\{b \right\}$). This diagram was modified from Perkins TA, *et al.* 2013.^[Perkins TA, Scott TW, Le Menach A, Smith DL (2013). Heterogeneity, mixing, and the spatial scales of mosquito-borne pathogen transmission. PLoS Comput Biol 9:e1003327, https://doi.org/10.1371/journal.pcbi.1003540]](./vignettes/DynamicsOnPoints.png)

*** 


## Microsimulation & Behavioral States 

This software was developed to explore mosquito population dispersal, the spatial dynamics of mosquito populations, and malaria transmission dynamics on *point sets,* which we call *micro-simulation.*
The idea of micro-simulation was first described for malaria by the late Richard Carter (Carter, 2002)^[Carter R (2002) Spatial simulation of malaria transmission and its control by malaria transmission blocking vaccination. International Journal for Parasitology 32: 1617–1624. doi:10.1016/S0020-7519(02)00190-X]. 

This software also implements several *behavioral state models* for mosquitoes, which are a natural complement to micro-simulation. 
The basic premise of is that mosquitoes are in a physiological state that triggers behaviors, including behavioral algorithms to *search* for resources. 
The search algorithms have evolved to accomplish tasks, and they may have evolved to increase a mosquito's fitness in context. 
To survive and lay eggs, mosquitoes must find resources and accomplish several tasks: find mates; find and blood feed on vertebrate hosts; find aquatic habitats and lay eggs; and find and feed on sugar sources for energy. 
Other states could include a resting period, such as the post-prandial rest to concentrate a blood meal. 
The mosquito is, meanwhile, responding to avoiding hazards by responding to environmental cues such as wind, temperature, and humidity. 

In these models, 
Behavioral state models for mosquitoes are a type of compartmental model where mosquito populations are sub-divided by their physiological / behavioral states, and changes in behavioral states can occur as a result of successfully blood feeding, egg laying, sugar feeding, mating, and resting. 
These behavioral states are different than the Ross-Macdonald model's states that represent infection status for parasite transmission dynamics: uninfected, infected but not yet infective, and infective. By considering *both* the physiological/behavioral state and infection states, it might be possible to understand how local features of malaria transmission could be the result mosquito behavioral ecology and the heterogeneous distribution of resources, an idea pioneered by Arnaud Le Menach, *et al.* (2005)^[Le Menach A, McKenzie FE, Flahault A, Smith DL (2005) The unexpected importance of mosquito oviposition behaviour for malaria: Non-productive larval habitats can be sources for malaria transmission. Malar J. 4: 23, doi:10.1186/1475-2875-4-23].

A mosquito search for resources ends at a location where the resources can be found. 
The point sets in these micro-simulation models thus represent the locations of different resources that mosquitoes need. 
The first paper to combine *micro-simulation* with *behavioral state modeling* was an agent based model published by Weidong Gu and Robert J Novak (2009 a,^[Gu W,  Novak RJ (2009). Agent-based modelling of mosquito foraging behaviour for malaria control, Trans R Soc Trop Med Hyg 103: 1105–1112, https://doi.org/10.1016/j.trstmh.2009.01.006] b^[Gu W, Novak RJ (2009). Predicting the impact of insecticide-treated bed nets on malaria transmission: the devil is in the detail. Malar J 8:256, https://doi.org/10.1186/1475-2875-8-256]).
A rigorous mathematical framework to describe mosquito behavioral state micro-simulation was developed by Alex Perkins, *et al.*, (2013)^[Perkins TA, Scott TW, Le Menach A, Smith DL (2013). Heterogeneity, mixing, and the spatial scales of mosquito-borne pathogen transmission. PLoS Comput Biol 9:e1003327, https://doi.org/10.1371/journal.pcbi.1003540]. An individual-based model, 
called MBITES (Mosquito Bout-based and Individual-based Transmission Ecology Simulator), 
was developed by Sean L Wu, *et al.* (2020)^[Wu SL, *et al.*, (2020). Vector bionomics and vectorial capacity as emergent properties of mosquito behaviors and ecology. PLoS Comput Biol 16:e1007446, doi:10.1371/journal.pcbi.1007446]. 
This list of papers is by no means exhaustive, and a systematic review of such models is needed.

## Complexity 


These behavioral state models are designed to be highly mimetic. 
While searching, mosquitoes move around until they find a resource and accomplish the task(s).
The physiological state then changes, and the mosquitoes must find a different resource.
In searching for resources, the wind speed and direction are major concerns that could affect mosquito searching efficiency and the distance traveled during a single flight.
Mosquito movement and mosquito population dynamics are thus determined by wind, behavior, and the distribution of resources.
Mosquitoes are moving among resource point sets distributed on landscapes, so while standard approaches based on the mathematics of *diffusion* might seem reasonable (for example, see Lutambi ML, *et al*, 2013)^[Lutambi ML, *et al.* (2013).
Mathematical modelling of mosquito dispersal in a heterogeneous environment.
Mathematical Biosciences 241:198-216, https://doi.org/10.1016/j.mbs.2012.11.013], the underlying mathematics governing aggregate movement patterns of mosquito populations are related to both search and to the co-distribution of resources. 
The emerging patterns are highly structured and complex. 
This is not to say that diffusion-based models are not *useful,* but that we might learn something new by describing and analyzing mosquito movement in models that are highly realistic. 
*Spatial Dynamics of Malaria Transmission,* for example, introduces a meta-population model for mosquito spatial ecology, where the patch-based diffusion model is motivated by the ideas of *search* and heterogeneous resource availability (Wu SL, *et al.*, 2023)^[Wu SL, *et al.* (2023) Spatial dynamics of malaria transmission. PLoS Comput Biol 19(6): e1010684. https://doi.org/10.1371/journal.pcbi.1010684]; the approach was motivated by MBITES (Wu SL, *et al.*, 2020)^[Sean L Wu, *et al.* (2020). Vector bionomics and vectorial capacity as emergent properties of mosquito behaviors and ecology. PLoS Comput Biol 16:e1007446, doi:10.1371/journal.pcbi.1007446]

It is challenging to set up, simulate, analyze, and visualize the outputs of micro-simulation models. This software, `ramp.micro`, was developed lower the human costs of model development with the hope of encouraging new research. 
We hope this software can be used to advance our understanding of mosquito spatial ecology and the spatial dynamics of malaria and other mosquito-transmitted diseases. There are several promising areas of research: understanding intervention *coverage* as a spatial concept; the *effect sizes* understood spatially, especially in relation to coverage when those interventions *repel* mosquitoes from an area; spatial *area effects* of vector control; mosquito spatial dynamics; *spatial targeting* for malaria control; sampling mosquito populations and the robustness of the metrics used to measure mosquito bionomic parameters; theory to support *larval source management;* the spread of genetically modified mosquitoes and effect sizes; and finally, the fundamental concept of a mosquito *niche.*

The software has a modular design, and there are plans to integrate some of this functionality with [`exDE`](https://dd-harp.github.io/exDE/). Please contact us if you would like to  be involved with the project, or if you have any questions or suggestions.  


***

