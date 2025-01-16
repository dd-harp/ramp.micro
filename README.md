# **`ramp.micro`**  -- Mosquito Microsimulation 

![**Figure 1:** Modified from Perkins TA, *et al.* 2013.^[Perkins TA, Scott TW, Le Menach A, Smith DL (2013). Heterogeneity, mixing, and the spatial scales of mosquito-borne pathogen transmission. PLoS Comput Biol 9:e1003327, https://doi.org/10.1371/journal.pcbi.1003540]](./vignettes/DynamicsOnPoints.png)

*** 

## Download & Install

**`ramp.micro`** is a GitHub R Package.
It was developed to lower the costs of building and analyzing micro-simulation models for mosquito ecology and for mosquito-borne pathogen transmission dynamics and control.
It can be downloaded and installed using `devtools` by running this code from an R prompt:

```
library(devtools)
devtools::install_github("dd-harp/ramp.micro")
```

After it has been installed, it can be loaded in the normal way: 

```
library(ramp.micro)
```

## Microsimulation & Behavioral States 

This software was developed to explore complexity in mosquito ecology through the development and analysis of simulation models describing mosquito behavioral & spatial dynamics on *point sets* representing resources, which we call *micro-simulation.* In this software, we develop models that combine micro-simulation with behavioral state models for mosquitoes. 


Behavioral state models for mosquitoes are a type of compartmental model where mosquito population states include their *behavioral* states in addition to a mosquito's *infection* states. Behavioral dynamics describe changes that occur as mosquitoes search for resources and then change behavioral states after they successfully blood feed, lay eggs, sugar feed, mate, or rest. 
These behavioral states are distinct from the *infection* state space; infection occurs after mosquitoes successfully blood feeds on an infectious host and gets infected. They are a complement to the Ross-Macdonald model and other models that track infection dynmics. 
The first behavioral state models for malaria were developed as patch-based / metapopulation models by Arnaud Le Menach, *et al.* (2005)^[Le Menach A, McKenzie FE, Flahault A, Smith DL (2005) The unexpected importance of mosquito oviposition behaviour for malaria: Non-productive larval habitats can be sources for malaria transmission. Malar J. 4: 23, doi:10.1186/1475-2875-4-23]. 
By considering *both* the physiological/behavioral state and infection states, it might be possible to understand how local features of malaria transmission could be the result mosquito behavioral ecology and the heterogeneous distribution of resources.


The idea of micro-simulation was first described for malaria by the late Richard Carter (Carter, 2002)^[Carter R (2002) Spatial simulation of malaria transmission and its control by malaria transmission blocking vaccination. International Journal for Parasitology 32: 1617–1624. doi:10.1016/S0020-7519(02)00190-X]. 
The first behavioral state microsimulation model was an individual-based simulation model 
described in a pair of papers by
Weidong Gu and Robert J Novak (2009 a,^[Gu W,  Novak RJ (2009). Agent-based modelling of mosquito foraging behaviour for malaria control, Trans R Soc Trop Med Hyg 103: 1105–1112, https://doi.org/10.1016/j.trstmh.2009.01.006] b^[Gu W, Novak RJ (2009). Predicting the impact of insecticide-treated bed nets on malaria transmission: the devil is in the detail. Malar J 8:256, https://doi.org/10.1186/1475-2875-8-256]).
A rigorous mathematical framework to describe mosquito behavioral state micro-simulation was developed by Alex Perkins, *et al.*, (2013)^[Perkins TA, Scott TW, Le Menach A, Smith DL (2013). Heterogeneity, mixing, and the spatial scales of mosquito-borne pathogen transmission. PLoS Comput Biol 9:e1003327, https://doi.org/10.1371/journal.pcbi.1003540]. An individual-based model, 
called MBITES (Mosquito Bout-based and Individual-based Transmission Ecology Simulator), 
was developed by Sean L Wu, *et al.* (2020)^[Wu SL, *et al.*, (2020). Vector bionomics and vectorial capacity as emergent properties of mosquito behaviors and ecology. PLoS Comput Biol 16:e1007446, doi:10.1371/journal.pcbi.1007446]. 

This software implements several *behavioral state, micro-simulation models* for mosquitoes. The basic premise of is that mosquitoes are in a physiological state that triggers behaviors, including behavioral algorithms to *search* for resources. 
The resources are distributed as points on a landscape. 
Mosquito search algorithms have evolved to accomplish tasks, and they may have evolved to increase a mosquito's fitness in response to specific features of a landscape and the local ecology. 
To survive and lay eggs, mosquitoes must find resources and accomplish several tasks: find mates; find and blood feed on vertebrate hosts; find aquatic habitats and lay eggs; and find and feed on sugar sources for energy. 
Other states could include a resting period, such as the post-prandial rest to concentrate a blood meal. 

## Complexity & Robust Analytics

In developing this software, we have two goals. 

1. We hope to develop a platform to facilitate an understanding of 
mosquito ecology as an emergent feature of individual mosquito behaviors 
constrained by a biotic and abiotic environment. In taking an approach to mosquito 
ecology that focuses on complexity, we hope to improve our ability to understand and model mosquito biology, 
the role they play in pathogen transmission, and the prospects for vector control.

2. We hope this platform can facilitate studies to explore uncertainty about mosquito ecology and its role in evaluating policies to manage mosquito-borne pathogen transmission.
By taking a deep dive into mosquito ecology and behavior, we seek to 
develop *robust analytics for malaria policy* (RAMP) as well as robust analytics to support other 
mosquito-transmitted pathogens. 
We thus see **`ramp.micro`** as integral to [**RAMP** and **Adaptive Malaria Control**](https://faculty.washington.edu/smitdave/ramp/), especially **Adaptive Vector Control.**
It is one of several R packages developed for *Sim*ulation *B*ased *A*analytics
[(**SimBA**)](https://faculty.washington.edu/smitdave/simba/index.html).

It is challenging to set up, simulate, and analyze micro-simulation models. 
It is also difficult to visualize the outputs. 
This software, **`ramp.micro,`** aims developed lower the human costs of model developing, analysis, and visualization,
and to serve as a platform for developing studies that could be replicated by others. 
To avoid constraining research, the software was designed to be modular, flexible, and extensible.  

We hope this software can be used to advance our understanding of mosquito spatial ecology and the spatial dynamics of malaria and other mosquito-transmitted diseases. 
There are several promising areas of research: understanding intervention *coverage* as a spatial concept; the *effect sizes* understood spatially, especially in relation to coverage when those interventions *repel* mosquitoes from an area; spatial *area effects* of vector control; mosquito spatial dynamics; *spatial targeting* for malaria control; sampling mosquito populations and the robustness of the metrics used to measure mosquito bionomic parameters; theory to support *larval source management;* the spread of genetically modified mosquitoes and effect sizes; and finally, the fundamental concept of a mosquito *niche.*

If you are interested in collaborating, please contact [Professor David L Smith, University of Washington](https://faculty.washington.edu/smitdave/).
