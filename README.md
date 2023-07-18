# A gene-regulatory network model for density-dependent and sex-biased dispersal evolution during range expansions.
Code for:

Jhelam N. Deshpande and Emanuel A. Fronhofer: A gene-regulatory network model for density-dependent and sex-biased dispersal evolution during range expansions.

Abstract: 

Eco-evolutionary dynamics of range expansions are of great interest in the context of global change. Range expansions are spatial phenomena and critically depend on organisms' dispersal abilities. Dispersal has a genetic basis, it can evolve, but also be plastic to external environmental conditions and the internal state of an organism. Importantly, dispersal plasticity itself can evolve which has most often been studied theoretically with a focus on optimal reaction norms under equilibrium conditions. However, under rapidly changing conditions, the rate of dispersal plasticity evolution will impact eco-evolutionary dynamics of range expansions. Rates of evolution in turn depend on the genetic architecture of underlying traits. To elucidate this interplay, we develop an individual-based metapopulation model of the evolution of density-dependent and sex-biased dispersal during range expansions. We represent the dispersal trait as a gene-regulatory network (GRN), which can take population density and an individual's sex as an input and analyse emergent context- and condition-dependent dispersal responses. We compare dispersal evolution and ecological dynamics in this GRN model to a standard reaction norm (RN) approach under equilibrium metapopulation conditions and during range expansions. We find that under equilibrium metapopulation conditions, the GRN model produces emergent density-dependent and sex-biased dispersal plastic response shapes that match theoretical expectation of the RN model. However, during range expansion, the GRN model leads to faster range expansion because GRNs can maintain higher adaptive potential. Our results imply that, in order to understand eco-evolutionary dynamics in contemporary time, the genetic architecture of traits must be taken into account.

Description:

In each folder there is the associated cpp code range_expansion.cpp, the parameter input file and the Makefile corresponding to the models described below. The cpp code is for individual-based simulations of dispersal plasticity evolution and range expansions in a metapopulation model.

The code range_expansion.cpp in DDD/GRN assumes a GRN model for density-dependent dispersal
The code range_expansion.cpp in DDD/PH assumes a reaction norm according to Poethke & Hovestadt (2002) for density-dependent dispersal
The code range_expansion.cpp in DDD_sex_bias/GRN assumes a GRN for density-dependent and sex-biasd dispersal
The code range_expansion.cpp in DDD_sex_bias/PH assumes a reaction norm according to Poethke & Hovestadt (2002) for density-dependent and sex-biased dispersal
