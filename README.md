# tjmatrix.m is the core code which provides a matrix of 10j symbol. Input is the configuration of initial state and final state which consist of 9 objects and 6 1-morphisms. 

# tj provides an element of the matrix, which is a 10j symbol. Input is the complete configuration of a 4-simplex consists of 10 objects and 10 1-morphisms.

# In lib_AssociBmm, we list all retraction associator bimodule maps. That can be regarded as part of the gauge choice.

# Interchange is the interchanger 2-morphism. That can be regarded as part of the gauge choice.

# BRbmm1, BRbmm2, FRbmm determines how associator bimodule map induce a “bigger” bimodule map, this is also related the relative tensor product between bimodules. In this case, we make a choice such that the assciator of bimodules are all trivial.

# hexagon024 and hexagon135 are two paths of the hexagon, for a given 5-simplex with fixed initial and final (all 15 objects and 18 1-morphisms are fixed and morphisms 024 and 135 non-uniquely determined), the two paths should give same map.

# config18 records all legal configurations of 5-simplex with fixed initial and final states. For each such configuration, config024 and config135 records legal configurations of morphisms 024 and 135. Once the 18 1-morphisms are fixed, the 15 objects are according fixed up to 32 choices which can be selected in possibleobj.

# In hexagon_check_total, One can set "i=1:num_high_dim" to show the two paths give same map (Sometimes differs by an amplitude sqrt(2), this is result from the choice of normalization factor) for cases that morphisms 024 and 135 have more that one choices. For the cases that morphisms 024 and 135 have only one choice, one can esaily verify the two paths give same map by changing the code hexagon_check_total.

