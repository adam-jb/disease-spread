# Simulations of the efficacy of herd immunity with various parameters

Functions to model a network and the spread of disease through that network over time. 

This allows the user to output results for different combinations of parameters:
-how clustered the immune population are
-proportion of population immune
-number of initial cases of the disease
-r value of disease (transmission)
-number of days someone is infectious for

Source this code for two primary functions:

1) initialise_network() creates a network with 4 levels of connectively (dunbar layers). The user can set the no. of connections at each of the 4 levels, with the default values roughly reflecting Dunbar's hypothesis (5, 15, 40, 90). Everyone in the network will have the same number of connections at each level (this is clearly not a reflection of the real world!)

This runs in under 5 minutes for a 10,000 person network (2012 Macbook Pro), however runtime increases exponentially with size. The network only needs to be initialised once, so a large network could be saved for all disease simulations. 

Inputs are all single values.

Returns a dictionary object of all network connections, which an the input for the next function:


2) sim_spread_mult() determines who is immune under herd immunity, then models the spread of the disease multiple times. 

Inputs are initialise_network()'s dictionary object and a mix of 1D arrays and single values.

Returns an array with a row for each combination and iteration, where the rightmost column is the proportion of the non-immune population who caught the disease.



The original motivation for this was to understand how the level of clustering in herd immunity would effect it's efficacy as a strategy.



Requires InvertedIndices and StatsBase packages
