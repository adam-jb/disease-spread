
### setting parameters for network creation
net_size = 100000

# size and weight of infection chance for each level of someone's network
dunbar1 = 5
dunbar2 = 15
dunbar3 = 40
dunbar4 = 90
dweight1 = 10.0
dweight2 = 4.0
dweight3 = 2.0
dweight4 = 1.0

# how many connections will the average person share with all their
# closest (dunbar1) connections?
proportion_shared = 0.50



# setting parameters for disease spread
cluster_strength = [0.1]
initial_cases = [50]
days_vals = [6]
population_start_immune = [0.6]
r_values = [3.2]
iters = 50
cluster_link_multiplier = [0.05, 0.3, 1.5]


# initialise network
net = initialise_network(net_size = net_size,
                        dunbar1 = dunbar1,
                        dunbar2 = dunbar2,
                        dunbar3 = dunbar3,
                        dunbar4 = dunbar4,
                        dweight1 = dweight1,
                        dweight2 = dweight2,
                        dweight3 = dweight3,
                        dweight4 = dweight4,
                        proportion_shared = proportion_shared,
                        verbose = true)



# running simulation
output, days_mat = sim_spread_mult(net = net,
                cluster_strength = cluster_strength,
                initial_cases = initial_cases,
                days_vals = days_vals,
                population_start_immune = population_start_immune,
                iters = iters,
                r_values = r_values,
                dunbar1 = dunbar1,
                dunbar2 = dunbar2,
                dunbar3 = dunbar3,
                dunbar4 = dunbar4,
                dweight1 = dweight1,
                dweight2 = dweight2,
                dweight3 = dweight3,
                dweight4 = dweight4,
                verbose = true,
                cluster_link_multiplier = cluster_link_multiplier,
                return_days = true)
