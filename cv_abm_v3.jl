



### defining functions

# finds which existing dictionaries have a reference to this one
function find_conns(net::Dict, iter::Int64)

idx = Int64[]
strength = rand(0)

for y in 1:(iter - 1)
    if iter in net["$y"]["conn"]
        println(y)
        bt = [net["$y"]["conn"] .== iter][1]
        push!(strength, net["$y"]["strength"][bt][1])  # not subsetting on iter 331
        push!(idx, y)
    end
end

return idx, strength

end # close function






# add random connections from existing dunbar1 connections
function durbar1_conn(net::Dict, next_vals::Vector{Int64},
    next_strength::Vector{Float64}, dweight1::Float64,
    d1_conns_to_add::Int64)

d1_conns_to_add = trunc(Int64, (150 - length(next_vals)) * proportion_shared)

idx = next_strength .== dweight1

all_dunbar1_conns = next_vals[idx]

if length(all_dunbar1_conns) < 1
    return next_vals, next_strength
else

for i in 1:length(all_dunbar1_conns)
k = all_dunbar1_conns[i]

all_conn = Int64[]
all_strength = Float64[]

all_conn = vcat(all_conn, net["$k"]["conn"])
all_strength = vcat(all_strength, net["$k"]["strength"])
end

# random sample
idx = unique(wsample(1:length(all_strength), all_strength, d1_conns_to_add))


return vcat(next_vals, all_conn[idx]), vcat(next_strength, all_strength[idx])
end # end of if statement

end # end of function






# function to count each occurance
function counter(next_strength::Array{Float64})

    l1 = 0
    l2 = 0
    l3 = 0
    l4 = 0

    for i in 1:length(next_strength)
        if next_strength[i] == dweight1
            l1 += 1
        elseif next_strength[i] == dweight2
            l2 += 1
        elseif next_strength[i] == dweight3
            l3 += 1
        elseif next_strength[i] == dweight4
            l3 += 1
        end
    end
return l1, l2, l3, l4
end    # end of function




# function to return index of first K rows which are a certain value in an array
function index_k(inmain::Array{Float64}, key::Float64, K::Int64)

    out = Int64[]
    counter_val = 0

    for i in 1:length(inmain)
        if inmain[i] == key

        counter_val = counter_val + 1
        push!(out, i)
        end

        if counter_val == K
            break
        end


    end

    return out

end





using InvertedIndices
using Distributions
using StatsBase

net_size = 1000
r_value = 2.5
days_infectious = 6.0

# agent based modelling for spread of CV19, and extent to which 60% immunity
# prevents more spread, depending on the network characteristics


# setting parameters
dunbar1 = 5
dunbar2 = 15
dunbar3 = 40
dunbar4 = 90
dweight1 = 10.0
dweight2 = 4.0
dweight3 = 2.0
dweight4 = 1.0

daily_spread_base =  r_value /
            (((dunbar1 * dweight1) + (dunbar2 * dweight2) +
            (dunbar3 * dweight3) + (dunbar4 * dweight4)) *
            days_infectious)

# how many connections will the average person share with another?
proportion_shared = 0.75



# initialising network with 150 random connections
net = Dict()
net["1"] = Dict()
net["1"]["conn"] = rand(1:net_size, 150)
net["1"]["strength"] = vcat(fill(dweight1, dunbar1),
                fill(dweight2, dunbar2),
                fill(dweight3, dunbar3),
                fill(dweight4, dunbar4))
net["1"]["status"] = "healthy"


i = 1
println(length(net["$i"]["strength"]))



# Main loop to look through all dicts to find which reference it
for iter in 2:net_size

next_vals, next_strength = find_conns(net, iter)



# at this point skip iteration if the network size is over 150
if length(next_vals) < 150

println(length(next_vals))
next_vals, next_strength = durbar1_conn(net, next_vals, next_strength,
            dweight1, d1_conns_to_add)

println(length(next_vals))
end




# removing dodge values
if length(next_strength) > 1
l1, l2, l3, l4 = counter(next_strength)
l1 = max(l1 - dunbar1, 0)
l2 = max(l2 - dunbar2, 0)
l3 = max(l3 - dunbar3, 0)
l4 = max(l4 - dunbar4, 0)
else
l1, l2, l3, l4 = 0, 0, 0, 0
end


if l1 > 0
idx1 = index_k(next_strength, dweight1, l1)
else
idx1 = Int64[]
end

if l2 > 0
idx2 = index_k(next_strength, dweight2, l2)
else
idx2 = Int64[]
end

if l3 > 0
idx3 = index_k(next_strength, dweight3, l3)
else
idx3 = Int64[]
end

if l4 > 0
idx4 = index_k(next_strength, dweight4, l4)
else
idx4 = Int64[]
end

idx = vcat(idx1, idx2, idx3, idx4)

# negated subsetting
next_vals = next_vals[Not(idx)]
next_strength = next_strength[Not(idx)]














# another jump point: skip iteration if the network size is over 150
if length(next_vals) < 150

# create random connections to get to 150 total
final_vals = vcat(next_vals, rand(1:net_size, 150 - length(next_vals)))


# setting strengths of connections
d1 = fill(dweight1, max(dunbar1 - sum(next_strength .== dweight1), 0))
d2 = fill(dweight2, max(dunbar2 - sum(next_strength .== dweight2), 0))
d3 = fill(dweight3, max(dunbar3 - sum(next_strength .== dweight3), 0))
d4 = fill(dweight4, max(dunbar4 - sum(next_strength .== dweight4), 0))


final_strength = vcat(next_strength, d1, d2, d3, d4)
end





# saving in dictionary
net["$iter"] = Dict()
net["$iter"]["conn"] = final_vals
net["$iter"]["strength"] = final_strength
net["$iter"]["status"] = "healthy"


# ending loop
println("Connection no. $iter grown")
end





for iter=1:length(net)
println(countmap(net["$iter"]["strength"]))
end

for iter=1:length(net)
println(countmap(net["$iter"]["conn"]))
end


for iter=1:length(net)
println(length(net["$iter"]["strength"]))
println(length(net["$iter"]["conn"]))
end



function ia(net::Dict)
x = []
for i=1:length(net)
x = vcat(x, net["$i"]["conn"])
end
return x
end
d = ia(net)
println(countmap(d))

#### the one outstanding problem is getting strength added: it's not
#### working somewhere along the way






# will having yourself listed in your network break the code? Who can say







#### Spreading the virus: function to update someone's status


# each timestep can represent a day, with the likelihood of passing it
# on determined by (2.5 / no. of days one is infected for on average), with
# some perturbation


# who is immune can be set at random at the start, or for clustering can
# "spread" immunity using the same contagion model as the virus, starting
# with one person and running the spread until 60% of the population is
# immune (this would represent the most clustered version). A halfway
# point between these two would be to pick multiple random points in the
# population from which to "spread immunity" (at the other extreme, starting
# with 60% of the population to spread immunity from would stop the modelling
# right away, and therefore be the same as assigning immunity randomly)
