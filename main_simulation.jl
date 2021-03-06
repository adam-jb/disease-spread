

# agent based modelling for spread of a contagion, and extent to which a certain
# level of existing immunity prevents more spread, depending on the
# network characteristics and disease parameters

# could add some perturbation in the network as it's very uniform




using InvertedIndices
using StatsBase



### Functions to create network

"finds which people already have a connection to this one"
function find_conns(net::Dict, iter::Int64)

idx = Int64[]
strength = rand(0)

for y in 1:(iter - 1)
    if iter in net["$y"]["conn"]
        # println(y)
        bt = [net["$y"]["conn"] .== iter][1]
        push!(strength, net["$y"]["strength"][bt][1])
        push!(idx, y)
    end
end

return idx, strength

end # close function




"removing any new connections to people with max connections already"
function looker(vals_out::Array{Int64}, strengths_out::Array{Float64},
    mentions::Dict, total_dun::Int64)
idx = Bool[]
len = length(vals_out)
counting = 0
m = 1
while counting < len
    if mentions[vals_out[m]] > (total_dun - 1)
        deleteat!(vals_out, m)
        deleteat!(strengths_out, m)
        counting = counting + 1
    else m = m + 1
        counting = counting + 1
    end
end
return vals_out, strengths_out

end # close function



"removing any new connections to people with max connections already"
function looker_val_only(vals_out::Array{Int64},
    mentions::Dict, total_dun::Int64)
idx = Bool[]
len = length(vals_out)
counting = 0
m = 1
while counting < len
    if mentions[vals_out[m]] > (total_dun - 1)
        deleteat!(vals_out, m)
        counting = counting + 1
    else m = m + 1
        counting = counting + 1
    end
end
return vals_out

end # close function




"count frequency of each dunbar level"
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




"returns index of first K values which have a specified (key) value"
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





"concatenate all dunbar1 connection's connections into an array"
function concat_opts(net::Dict, all_dunbar1_conns::Array)


all_conn = Int64[]
all_strength = Float64[]

for i in 1:length(all_dunbar1_conns)
k = all_dunbar1_conns[i]
all_conn = vcat(all_conn, net["$k"]["conn"])
all_strength = vcat(all_strength, net["$k"]["strength"])
end

return all_conn, all_strength

end # close function





"add random connections from existing dunbar1 connections"
function durbar1_conn(net::Dict, next_vals::Vector{Int64},
    next_strength::Vector{Float64}, dweight1::Float64,
    mentions::Dict, total_dun::Int64)

idx = next_strength .== dweight1

all_dunbar1_conns = next_vals[idx]

d1_conns_to_add = trunc(Int64,
    ((total_dun - length(next_vals)) * proportion_shared * length(all_dunbar1_conns)) /
        5)

if length(all_dunbar1_conns) < 1
    return next_vals, next_strength
else


# concatenating all connections of dunbar1 friends
all_conn, all_strength = concat_opts(net, all_dunbar1_conns)


# remove connections with people who already have threshold no. of connections
all_conn, all_strength = looker(all_conn, all_strength, mentions, total_dun)



# weighted random sample of these
idx = unique(wsample(1:length(all_strength), all_strength, d1_conns_to_add))
all_conn = all_conn[idx]
all_strength = all_strength[idx]



## removing new connections where the dunbar layer count is too high
l1, l2, l3, l4 = counter(next_strength)
l1 = dunbar1 - l1
l2 = dunbar2 - l2
l3 = dunbar3 - l3
l4 = dunbar4 - l4

n1, n2, n3, n4 = counter(all_strength)
n1 = n1 - l1
n2 = n2 - l2
n3 = n3 - l3
n4 = n4 - l4

if n1 > 0
idx1 = index_k(all_strength, dweight1, n1)
else
idx1 = Int64[]
end

if n2 > 0
idx2 = index_k(all_strength, dweight2, n2)
else
idx2 = Int64[]
end

if n3 > 0
idx3 = index_k(all_strength, dweight3, n3)
else
idx3 = Int64[]
end

if n4 > 0
idx4 = index_k(all_strength, dweight4, n4)
else
idx4 = Int64[]
end

idx = vcat(idx1, idx2, idx3, idx4)

# negated subsetting to remove the excess layers
all_conn = all_conn[Not(idx)]
all_strength = all_strength[Not(idx)]



# combining results
vals_out = vcat(next_vals, all_conn)
strength_out = vcat(next_strength, all_strength)


# removing duplicate connections
idx = indexin(unique(vals_out), vals_out)
vals_out = vals_out[idx]
strength_out = strength_out[idx]



return vals_out, strength_out
end # end of if statement

end # end of function





"creating new connections to people not yet considered"
function find_forward_newvals(next_vals::Array{Int64}, total_dun::Int64, iter::Int64,
    net_size::Int64, mentions::Dict)
countinga = 0
while (length(next_vals) < (total_dun)) & (countinga < 5)

max_new = net_size - (iter + 1)

cand_vals = sample((iter + 1):net_size,
                                min(total_dun - length(next_vals), max_new),
                                replace = false)

next_vals = unique(vcat(next_vals, looker_val_only(cand_vals, mentions, total_dun)))
countinga += 1

end  # close while loop

return next_vals
end # close function






#### diagnostics of network creation

# checking that all folks get their 150 connections from both directions
function ia(net::Dict)
x = []
for i=1:length(net)
x = vcat(x, net["$i"]["conn"])
end
return x
end

function histogram(s)
    d = Dict()
    for c in s
        if c ∉ keys(d)
            d[c] = 1
        else
            d[c] += 1
        end
    end
    d
end


# finding where the lengths of any vectors aren't equal to no. of desired connections
# that the results get the same columns as println(ih) above, shows the connections
# are missing in both directions
function tester(net::Dict)
out1 = []
out2 = []
for iter=1:length(net)
out1 = push!(out1, length(net["$iter"]["strength"]))
out2 = push!(out2, length(net["$iter"]["conn"]))
end
return out1, out2
end


"checking that every connection is logged for both parties"
function check_net(;net::Dict, total_dun::Int64)
dc = histogram(ia(net))
ih = Int64[]
for i in 1:length(net)
    if dc[i] != total_dun
        push!(ih, i)
    end
end

o1, o2 = tester(net)
ix = findall(map(x -> x != total_dun, o1))

g = findall(ih != ix)

# whether these are the same checks that every connection is reciprocated
    if ih == ix
    println("All connections in network reciprocated")
    else println("Error in network creations: some connections are only
                logged for one party. Discrepency in these rows: $g")
    end
end





function initialise_network(;net_size::Int64 = 500,
    proportion_shared = 0.75,
    dunbar1::Int64 = 5,
    dunbar2::Int64 = 15,
    dunbar3::Int64 = 40,
    dunbar4::Int64 = 90,
    dweight1::Float64 = 10.0,
    dweight2::Float64 = 4.0,
    dweight3::Float64 = 2.0,
    dweight4::Float64 = 1.0,
    verbose::Bool = true,
    check_network::Bool = false)


total_dun = dunbar1 + dunbar2 + dunbar3 + dunbar4

# initialising first member of network with 150 random connections
net = Dict()
for i in 1:net_size
net["$i"] = Dict()
net["$i"]["conn"] = Int64[]
net["$i"]["strength"] = zeros(0)
net["$i"]["status"] = "healthy"
end


net["1"]["conn"] = sample(2:net_size, total_dun, replace = false)
net["1"]["strength"] = vcat(fill(dweight1, dunbar1),
                fill(dweight2, dunbar2),
                fill(dweight3, dunbar3),
                fill(dweight4, dunbar4))


# initalising record of how many times a value is mentioned
mentions = Dict()
for i in 1:net_size
    mentions[i] = 0
end

# capturing mentions for first person's network
for m in 1:length(net["1"]["conn"])
    mentions[net["1"]["conn"][m]] += 1
end
mentions[1] = total_dun

# pushing values to other people's dictionaries
for m in 1:length(net["1"]["conn"])
    idx = net["1"]["conn"][m]
    push!(net["$idx"]["conn"], 1)
    push!(net["$idx"]["strength"], net["1"]["strength"][m])
end





### Main loop to create network
for iter in 2:net_size

# finding existing connections, and hold these as a separate
# vector as they won't need adding to mentions

if iter == net_size
    continue
end



next_vals = net["$iter"]["conn"]
next_strength = net["$iter"]["strength"]


score0 = length(next_vals)




# creating connections from existing dunbar1 level connections
if length(next_vals) < total_dun

next_vals, next_strength = durbar1_conn(net, next_vals, next_strength,
            dweight1, mentions, total_dun)
end


store1 = length(next_vals)




# removing any self-reference
if sum(next_vals .== iter) > 0
    next_vals = next_vals[Not(next_vals .== iter)]
end


# creating random connections to bring to maximum
if length(next_vals) < total_dun
final_vals = find_forward_newvals(next_vals, total_dun, iter, net_size, mentions)
else
final_vals = next_vals
end



# setting strengths of new connections
d1 = fill(dweight1, max(dunbar1 - sum(next_strength .== dweight1), 0))
d2 = fill(dweight2, max(dunbar2 - sum(next_strength .== dweight2), 0))
d3 = fill(dweight3, max(dunbar3 - sum(next_strength .== dweight3), 0))
d4 = fill(dweight4, max(dunbar4 - sum(next_strength .== dweight4), 0))

final_strength = vcat(next_strength, d1, d2, d3, d4)


store2 = length(final_vals)


# removing excess strength values
if length(final_strength) > length(final_vals)
    final_strength = final_strength[1:length(final_vals)]
end

if length(final_strength) < length(final_vals)
    final_strength = vcat(final_strength,
        fill(dweight4, length(final_vals) - length(final_strength)))
end



store3 = length(final_strength)



# determining which connections are new, dropping the existing onea
function find_new_idx(net::Dict, final_vals::Array{Int64})
idx = Int64[]
for m in 1:length(final_vals)
    if final_vals[m] ∉ net["$iter"]["conn"]
        push!(idx, m)
    end
end
return idx
end # close function

idx = find_new_idx(net, final_vals)
final_vals = final_vals[idx]
final_strength = final_strength[idx]


# pushing new values to other people's dictionaries
for m in 1:length(final_vals)
    idx = final_vals[m]
    push!(net["$idx"]["conn"], iter)
    push!(net["$idx"]["strength"], final_strength[m])
end


# capturing mentions of upcoming connections
final_vals = final_vals[final_vals .> iter]
for m in 1:length(final_vals)
    mentions[final_vals[m]] += 1
end

# saving new  connections in dictionary
net["$iter"]["conn"] = vcat(net["$iter"]["conn"], final_vals)
net["$iter"]["strength"] = vcat(net["$iter"]["strength"], final_strength)


# updating this person's mentions
mentions[iter] = total_dun




if verbose
println("Connection no. $iter grown: $score0, $store1, $store2, $store3")
end

# ending loop
end

println(string("Network created: object size = ",
            trunc(Int64, Base.summarysize(net)/10^6), "Mb; has $net_size people"))


if check_network == true
 check_net(net = net, total_dun = total_dun)
end

return net

end # close function





## only need to make a network once: can reset the statuses of all members
## to healthy for each run of the disease spread



##### now the network is made: immunity and disease spread section

"Summary stats on the status of the population in a network"
function get_status(net::Dict)

        d = Dict()
    for i in 1:length(net)
        if net["$i"]["status"] ∉ keys(d)
            d[net["$i"]["status"]] = 1
        else
            d[net["$i"]["status"]] += 1
        end
    end
    return d
end # closing function



"Resetting everyone's status to healthy"
function reset_net(net::Dict)
for i in 1:length(net)
    net["$i"]["status"] = "healthy"
end
return net
end




"initialising immuniity to be spread from init_k peaple"
function immune_init(net::Dict, init_k::Int64)
start = sample(1:length(net), init_k, replace = false)
for ref in start
net["$ref"]["status"] = "already immune"
end
return net
end




"Spreading immunity"
function immunity_spread(net::Dict, limit::Float64, net_size::Int64,
                        daily_spread_base::Float64, imspread_mult::Float64)

daily_spread_base = daily_spread_base * imspread_mult

while get_status(net)["already immune"] < limit

    i = 0
    while i < net_size
        i += 1

        if net["$i"]["status"] == "already immune"

            idx = (net["$i"]["strength"] * daily_spread_base) .>
                        rand(length(net["$i"]["strength"]))
            idx2 = net["$i"]["conn"][idx]

            if length(idx2) .> 0
                for z in 1:length(idx2)
                    id = idx2[z]
                    if net["$id"]["status"] == "healthy"
                        net["$id"]["status"] = "already immune"
                    end
                end
            end

        end

    end # inner while

end # close while loop


# pruning to ensure exactly $limit are made immune
nprune = trunc(Int64, get_status(net)["already immune"] - limit)
if nprune > 0
    immuni = Int64[]
    for i in 1:length(net)
        if net["$i"]["status"] == "already immune"
            push!(immuni, i)
        end
    end

    idx = sample(immuni, nprune, replace = false)

    for i in 1:length(idx)
        id = idx[i]
        net["$id"]["status"] = "healthy"
    end
end # close prune if statement


return net
end # close function




"initialising illness with K randomly chosen peaple"
function illness_init(net::Dict, k_illness::Int64)

for i in 1:k_illness

    found_one = 0
    while found_one < 1
        ref = sample(1:length(net))
        if net["$ref"]["status"] == "healthy"
        net["$ref"]["status"] = 0
            found_one += 1
        end
    end
end

return net
end # closing function




"Seeing there's at least 1 infectious person in the population"
function current_cases(days_infectious::Int64, stat_temp::Dict)
a = collect(0:days_infectious)
b = collect(keys(stat_temp))
counter = 0
for a in a
    for b in b
        if string(a) == string(b)
            counter += 1
        end
    end
end

if counter == 0
    return false
else
    return true
end

end # close function





"Simulating the spread of the disease"
function disease_spread(net::Dict, max_days::Int64, daily_spread_base::Float64,
    days_infectious::Int64, verbose::Bool, return_days::Bool = false,
    init_imu::Float64 = 0.6)

stat_temp = Dict("1" => 0)
days = 0
infect_store = Int64[]

while days < max_days && current_cases(8, stat_temp)

days += 1

today_infections = 0

for i in 1:length(net)

    if net["$i"]["status"] in 0:days_infectious-1

        idx = (net["$i"]["strength"] * daily_spread_base) .>
                    rand(length(net["$i"]["strength"]))
        idx2 = net["$i"]["conn"][idx]
        if length(idx2) .> 0
            for i in 1:length(idx2)
                id = idx2[i]
                if net["$id"]["status"] == "healthy"
                    net["$id"]["status"] = 0
                end

                today_infections += 1

            end
        end

    net["$i"]["status"] = net["$i"]["status"] += 1

elseif net["$i"]["status"] == days_infectious
        net["$i"]["status"] = "had disease"
end

end

stat_temp = get_status(net)

push!(infect_store, today_infections)

end  # close while loop

if verbose
println(string("Simulation resolved in $days days"))
end


result = get_status(net)
outs = result["had disease"] / (net_size * (1 - init_imu))


if (return_days == false)
    return outs, days
else
    return outs, days, infect_store
end


end # close function



"Modelling spread for a given set of parameters"
function sim_spread(;net::Dict,
                    r_value::Float64 = 3.2,
                    init_imu::Float64 = 0.6,
                    cluspa::Float64 = 0.1,
                    k_illness::Int64 = 50,
                    days_infectious::Int64 = 6,
                    dunbar1::Int64 = 5,
                    dunbar2::Int64 = 15,
                    dunbar3::Int64 = 40,
                    dunbar4::Int64 = 90,
                    dweight1::Float64 = 10.0,
                    dweight2::Float64 = 4.0,
                    dweight3::Float64 = 2.0,
                    dweight4::Float64 = 1.0,
                    max_days::Int64 = 1000,
                    verbose::Bool = true,
                    cluster_link_multiplier::Float64 = 1.0,
                    return_days::Bool =false)


# this is proportion of someone infecting someone in the dunbar4 level of their
# network in a single given day
daily_spread_base =  r_value /
    (((dunbar1 * dweight1) + (dunbar2 * dweight2) +
    (dunbar3 * dweight3) + (dunbar4 * dweight4)) *
    days_infectious)

dunbar_total = dunbar1 + dunbar2 + dunbar3 + dunbar4

net_size = length(net)

limit = net_size * init_imu
init_k = max(trunc(Int64, cluspa * limit), 1)


# resetting everyone to healthy
net = reset_net(net)

# classifying immunity
a = time()
net = immune_init(net, init_k)
imspread_mult = (dunbar_total / r_value) * 0.5 * cluster_link_multiplier
net = immunity_spread(net, limit, net_size, daily_spread_base, imspread_mult)
b = time() - a
#println("Classified who's already immune in $b seconds")

# spreading disease
a = time()
net = illness_init(net, k_illness)

if return_days == false
    outs, days_len = disease_spread(net, max_days, daily_spread_base,
                                    days_infectious,
                                    verbose, return_days, init_imu)
else
    outs, days_len, days_vec = disease_spread(net, max_days, daily_spread_base,
                                        days_infectious,
                                        verbose, return_days, init_imu)
    days_vec = vcat(days_vec, fill(0, max_days - length(days_vec)))
end



b = time() - a
#println("Disease spread simulated in $b seconds")

# println(get_status(net))
if return_days == false
    return outs, (days_len * 1.0)
else
    return outs, (days_len * 1.0), days_vec
end

end # end of function




"Simulating many times over for each combination of inputs"
function sim_spread_mult(;net::Dict,
    cluster_strength::Array{Float64},
    initial_cases::Array{Int64},
    days_vals::Array{Int64},
    population_start_immune::Array{Float64},
    iters::Int64,
    r_values::Array{Float64},
    dunbar1::Int64 = 5,
    dunbar2::Int64 = 15,
    dunbar3::Int64 = 40,
    dunbar4::Int64 = 90,
    dweight1::Float64 = 10.0,
    dweight2::Float64 = 4.0,
    dweight3::Float64 = 2.0,
    dweight4::Float64 = 1.0,
    max_days::Int64 = 1000,
    verbose::Bool = true,
    cluster_link_multiplier::Array{Float64} = 1.0,
    return_days::Bool =false)



holder = reshape(Float64[], 0, 8)

days_holder = reshape(Int64[], 0, max_days)

    for i in cluster_strength
        for j in initial_cases
            for k in days_vals
                for l in population_start_immune
                    for m in r_values
                        for a in cluster_link_multiplier
                            for z in 1:iters

                                if return_days == false
                                g, r = sim_spread(net = net,
                                cluspa = i,
                                r_value = m,
                                init_imu = l,
                                k_illness = j,
                                days_infectious = k,
                                dunbar1 = dunbar1,
                                dunbar2 = dunbar2,
                                dunbar3 = dunbar3,
                                dunbar4 = dunbar4,
                                dweight1 = dweight1,
                                dweight2 = dweight2,
                                dweight3 = dweight3,
                                dweight4 = dweight4,
                                max_days = max_days,
                                verbose = verbose,
                                cluster_link_multiplier = a)

                                holder = vcat(holder,
                                    hcat(i, m, l, j, k, a, g, r))


                                else
                                    g, r, days_vec = sim_spread(net = net,
                                    cluspa = i,
                                    r_value = m,
                                    init_imu = l,
                                    k_illness = j,
                                    days_infectious = k,
                                    dunbar1 = dunbar1,
                                    dunbar2 = dunbar2,
                                    dunbar3 = dunbar3,
                                    dunbar4 = dunbar4,
                                    dweight1 = dweight1,
                                    dweight2 = dweight2,
                                    dweight3 = dweight3,
                                    dweight4 = dweight4,
                                    max_days = max_days,
                                    verbose = verbose,
                                    cluster_link_multiplier = a,
                                    return_days = true)

                                    holder = vcat(holder,
                                        hcat(i, m, l, j, k, a, g, r))

                                    days_holder = vcat(days_holder, days_vec')
                                end # end of if/else for days return

                            end
                        end
                    end
                end
            end
        end
    end

if return_days == false
    return holder
else
    return holder, days_holder
end

end
