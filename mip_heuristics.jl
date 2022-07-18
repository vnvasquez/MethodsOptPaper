using DataFrames
using CSV

###################################
# Post-processing functions (A & B)
###################################

function post_process_minval(
    release_size_min_per_timestep::Float64,
    control_vec::Vector{Float64},
    node_number::Int64,
)
    times = collect(1:length(control_vec))
    release = zeros(length(control_vec))

    for (i, v) in enumerate(control_vec)
        if v >= release_size_min_per_timestep
            release[i] = v
        end
    end

    node = fill(node_number, length(times))
    new_control_mat = hcat(Int64.(node), times, release)

    return new_control_mat
end

function post_process_intval(
    rel_vec::Vector{Float64},
    control_vec::Vector{Float64},
    node_number::Int64,
)
    times = collect(1:length(control_vec))
    release = zeros(length(control_vec))

    for (i, v) in enumerate(rel_vec)
        if isapprox(1.0, v; rtol=0.1) # NB 
            release[i] = control_vec[i]
        end
    end

    node = fill(node_number, length(times))
    new_control_mat = hcat(Int64.(node), times, round.(release; digits=3))

    return new_control_mat
end

###################################
# RIDL MIP
###################################

# Data 
ridl_binloc = CSV.read("mip_bonmin1_RIDL_spatial_releaselocation.csv", DataFrame)
ridl_cm4 = CSV.read("mip_bonmin1_RIDL_spatial_node4_controlM.csv", DataFrame)
ridl_cm5 = CSV.read("mip_bonmin1_RIDL_spatial_node5_controlM.csv", DataFrame);

# BINARY rule 
int_4 = post_process_intval(ridl_binloc.x4, ridl_cm4.control_M_G3, 4)
int_5 = post_process_intval(ridl_binloc.x5, ridl_cm5.control_M_G3, 5)

# MINIMUM rule
min_4 = post_process_minval(1000.0, ridl_cm4.control_M_G3, 4)
min_5 = post_process_minval(1000.0, ridl_cm5.control_M_G3, 5)

# Policy 
RIDL_int_policy = Dict("RIDL_int_n4" => int_4, "RIDL_int_n5" => int_5)
RIDL_min_policy = Dict("RIDL_min_n4" => min_4, "RIDL_min_n5" => min_5)

###################################
# WOLB MIP
###################################

# Data
wolb_binloc = CSV.read("mip_bonmin1_WOLB_spatial_releaselocation.csv", DataFrame)
wolb_cm1 = CSV.read("mip_bonmin1_WOLB_spatial_node1_controlM.csv", DataFrame)
wolb_cm2 = CSV.read("mip_bonmin1_WOLB_spatial_node2_controlM.csv", DataFrame)
wolb_cm3 = CSV.read("mip_bonmin1_WOLB_spatial_node3_controlM.csv", DataFrame)
wolb_cm4 = CSV.read("mip_bonmin1_WOLB_spatial_node4_controlM.csv", DataFrame)
wolb_cm5 = CSV.read("mip_bonmin1_WOLB_spatial_node5_controlM.csv", DataFrame);

# BINARY rule 
int_1m = post_process_intval(wolb_binloc.x1, wolb_cm1.control_M_G1, 1)
int_2m = post_process_intval(wolb_binloc.x2, wolb_cm2.control_M_G1, 2)
int_3m = post_process_intval(wolb_binloc.x3, wolb_cm3.control_M_G1, 3)
int_4m = post_process_intval(wolb_binloc.x4, wolb_cm4.control_M_G1, 4)
int_5m = post_process_intval(wolb_binloc.x5, wolb_cm5.control_M_G1, 5)

# MINIMUM rule 
min_1m = post_process_minval(1.0, wolb_cm1.control_M_G1, 1)
min_2m = post_process_minval(1.0, wolb_cm2.control_M_G1, 2)
min_3m = post_process_minval(1.0, wolb_cm3.control_M_G1, 3)
min_4m = post_process_minval(1.0, wolb_cm4.control_M_G1, 4)
min_5m = post_process_minval(1.0, wolb_cm5.control_M_G1, 5)

# Policy 
WOLB_int_policy_M = Dict(
    "WOLB_int_n1" => int_1m,
    "WOLB_int_n2" => int_2m,
    "WOLB_int_n3" => int_3m,
    "WOLB_int_n4" => int_4m,
    "WOLB_int_n5" => int_5m,
)
WOLB_min_policy_M = Dict(
    "WOLB_min_n1" => min_1m,
    "WOLB_min_n2" => min_2m,
    "WOLB_min_n3" => min_3m,
    "WOLB_min_n4" => min_4m,
    "WOLB_min_n5" => min_5m,
)
