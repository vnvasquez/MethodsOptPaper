
using GeneDrive
using Ipopt
using JuMP
using Juniper
using Cbc
using Gurobi
using AmplNLWriter
using Bonmin_jll
using DataStructures

include("objectives.jl");

###################################
# Data 
###################################

species = AedesAegypti
genetics = genetics_ridl();
enviro_response = stages_rossi();
update_population_size(enviro_response, 500);
organisms = make_organisms(species, genetics, enviro_response);
temperature = example_temperature_timeseries;
length(temperature.values)
coordinates = (16.0, 146.0);

node1 = Node(:Node1, organisms, temperature, coordinates);
node2 = Node(:Node2, organisms, temperature, coordinates);
node3 = Node(:Node3, organisms, temperature, coordinates);
node4 = Node(:Node4, organisms, temperature, coordinates);
node5 = Node(:Node5, organisms, temperature, coordinates);

ridlnet2 = Network(:RIDLNet, node1, node2, node3, node4, node5);
move_rate = 0.002

migration_data = Dict( # 1 <-> 2
    ("Male", "WW") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Male", "WR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Male", "RR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Female", "WW") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Female", "WR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Female", "RR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    # 1 <-> 3
    ("Male", "WW") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Male", "WR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Male", "RR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Female", "WW") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Female", "WR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Female", "RR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    # 3 <-> 2
    ("Male", "WW") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Male", "WR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Male", "RR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Female", "WW") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Female", "WR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Female", "RR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    # 3 <-> 4
    ("Male", "WW") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Male", "WR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Male", "RR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Female", "WW") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Female", "WR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Female", "RR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    # 4 <-> 5
    ("Male", "WW") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Male", "WR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Male", "RR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Female", "WW") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Female", "WR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Female", "RR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
);

assign_migration!(ridlnet, migration_data2, species);
release_gene = get_homozygous_modified(node1, species)
wild_gene = get_wildtype(node1, species)

tspan = (1, length(temperature.values))

op_constraints1 = ReleaseStrategy(
    release_this_gene_index=release_gene,
    release_this_life_stage=Male,
    release_time_interval=7,
    release_start_time=10,
    release_size_max_per_timestep=50000.0,
);

mystrategy = Dict(
    1 => op_constraints1,
    2 => op_constraints1,
    3 => op_constraints1,
    4 => op_constraints1,
    5 => op_constraints1,
);

###################################
# NLP 
###################################

i = JuMP.optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level" => 1,
    "linear_solver" => "pardiso",
);
#"linear_solver" =>  "ma86");

prob = GeneDrive.create_decision_model(
    ridlnet,
    tspan;
    node_strategy=mystrategy,
    species=species,
    optimizer=i,
    slack_small=true,
);

# TargetPercentageByDate:
sol_targperc = solve_decision_model(
    prob,
    TargetPercentageByDate;
    wildtype=wild_gene,
    percent_suppression=0.20,
    target_timestep=200,
)

res_targperc = format_decision_model_results(sol_targperc);

# SpatialReduction:
sol_spatial = solve_decision_model(
    prob,
    SpatialReduction;
    wildtype=wild_gene,
    node_to_suppress=4,
    percent_to_suppress=0.98,
    do_binary=false,
)

res_spatial = format_decision_model_results(sol_spatial)

###################################
# MIP 
###################################

op_constraints = ReleaseStrategy(
    release_this_gene_index=release_gene,
    release_this_life_stage=Male,
    release_time_interval=7,
    release_start_time=10,
    release_size_min_per_timestep=1000.0,
    release_size_max_per_timestep=50000.0,
);

mystrategy = Dict(
    1 => op_constraints,
    2 => op_constraints,
    3 => op_constraints,
    4 => op_constraints,
    5 => op_constraints,
);

i_mip = JuMP.optimizer_with_attributes(
    Juniper.Optimizer,
    "nl_solver" => i,
    #"mip_solver" => HiGHS.Optimizer,
    #"mip_solver" => JuMP.optimizer_with_attributes(Cbc.Optimizer))
    "mip_solver" => JuMP.optimizer_with_attributes(
        Gurobi.Optimizer,
        "OutputFlag" => 1,
        "NumericFocus" => 1,
    ),
    "log_levels" => [:Table, :Info, :Option],
    "feasibility_pump" => true,
    "feasibility_pump_time_limit" => 3600,
    "strong_branching_time_limit" => Inf,
    "branch_strategy" => :MostInfeasible,
    "feasibility_pump_tolerance_counter" => 1000,
    "num_resolve_nlp_feasibility_pump" => 1000,
);

prob_mip = GeneDrive.create_decision_model(
    ridlnet,
    tspan;
    node_strategy=mystrategy,
    species=AedesAegypti,
    do_binary=true,
    #optimizer = i_mip,
    optimizer=JuMP.optimizer_with_attributes(
        () -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe),
        "bonmin.nlp_solver" => "Ipopt",
        "linear_solver" => "ma86",
        "honor_original_bounds" => "yes",
    ),
    slack_small=true,
);

# SpatialReduction:
sol_spatial_mip = solve_decision_model_test(
    prob_mip,
    SpatialReduction;
    wildtype=wild_gene,
    node_to_suppress=4,
    percent_to_suppress=0.98,
    do_binary=true,
)

res_spatial_mip = format_decision_model_results(sol_spatial_mip)
