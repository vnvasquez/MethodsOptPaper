
using GeneDrive
using Ipopt
using JuMP
using Juniper
using Cbc
using Gurobi
using DataStructures

include("objectives.jl");

###################################
# Data 
###################################

species = AedesAegypti
genetics = genetics_mcr();
enviro_response = stages_rossi();
update_population_size(enviro_response, 500);
organisms = make_organisms(species, genetics, enviro_response);
temperature = example_temperature_timeseries;
coordinates = (16.0, 146.0);

node1 = Node(:Node1, organisms, temperature, coordinates);
node2 = Node(:Node2, organisms, temperature, coordinates);
node3 = Node(:Node3, organisms, temperature, coordinates);
node4 = Node(:Node4, organisms, temperature, coordinates);
node5 = Node(:Node5, organisms, temperature, coordinates);

mcrnet = Network(:MCRNet, node1, node2, node3, node4, node5);

move_rate = 0.002

migration_data = Dict( # 1 <-> 2
    ("Male", "HH") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Male", "Hh") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Male", "HR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Male", "hh") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Male", "hR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Male", "RR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Female", "HH") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Female", "Hh") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Female", "HR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Female", "hh") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Female", "hR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    ("Female", "RR") =>
        Dict((:Node1, :Node2) => move_rate, (:Node2, :Node1) => move_rate),
    # 1 <-> 3
    ("Male", "HH") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Male", "Hh") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Male", "HR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Male", "hh") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Male", "hR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Male", "RR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Female", "HH") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Female", "Hh") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Female", "HR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Female", "hh") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Female", "hR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    ("Female", "RR") =>
        Dict((:Node1, :Node3) => move_rate, (:Node3, :Node1) => move_rate),
    # 3 <-> 2
    ("Male", "HH") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Male", "Hh") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Male", "HR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Male", "hh") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Male", "hR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Male", "RR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Female", "HH") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Female", "Hh") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Female", "HR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Female", "hh") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Female", "hR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    ("Female", "RR") =>
        Dict((:Node2, :Node3) => move_rate, (:Node3, :Node2) => move_rate),
    # 3 <-> 4
    ("Male", "HH") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Male", "Hh") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Male", "HR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Male", "hh") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Male", "hR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Male", "RR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Female", "HH") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Female", "Hh") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Female", "HR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Female", "hh") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Female", "hR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    ("Female", "RR") =>
        Dict((:Node4, :Node3) => move_rate, (:Node3, :Node4) => move_rate),
    # 4 <-> 5
    ("Male", "HH") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Male", "Hh") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Male", "HR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Male", "hh") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Male", "hR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Male", "RR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Female", "HH") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Female", "Hh") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Female", "HR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Female", "hh") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Female", "hR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
    ("Female", "RR") =>
        Dict((:Node4, :Node5) => move_rate, (:Node5, :Node4) => move_rate),
);

assign_migration!(mcrnet, migration_data, species);
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
    "print_level" => 3,
    #"linear_solver" =>  "pardiso");
    "linear_solver" => "ma86",
);

prob = GeneDrive.create_decision_model(
    mcrnet,
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

res_targperc = GeneDrive.format_decision_model_results(sol_targperc)

# SpatialReduction:
sol_spatial = solve_decision_model(
    prob,
    SpatialReduction;
    wildtype=wild_gene,
    node_to_suppress=4,
    percent_to_suppress=0.98,
)

res_spatial = GeneDrive.format_decision_model_results(sol_spatial)
