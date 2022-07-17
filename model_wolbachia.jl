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
genetics = genetics_wolbachia();
enviro_response = stages_rossi();
update_population_size(enviro_response, 500);
organisms = make_organisms(species, genetics, enviro_response);
temperature = example_temperature_timeseries
coordinates = (16.0, 146.0);

node1 = Node(:Node1, organisms,
    temperature, coordinates);
node2 = Node(:Node2, organisms,
    temperature, coordinates);
node3 = Node(:Node3, organisms,
    temperature, coordinates);
node4 = Node(:Node4, organisms,
    temperature, coordinates);
node5 = Node(:Node5, organisms,
    temperature, coordinates);

wolbnet = Network(:WOLBNet, node1, node2, node3, node4, node5);
move_rate = 0.002 

migration_data = Dict( # 1 <-> 2
("Male", "WW") => Dict((:Node1, :Node2) => move_rate,
                    (:Node2, :Node1) => move_rate),
("Male", "ww") => Dict((:Node1, :Node2) => move_rate,
                    (:Node2, :Node1) => move_rate),
("Female", "WW") => Dict((:Node1, :Node2) => move_rate,
                        (:Node2, :Node1) => move_rate),
("Female", "ww") => Dict((:Node1, :Node2) => move_rate,
                        (:Node2, :Node1) => move_rate),
                        # 1 <-> 3
("Male", "WW") => Dict((:Node1, :Node3) => move_rate,
                        (:Node3, :Node1) => move_rate),
("Male", "ww") => Dict((:Node1, :Node3) => move_rate,
                        (:Node3, :Node1) => move_rate),
("Female", "WW") => Dict((:Node1, :Node3) => move_rate,
                        (:Node3, :Node1) => move_rate),
("Female", "ww") => Dict((:Node1, :Node3) => move_rate,
                        (:Node3, :Node1) => move_rate),
                            # 3 <-> 2
("Male", "WW") => Dict((:Node2, :Node3) => move_rate,
                        (:Node3, :Node2) => move_rate),
("Male", "ww") => Dict((:Node2, :Node3) => move_rate,
                        (:Node3, :Node2) => move_rate),
("Female", "WW") => Dict((:Node2, :Node3) => move_rate,
                        (:Node3, :Node2) => move_rate),
("Female", "ww") => Dict((:Node2, :Node3) => move_rate,
                        (:Node3, :Node2) => move_rate),
                            # 3 <-> 4
("Male", "WW") => Dict((:Node4, :Node3) => move_rate,
                        (:Node3, :Node4) => move_rate),
("Male", "ww") => Dict((:Node4, :Node3) => move_rate,
                        (:Node3, :Node4) => move_rate),
("Female", "WW") => Dict((:Node4, :Node3) => move_rate,
                        (:Node3, :Node4) => move_rate),
("Female", "ww") => Dict((:Node4, :Node3) => move_rate,
                        (:Node3, :Node4) => move_rate),
                            # 4 <-> 5
("Male", "WW") => Dict((:Node4, :Node5) => move_rate,
                        (:Node5, :Node4) => move_rate),
("Male", "ww") => Dict((:Node4, :Node5) => move_rate,
                        (:Node5, :Node4) => move_rate),
("Female", "WW") => Dict((:Node4, :Node5) => move_rate,
                        (:Node5, :Node4) => move_rate),
("Female", "ww") => Dict((:Node4, :Node5) => move_rate,
                        (:Node5, :Node4) => move_rate),
);

assign_migration!(wolbnet2, migration_data, species);
release_gene = get_homozygous_modified(node1, species)
wild_gene = get_wildtype(node1, species)

tspan = (1, length(temperature.values))

op_constraints1 = ReleaseStrategy(
    release_this_gene_index=release_gene,
    release_this_life_stage= (Male,Female), 
    release_time_interval=7,
    release_start_time= 10,
    release_size_max_per_timestep= 50000.0,
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

i = JuMP.optimizer_with_attributes(Ipopt.Optimizer,
    #"linear_solver" =>  "pardiso");
    "linear_solver" =>  "ma86");

prob = GeneDrive.create_decision_model(
    wolbnet,
    tspan;
    node_strategy = mystrategy,
    species = species,
    optimizer = i,
    slack_small = true
    );

# TargetPercentageByDate:
sol_targperc = solve_decision_model_wolb(prob,
    TargetPercentageByDate;
    wildtype=wild_gene,
    releasetype=release_gene, 
    percent_suppression=.20,
    target_timestep= 200) 

res_targperc = format_decision_model_results(sol_targperc)

# SpatialReduction:
sol_spatial = solve_decision_model_wolb(prob,
    SpatialReduction;
    wildtype=wild_gene,
    node_to_suppress = 4,
    percent_to_suppress=.98,
    do_binary=false)
    
res_spatial = format_decision_model_results(sol_spatial)

###################################
# MIP 
###################################

op_constraints = ReleaseStrategy(
    release_this_gene_index=release_gene,
    release_this_life_stage= (Male,Female),  
    release_time_interval=7,
    release_start_time= 10,
    release_size_max_per_timestep= 50000.0,
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
    "mip_solver" => JuMP.optimizer_with_attributes(Gurobi.Optimizer,
        "NumericFocus" => 1)) 

prob_mip = GeneDrive.create_decision_model(
    wolbnet2,
    tspan;
    node_strategy = mystrategy3, # no forcing, yes min
    species = AedesAegypti,
    do_binary = true,
    #optimizer = i_mip,
    optimizer = JuMP.optimizer_with_attributes(
        () -> AmplNLWriter.Optimizer(Bonmin_jll.amplexe),
        "bonmin.nlp_solver" => "Ipopt",
        "linear_solver" => "ma86",
        "honor_original_bounds" => "yes",
        ),
    slack_small = true
    );

# SpatialReduction:
sol_spatial_mip = solve_decision_model_wolb(prob_mip,
    SpatialReduction;
    wildtype=wild_gene,
    node_to_suppress = 4,
    percent_to_suppress=.98,
    do_binary=true)

res_spatial_mip = format_decision_model_results(sol_spatial_mip)
