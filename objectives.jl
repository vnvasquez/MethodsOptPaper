abstract type ObjectiveFunction end
struct TargetPercentageByDate <: ObjectiveFunction end
struct SpatialReduction <: ObjectiveFunction end

###################################
# Objective 2.2
###################################

function solve_decision_model(
    model::JuMP.Model,
    objective_function::Type{<:TargetPercentageByDate};
    wildtype=nothing,
    percent_suppression=nothing,
    target_timestep=nothing,
)
    if isa(wildtype, Nothing) ||
       isa(percent_suppression, Nothing) ||
       isa(target_timestep, Nothing) #|| isa(which_node, Nothing)
        @warn(
            "Model not solved because the objective function is missing information. Values must be supplied for all keyword arguments: `wildtype` (Int64), `percent_suppression (Float64), and `target_timestep` (Int64)."
        )
    else
        # Variables
        control_M = model[:control_M]
        F = model[:F]
        control_F = model[:control_F]
        JuMP.fix.(control_F, 0.0; force=true)

        # Sets
        sets = model.obj_dict[:Sets]
        N = sets[:N]
        O = sets[:O]
        SF = sets[:SF]
        G = sets[:G]
        T = sets[:T]

        # Expression
        JuMP.@objective(
            model,
            Min,
            (
                sum(
                    (
                        F[n, o, s, wildtype, t] -
                        percent_suppression * F[n, o, s, wildtype, 1]
                    ) .^ 2 for n in N, o in O, s in SF, t in T[target_timestep:end]
                ) + 1e-8 * sum(control_M)
            )
        )

        # Run
        JuMP.optimize!(model)

        # Output
        if termination_status(model) != OPTIMAL
            @info("Termination status: $(termination_status(model))")
            return model
        else
            println("Objective value:", objective_value(model))
        end

        return model
    end
end

function solve_decision_model_wolb(
    model::JuMP.Model,
    objective_function::Type{<:TargetPercentageByDate};
    wildtype=nothing,
    releasetype=release_gene,
    percent_suppression=nothing,
    target_timestep=nothing,
)
    if isa(wildtype, Nothing) ||
       isa(percent_suppression, Nothing) ||
       isa(target_timestep, Nothing) #|| isa(which_node, Nothing)
        @warn(
            "Model not solved because the objective function is missing information. Values must be supplied for all keyword arguments: `wildtype` (Int64), `percent_suppression (Float64), and `target_timestep` (Int64)."
        )
    else

        # Variables
        control_M = model[:control_M]
        F = model[:F]
        control_F = model[:control_F]

        # Sets
        sets = model.obj_dict[:Sets]
        N = sets[:N]
        O = sets[:O]
        SF = sets[:SF]
        G = sets[:G]
        T = sets[:T]

        # Expression
        JuMP.@objective(
            model,
            Min,
            (
                sum(
                    (
                        F[n, o, s, wildtype, t] -
                        percent_suppression * F[n, o, s, wildtype, 1]
                    ) .^ 2 for n in N, o in O, s in SF, t in T[target_timestep:end]
                ) +
                1e-8 * sum(control_M) +
                1e-8 * sum(control_F)
            )
        )

        # Run
        JuMP.optimize!(model)

        # Output
        if termination_status(model) != OPTIMAL
            @info("Termination status: $(termination_status(model))")
            return model
        else
            println("Objective value:", objective_value(model))
        end
        return model
    end
end

###################################
# Objective 2.3
###################################

function solve_decision_model(
    model::JuMP.Model,
    objective_function::Type{<:SpatialReduction};
    wildtype=nothing,
    homozygous_modified=nothing,
    node_to_suppress=nothing,
    percent_to_suppress=nothing,
    do_binary=false,
    policy=nothing,
)

    # Variables
    control_M = model[:control_M]
    F = model[:F]
    control_F = model[:control_F]
    JuMP.fix.(control_F, 0.0; force=true)

    # Sets
    sets = model.obj_dict[:Sets]
    N = sets[:N]
    O = sets[:O]
    SF = sets[:SF]
    SM = sets[:SM]
    G = sets[:G]
    T = sets[:T]

    # Expression
    costquant = 1e-8
    costloc = 1e-8
    release_location = model[:release_location]

    if do_binary == true
        @info "Running spatial objective function as an MINLP"
    end

    JuMP.@objective(
        model,
        Min,
        sum(
            (F[n, o, s, wildtype, t] - F[n, o, s, wildtype, 1] * (percent_to_suppress)^t) .^
            2 for n in N[node_to_suppress], o in O, s in SF, t in T
        ) +
        costquant * sum(control_M) +
        costloc * sum(release_location[n, t] for n in N, t in T)
    )

    # Run
    JuMP.optimize!(model)

    # Output
    if termination_status(model) != OPTIMAL
        @info("Termination status: $(termination_status(model))")
        return model
    else
        println("Objective value:", objective_value(model))
    end

    return model
end

function solve_decision_model_wolb(
    model::JuMP.Model,
    objective_function::Type{<:SpatialReduction};
    wildtype=nothing,
    homozygous_modified=nothing,
    node_to_suppress=nothing,
    percent_to_suppress=nothing,
    do_binary=false,
    policy=nothing,
)

    # Variables
    control_M = model[:control_M]
    F = model[:F]
    control_F = model[:control_F]

    # Sets
    sets = model.obj_dict[:Sets]
    N = sets[:N]
    O = sets[:O]
    SF = sets[:SF]
    SM = sets[:SM]
    G = sets[:G]
    T = sets[:T]

    # Expression 
    costquant = 1e-8
    costloc = 1e-8
    release_location = model[:release_location]

    if do_binary
        @info "Running spatial objective function as an MINLP"
    end

    JuMP.@objective(
        model,
        Min,
        sum(
            (F[n, o, s, wildtype, t] - F[n, o, s, wildtype, 1] * (percent_to_suppress)^t) .^
            2 for n in N[node_to_suppress], o in O, s in SF, t in T
        ) +
        costquant * sum(control_M) +
        costquant * sum(control_F)
    )
    +costloc * sum(release_location[n, t] for n in N, t in T)

    # Run
    JuMP.optimize!(model)

    # Output
    if termination_status(model) != OPTIMAL
        @info("Termination status: $(termination_status(model))")
        return model
    else
        println("Objective value:", objective_value(model))
    end

    return model
end
