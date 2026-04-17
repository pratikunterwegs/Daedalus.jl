
module Events

using ..Constants
using ..DaedalusStructs
using ..Helpers

using DiffEqCallbacks: SavedValues, SavingCallback, PresetTimeCallback
using LinearAlgebra: eigen, norm
using OrdinaryDiffEq

export make_events, make_param_changer, make_param_reset,
       make_save_events, make_rt_logger, get_coef

"""
    make_param_changer(eff::Effect)::Function

Create a callback function that modifies a parameter when activated.

Handles both `ReactiveTrigger` (state-dependent) and `TimeTrigger` (time-based)
activation.

# Arguments
- `eff::Effect`: Effect with trigger conditions

# Returns
A `Function` that modifies the target parameter in-place
"""
function make_param_changer(eff::Effect)::Function
    effect! = function () end
    if isa(eff.trigger_on, TimeTrigger)
        effect! = function (integrator)
            original = getproperty(integrator.p, eff.target)
            new_val = eff.func(original)
            setproperty!(integrator.p, eff.target, new_val)
        end
        return effect!
    else
        effect! = function (integrator)
            if length(eff.saved_values.saveval) > 0
                value_on = eff.trigger_on.value
                u = eff.saved_values.saveval[end][1] # index 1 for comp_on
                if u > value_on && !eff.ison
                    eff.ison = true
                    original = getproperty(integrator.p, eff.target)
                    new_val = eff.func(original)
                    setproperty!(integrator.p, eff.target, new_val)
                end
            else
                nothing
            end
        end
        return effect!
    end
end

"""
    make_param_reset(eff::Effect)::Function

Create a callback function that resets a parameter when deactivated.

Handles both `ReactiveTrigger` (state-dependent) and `TimeTrigger` (time-based)
deactivation.

# Arguments
- `eff::Effect`: Effect with trigger conditions

# Returns
A `Function` that resets the target parameter in-place
"""
function make_param_reset(eff::Effect)::Function
    effect! = function () end
    if isa(eff.trigger_on, TimeTrigger)
        effect! = function (integrator)
            current = getproperty(integrator.p, eff.target)
            reset_val = eff.reset_func(current)
            setproperty!(integrator.p, eff.target, reset_val)
        end

        return effect!
    else
        effect! = function (integrator)
            if length(eff.saved_values.saveval) > 0
                value_off = eff.trigger_off.value
                u = eff.saved_values.saveval[end][2] # index 2 for trigger_off

                if u < value_off && eff.ison
                    eff.ison = false
                    current = getproperty(integrator.p, eff.target)
                    reset_val = eff.reset_func(current)
                    setproperty!(integrator.p, eff.target, reset_val)
                end
            else
                nothing
            end
        end

        return effect!
    end
end

"""
    make_save_events(eff::ParamEffect, savepoints)::Union{SavingCallback, Nothing}

Create a SavingCallback for an effect that records trigger compartment values.

Skips effects with all `TimeTrigger` conditions (no state tracking needed).
For effects with at least one `ReactiveTrigger`, captures trigger compartment
values at each savepoint.

# Arguments
- `eff::ParamEffect`: Effect with trigger conditions
- `savepoints`: Time points at which to save state values

# Returns
A `SavingCallback` for reactive effects, or `nothing` if both triggers are
    time-based
"""
function make_save_events(eff::ParamEffect, savepoints)
    # Skip if both triggers are time-based (this is a timed effect)
    if isa(eff.trigger_on, TimeTrigger) && isa(eff.trigger_off, TimeTrigger)
        return nothing
    elseif isa(eff.trigger_on, TimeTrigger)
        idx_off = Constants.get_indices(eff.trigger_off.name)
        savingcb = SavingCallback(
            (u, t, integrator) -> begin
                sum_off = sum(@view u[idx_off])
                return (t, sum_off)
            end,
            eff.saved_values, saveat = savepoints
        )
        return savingcb
    else
        # both are reactive
        idx_on = Constants.get_indices(eff.trigger_on.name)
        idx_off = Constants.get_indices(eff.trigger_off.name)
        savingcb = SavingCallback(
            (u, t, integrator) -> begin
                sum_on = sum(@view u[idx_on])
                sum_off = sum(@view u[idx_off])
                return (sum_on, sum_off)
            end,
            eff.saved_values, saveat = savepoints
        )
        return savingcb
    end
end

"""
    make_save_events(npi::Npi, savepoints)

Create SavingCallbacks for all effects that require state tracking.

Returns a vector of `SavingCallback`s, one per effect with at least one
`ReactiveTrigger`. Effects with only `TimeTrigger` conditions are skipped.

# Arguments
- `npi::Npi`: Npi struct with effects
- `savepoints`: Time points at which to save state values

# Returns
A `Vector{SavingCallback}` — one per effect that has reactive triggers
"""
function make_save_events(npi::Npi, savepoints)
    cbset = []
    for eff in npi.effects
        cb = make_save_events(eff, savepoints)
        if !isnothing(cb)
            push!(cbset, cb)
        end
    end
    return cbset
end

"""
    make_save_events(events::Vector{Event}, savepoints)

Create SavingCallbacks for all events that require state tracking.

Only Npi events produce save callbacks; Vaccination events do not.

# Arguments
- `events::Vector{Event}`: Vector of events
- `savepoints`: Time points at which to save state values

# Returns
A `Vector{SavingCallback}` for all reactive effects in Npi events
"""
function make_save_events(events::Vector{DaedalusStructs.Event}, savepoints)
    cbset = []
    for ev in events
        if isa(ev, DaedalusStructs.Npi)
            append!(cbset, make_save_events(ev, savepoints))
        end
    end
    return cbset
end

"""
    make_events(eff::ParamEffect, savepoints)::CallbackSet

Create a CallbackSet for an effect's activation and deactivation callbacks.

# Arguments
- `eff::ParamEffect`: Effect with triggers
- `savepoints`: Time points for state checking (used if trigger is 
    `ReactiveTrigger`)

# Returns
A `CallbackSet` with on/off callbacks
"""
function make_events(eff::ParamEffect, savepoints)::CallbackSet
    affect_on! = make_param_changer(eff)
    affect_off! = make_param_reset(eff)

    tpoints_on = isa(eff.trigger_on, ReactiveTrigger) ? savepoints : eff.trigger_on.value
    tpoints_off = isa(eff.trigger_off, ReactiveTrigger) ? savepoints : eff.trigger_off.value

    cb_on = PresetTimeCallback(tpoints_on, affect_on!)
    cb_off = PresetTimeCallback(tpoints_off, affect_off!)

    return CallbackSet(cb_on, cb_off)
end

"""
    make_events(npi::Npi, savepoints)::CallbackSet

Create a CallbackSet of callbacks for all effects in an NPI.

# Arguments
- `npi::Npi`: Npi struct with effects
- `savepoints`: Time points for state checking (for reactive triggers)

# Returns
A `CallbackSet` with all activation/deactivation callbacks
"""
function make_events(npi::Npi, savepoints)::CallbackSet
    cbset = CallbackSet()

    for eff in npi.effects
        cbset_eff = make_events(eff, savepoints)

        # cbset grown flexibly
        cbset = CallbackSet(cbset, cbset_eff)
    end

    return cbset
end

"""
    make_events(vax::Vaccination, savepoints)::CallbackSet

Create a CallbackSet for vaccination campaign callbacks.

Generates two callbacks:
1. Start callback: fires at vax.start_time, sets nu and psi parameters
2. Saturation callback: fires at savepoints, stops vaccination when coverage limit reached

# Arguments
- `vax::Vaccination`: Vaccination event with parameters
- `savepoints`: Time points for checking coverage saturation

# Returns
A `CallbackSet` with start and saturation callbacks
"""
function make_events(vax::DaedalusStructs.Vaccination, savepoints)::CallbackSet
    # Start callback: fires once at vaccination start time
    start_affect! = function (integrator)
        if !vax.ison
            vax.ison = true
            integrator.p.nu = vax.rate
            integrator.p.psi = 1.0 / vax.waning_period
        end
    end
    cb_start = PresetTimeCallback([vax.start_time], start_affect!)

    idx_vaxxed = get_indices("vax")

    # Saturation callback: fires at savepoints to check if coverage limit reached
    saturation_affect! = function (integrator)
        if vax.ison
            # sum all vaccinated including dead
            vax_pop = sum(@view integrator.u[idx_vaxxed])
            total_pop = sum(integrator.p.demography)
            coverage = vax_pop / total_pop

            if coverage >= vax.uptake_limit
                integrator.p.nu = 0.0
            end
        end
    end
    cb_sat = PresetTimeCallback(savepoints, saturation_affect!)

    return CallbackSet(cb_start, cb_sat)
end

"""
    make_events(events::Vector{Event}, savepoints)::CallbackSet

Create a unified CallbackSet for all events (NPI and Vaccination).

# Arguments
- `events::Vector{Event}`: Vector containing Npi and/or Vaccination events
- `savepoints`: Time points for state checking

# Returns
A `CallbackSet` combining all event callbacks
"""
function make_events(events::Vector{DaedalusStructs.Event}, savepoints)::CallbackSet
    cbset = CallbackSet()
    for ev in events
        if isa(ev, DaedalusStructs.Npi)
            cbset = CallbackSet(cbset, make_events(ev::DaedalusStructs.Npi, savepoints))
        elseif isa(ev, DaedalusStructs.Vaccination)
            cbset = CallbackSet(cbset, make_events(ev::DaedalusStructs.Vaccination, savepoints))
        end
    end
    return cbset
end

"""
    make_rt_logger(savepoints)

Create a callback that computes and logs the effective reproduction number Rt.

Uses power iteration on the susceptibility-adjusted NGM to compute Rt at each
savepoint. Maintains a warm-start eigenvector for improved convergence.

# Arguments
- `savepoints`: Time points at which to compute and log Rt

# Returns
A `PresetTimeCallback` that updates the Rt state variable at each savepoint
"""
function make_rt_logger(savepoints)
    # make a PresetTimeCallback that updates u at integerish times
    iRt = Constants.get_indices("Rt")

    # Preallocate initial vector for power iteration warm starts
    # Using warm starts can improve convergence speed between timesteps
    v_prev = randn(N_TOTAL_GROUPS)
    v_prev = v_prev / norm(v_prev)

    idx_S = Constants.get_indices("S")

    function affect!(integrator)
        S = @view integrator.u[idx_S] # hardcoded, as vax comp is empty

        p_susc = S ./ integrator.p.demography
        ngm_susc = integrator.p.ngm .* p_susc

        # Use power iteration to compute only the dominant eigenvalue
        # This is significantly faster than computing all eigenvalues
        rt = Helpers.dominant_eigenvalue(
            ngm_susc, v_init = v_prev, max_iter = 100, tol = 1e-5)

        integrator.u[iRt] = rt # only logged here, use explicit index rather than end
    end

    pstcb_rt = PresetTimeCallback(savepoints, affect!)

    return pstcb_rt
end

end
