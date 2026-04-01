
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
    make_param_changer(eff::ReactiveEffect)::Function

Create a callback function that modifies a parameter based on reactive state conditions.

# Arguments
- `eff::ReactiveEffect`: A ReactiveEffect with state-dependent trigger conditions

# Returns
A `Function` that takes an `integrator` and modifies the target parameter in-place
"""
function make_param_changer(eff::ReactiveEffect)::Function
    function effect!(integrator)
        if length(eff.saved_values.saveval) > 0
            value_on = eff.comp_on.value
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

"""
    make_param_changer(eff::TimedEffect)::Function

Create a callback function that modifies a parameter for a timed effect.

# Arguments
- `eff::TimedEffect`: A TimedEffect with fixed start/end times

# Returns
A `Function` that takes an `integrator` and modifies the target parameter in-place
"""
function make_param_changer(eff::TimedEffect)::Function
    function effect!(integrator)
        original = getproperty(integrator.p, eff.target)
        new_val = eff.func(original)
        setproperty!(integrator.p, eff.target, new_val)
    end
    return effect!
end

"""
    make_param_reset(eff::ReactiveEffect)::Function

Create a callback function that resets a parameter based on reactive state conditions.

# Arguments
- `eff::ReactiveEffect`: A ReactiveEffect with state-dependent trigger conditions

# Returns
A `Function` that takes an `integrator` and resets the target parameter in-place
"""
function make_param_reset(eff::ReactiveEffect)::Function
    function effect!(integrator)
        if length(eff.saved_values.saveval) > 0
            value_off = eff.comp_off.value
            u = eff.saved_values.saveval[end][2] # index 2 for comp_off

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

"""
    make_param_reset(eff::TimedEffect)::Function

Create a callback function that resets a parameter for a timed effect.

# Arguments
- `eff::TimedEffect`: A TimedEffect with fixed start/end times

# Returns
A `Function` that takes an `integrator` and resets the target parameter in-place
"""
function make_param_reset(eff::TimedEffect)::Function
    function effect!(integrator)
        current = getproperty(integrator.p, eff.target)
        reset_val = eff.reset_func(current)
        setproperty!(integrator.p, eff.target, reset_val)
    end
    return effect!
end

"""
    make_save_events(npi::Npi, savepoints)

Create SavingCallbacks that record compartment values at specified timepoints.

Returns a vector of `SavingCallback`s, one per reactive effect. Each callback captures the
values of that effect's on/off trigger compartments at each savepoint. TimedEffect
entries are skipped (they do not require state-based saving).

# Arguments
- `npi::Npi`: Npi struct with a vector of effect specifications (both ReactiveEffect and TimedEffect)
- `savepoints`: Time points at which to save state values

# Returns
A `Vector{SavingCallback}` — one callback per ReactiveEffect, each writing to `eff.saved_values`
"""
function make_save_events(npi::Npi, savepoints)
    callbacks = []
    for eff in npi.effects
        isa(eff, ReactiveEffect) || continue
        idx_on = Constants.get_indices(eff.comp_on.name)
        idx_off = Constants.get_indices(eff.comp_off.name)

        savingcb = SavingCallback(
            (u, t, integrator) -> begin
                sum_on = sum(@view u[idx_on])
                sum_off = sum(@view u[idx_off])
                return (sum_on, sum_off)
            end,
            eff.saved_values, saveat = savepoints
        )
        push!(callbacks, savingcb)
    end
    return callbacks
end

"""
    make_events(npi::Npi, savepoints)::CallbackSet

Create a CallbackSet of event callbacks for an NPI containing both reactive and timed effects.

For ReactiveEffect: generates paired callbacks (on/off) at savepoints, triggered by state.
For TimedEffect: generates paired callbacks at the effect's start_time and end_time.

# Arguments
- `npi::Npi`: Npi struct with a vector of effect specifications (ReactiveEffect and/or TimedEffect)
- `savepoints`: Time points at which to check state for reactive effect triggers

# Returns
A `CallbackSet` with all on/off callbacks for all effects
"""
function make_events(npi::Npi, savepoints)::CallbackSet
    cbset = CallbackSet()

    for eff in npi.effects
        affect_on! = make_param_changer(eff)
        affect_off! = make_param_reset(eff)

        if isa(eff, ReactiveEffect)
            # Reactive effects: trigger at savepoints based on state
            cb_on = PresetTimeCallback(savepoints, affect_on!)
            cb_off = PresetTimeCallback(savepoints, affect_off!)
        else  # TimedEffect
            # Timed effects: trigger at fixed times
            cb_on = PresetTimeCallback(eff.start_time, affect_on!)
            cb_off = PresetTimeCallback(eff.end_time, affect_off!)
        end

        cbset = CallbackSet(cbset, cb_on, cb_off)
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
