
module Events

using ..Constants
using ..DaedalusStructs
using ..Helpers

using DiffEqCallbacks: SavedValues, SavingCallback, PresetTimeCallback
using LinearAlgebra: eigen, norm
using OrdinaryDiffEq

export make_events, make_param_changer, make_param_reset,
       make_save_events, make_rt_logger, make_timed_npi_callbacks, get_coef

"""
    make_param_changer(npi::Npi)::Function

Create a callback function that modifies all parameters specified in an Npi.

Reads the parameter effects from `npi.effects` and generates a single function that
applies all modifications at once.

# Arguments
- `npi::Npi`: An Npi struct containing a list of ParamEffect specifications

# Returns
A `Function` that takes an `integrator` and modifies all affected parameters in-place
"""
function make_param_changer(eff::ParamEffect)::Function
    function effect!(integrator)
        if length(eff.saved_values.saveval) > 0
            value_on = eff.comp_on.value
            u = eff.saved_values.saveval[end][1] # index 1 for comp_on

            if u > value_on && !eff.ison
                eff.ison = true
                original = getproperty(integrator.p, eff.target)
                new_val = eff.func(original)
                # TODO: remove
                # time = integrator.t
                # println("t = $time Current value $(eff.target) = $original")
                # println("Setting value $(eff.target) = $new_val")
                setproperty!(integrator.p, eff.target, new_val)
            end
        else
            nothing
        end
    end
    return effect!
end

"""
    make_param_reset(npi::Npi)::Function

Create a callback function that resets all parameters specified in an Npi.

Reads the parameter names from `npi.effects` and generates a single function that
resets all of them to their original values.

# Arguments
- `npi::Npi`: An Npi struct containing a list of ParamEffect specifications

# Returns
A `Function` that takes an `integrator` and resets all affected parameters in-place
"""
function make_param_reset(eff::ParamEffect)::Function
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
    make_save_events(npi::Npi, savepoints)

Create SavingCallbacks that record compartment values at specified timepoints.

Returns a vector of `SavingCallback`s, one per effect. Each callback captures the
values of that effect's on/off trigger compartments at each savepoint.

# Arguments
- `npi::Npi`: Reactive NPI struct with a vector of ParamEffect specifications
- `savepoints`: Time points at which to save state values

# Returns
A `Vector{SavingCallback}` — one callback per effect, each writing to `eff.saved_values`
"""
function make_save_events(npi::Npi, savepoints)
    callbacks = []
    for eff in npi.effects
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

Create a CallbackSet of event callbacks for a reactive (state-dependent) NPI.

Generates paired callbacks (on/off) for each effect in the NPI. Each effect's callbacks
trigger independently based on that effect's own trigger conditions.

# Arguments
- `npi::Npi`: Reactive NPI struct with a vector of ParamEffect specifications
- `savepoints`: Time points at which to check state for triggers

# Returns
A `CallbackSet` with all on/off callbacks for all effects
"""
function make_events(npi::Npi, savepoints)::CallbackSet
    cbset = CallbackSet()

    for eff in npi.effects
        # make PSTCb for effect.on for this effect
        affect_on!::Function = make_param_changer(eff)
        cb_on = PresetTimeCallback(savepoints, affect_on!)
        cbset = CallbackSet(cbset, cb_on)

        # make PSTCb for effect.off for this effect
        affect_off!::Function = make_param_reset(eff)
        cb_off = PresetTimeCallback(savepoints, affect_off!)
        cbset = CallbackSet(cbset, cb_off)
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

"""
    make_timed_npi_callbacks(npi::TimedNpi)

Create a set of callbacks for time-limited NPIs that activate and deactivate
at specified times. Each intervention phase applies a transmission reduction
coefficient to the beta parameter.

# Arguments
- `npi::TimedNpi`: A TimedNpi struct containing intervention timing and coefficients

# Returns
- `CallbackSet`: A set of PresetTimeCallback objects that modify beta at the
  specified intervention times

# Details
The function creates paired callbacks for each intervention phase:
- At `start_times[i]`, beta is multiplied by `coefs[i]`
- At `end_times[i]`, beta is reset to its original value

Between phases, the original beta value is restored. Multiple phases can be
specified to create complex intervention scenarios.

# Example
```julia
timed_npi = TimedNpi([10.0, 30.0], [20.0, 40.0], [0.5, 0.3], "two_phase")
callbacks = make_timed_npi_callbacks(timed_npi)
```
"""
function make_timed_npi_callbacks(npi::DaedalusStructs.TimedNpi)::CallbackSet
    callbacks = []

    for (t_on, t_off, coef) in zip(npi.start_times, npi.end_times, npi.coefs)
        # Create effect function for NPI activation
        # Multiply beta by the coefficient
        affect_on! = function (integrator)
            integrator.p.beta_now = integrator.p.beta * coef
        end

        # Create callback for activation
        cb_on = PresetTimeCallback(t_on, affect_on!)
        push!(callbacks, cb_on)

        # Create effect function for NPI deactivation
        # Reset beta to original value
        affect_off! = function (integrator)
            integrator.p.beta_now = integrator.p.beta
        end

        # Create callback for deactivation
        cb_off = PresetTimeCallback(t_off, affect_off!)
        push!(callbacks, cb_off)
    end

    return CallbackSet(callbacks...)
end

end
