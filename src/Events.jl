
module Events

using ..Constants
using ..DaedalusStructs
using ..Helpers

using DiffEqCallbacks: SavedValues, SavingCallback, PresetTimeCallback
using LinearAlgebra: eigen, norm
using OrdinaryDiffEq

export make_state_condition,
       make_events, make_param_changer, make_param_reset,
       get_coef, make_save_events, make_rt_logger, make_timed_npi_callbacks

"""
    make_state_condition(threshold, idx, comparison)::Function

Factory function that creates a condition function for event callbacks.

Returns a function that checks whether a state sum crosses a threshold at the
specified indices. Used to trigger events when compartment totals reach certain values.

# Arguments
- `threshold::Float64`: The threshold value for comparison
- `idx`: Indices into the state vector to sum
- `comparison::Function`: Comparison function (e.g., `<`, `>`)

# Returns
A `Function` that takes `(u, t, integrator)` and returns 0.0 if condition is met, 1.0 otherwise
"""
function make_state_condition(threshold, idx, comparison::Function)::Function
    function fn_cond(u, t, integrator)
        state_sum = 0.0

        X = @view u[idx]
        state_sum = sum(X)

        return comparison(state_sum, threshold) ? 0.0 : 1.0
    end

    return fn_cond
end

"""
    make_param_changer(param_name, func, coef)::Function

Create a callback function that modifies a parameter during ODE integration.

Generates an effect function that applies `func` to the current parameter value
multiplied by `coef`, storing the result in `param_name_now`.

# Arguments
- `param_name::String`: Name of the parameter to modify (e.g., `"beta"`)
- `func::Function`: Function to apply (e.g., `.*` for multiplication)
- `coef`: Coefficient to apply via `func`

# Returns
A `Function` that takes an `integrator` and modifies the parameter in-place
"""
function make_param_changer(param_name, func, coef)::Function
    function effect!(integrator)
        new_value = func(getproperty(integrator.p, Symbol(param_name)), coef)
        setproperty!(integrator.p, Symbol(param_name * "_now"), new_value)
    end
end

"""
    make_param_reset(param_name)::Function

Create a callback function that resets a parameter to its original value.

# Arguments
- `param_name::String`: Name of the parameter to reset (e.g., `"beta"`)

# Returns
A `Function` that takes an `integrator` and resets `param_name_now` to the original `param_name`
"""
function make_param_reset(param_name)::Function
    function effect!(integrator)
        original_value = getproperty(integrator.p, Symbol(param_name))
        setproperty!(integrator.p, Symbol(param_name * "_now"), original_value)
    end
end

"""
    make_events(x::Npi, effect_on, effect_off, savepoints)::CallbackSet

Create a CallbackSet of event callbacks for a reactive (state-dependent) NPI.

Generates paired callbacks that trigger intervention activation and deactivation
based on epidemic state (e.g., hospitalizations reaching a threshold).

# Arguments
- `x::Npi`: Reactive NPI struct
- `effect_on::Function`: Callback to execute when intervention activates
- `effect_off::Function`: Callback to execute when intervention deactivates
- `savepoints`: Time points at which to check state for triggers

# Returns
A `CallbackSet` with on/off event callbacks
"""
function make_events(x::Npi, effect_on::Function, effect_off::Function, savepoints)::CallbackSet
    idx_on = DaedalusStructs.get_indices(x.params.comp_on)
    idx_off = DaedalusStructs.get_indices(x.params.comp_off)

    value_on = get_comp_threshold(x.params.comp_on)
    value_off = get_comp_threshold(x.params.comp_off)

    # pass appropriate fn `>` or `<` for comparison
    # checking against U can be replaced by checking x.saved_values
    function affect_on!(integrator)
        if length(x.saved_values.saveval) < 2 # need two values for slope
            return nothing
        end

        u = x.saved_values.saveval[end][2] # 1st value is flag
        _u = x.saved_values.saveval[end - 1][2]
        growing = u > _u

        if u > value_on && growing && !x.ison # only supporting gt for now
            x.ison = true
            effect_on(integrator)
        else
            nothing
        end
    end

    cb_on = PresetTimeCallback(savepoints, affect_on!)

    function affect_off!(integrator)
        if length(x.saved_values.saveval) == 0
            return nothing
        end

        u = x.saved_values.saveval[end][3]
        if u < value_off && x.ison # only supporting lt for now
            x.ison = false
            effect_off(integrator)
        else
            nothing
        end
    end

    cb_off = PresetTimeCallback(savepoints, affect_off!)

    return CallbackSet(cb_on, cb_off)
end

"""
    get_coef(x::Npi)::Float64

Extract the transmission reduction coefficient from an Npi struct.

# Arguments
- `x::Npi`: Reactive NPI struct

# Returns
The transmission coefficient (typically in [0, 1])
"""
function get_coef(x::Npi)::Float64
    return x.coefs.coef
end

"""
    make_save_events(x::Npi, savepoints)::SavingCallback

Create a SavingCallback that records state values at specified timepoints.

Captures NPI-on state and values of the NPI trigger and deactivation compartments
at each savepoint, storing results in the Npi's `saved_values` field.

# Arguments
- `x::Npi`: Reactive NPI struct with `saved_values` field
- `savepoints`: Time points at which to save state values

# Returns
A `SavingCallback` that writes to `x.saved_values`
"""
function make_save_events(x::Npi, savepoints)::SavingCallback
    # the saving callbacks save to x.saved_values
    savingcb = SavingCallback(
        (u, t, integrator) -> begin
            # handle comp on
            idx_on = DaedalusStructs.get_indices(x.params.comp_on)
            sum_on = sum(u[idx_on])

            # handle comp off, usually Rt
            idx_off = DaedalusStructs.get_indices(x.params.comp_off)
            sum_off = sum(u[idx_off])

            return (x.ison, sum_on, sum_off)
        end,
        x.saved_values, saveat = savepoints
    )

    return savingcb
end

"""
    make_rt_logger(savepoints)::PresetTimeCallback

Create a callback that computes and logs the effective reproduction number Rt.

Uses power iteration on the susceptibility-adjusted NGM to compute Rt at each
savepoint. Maintains a warm-start eigenvector for improved convergence.

# Arguments
- `savepoints`: Time points at which to compute and log Rt

# Returns
A `PresetTimeCallback` that updates the Rt state variable at each savepoint
"""
function make_rt_logger(savepoints)::PresetTimeCallback
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
function make_timed_npi_callbacks(npi::DaedalusStructs.TimedNpi)
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
