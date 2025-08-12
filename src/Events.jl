
module Events

using ..Constants
using ..DaedalusStructs
using DifferentialEquations

export make_state_condition, make_time_condition, start_vax!,
    make_events, get_times, make_param_changer, make_param_reset,
    get_coef, make_value_saver

"""
    make_cond(threshold)::Function

Factory function for conditions. Makes a function that checks whether a root is
    found at the sum of a state index. Passing `crossing = "up"` checks for an
    increasing root (-1 to +1), while `crossing = "down"` checks for a
    decreasing root (+1 to -1).
"""
function make_state_condition(threshold, idx, crossing)::Function
    function fn_cond(u, t, integrator)
        X = @view u[idx]
        state_sum = sum(X)

        if crossing == "up"
            return state_sum - threshold
        elseif crossing == "down"
            return threshold - state_sum
        else
            error("State condition maker: `crossing` must be 'up' or 'down'.")
        end

        return state_sum - threshold
    end

    return fn_cond
end

"""
    condition(u, t, integrator)

A condition function that triggers an event at a specific time.
"""
function make_time_condition(times)::Function
    function fn_cond(u, t, integrator)
        t in times
    end

    return fn_cond
end

"""
    affect!(integrator)

An event function that begins vaccination by setting a vaccination flag to true.
"""
function start_vax!(integrator)
    integrator.p.switch = true # initial value is 0.0
end

"""
    make_param_changer(param_name, func, coef)

A function factory to generate effect functions.
"""
function make_param_changer(param_name, func, coef)
    function effect!(integrator)
        new_value = func(getproperty(integrator.p, Symbol(param_name)), coef)
        setproperty!(integrator.p, Symbol(param_name * "_now"), new_value)
    end
end

function make_param_reset(param_name)
    function effect!(integrator)
        original_value = getproperty(integrator.p, Symbol(param_name))
        setproperty!(integrator.p, Symbol(param_name * "_now"), original_value)
    end
end

"""
    make_events(x::TimedNpi)

Make a CallbackSet from a TimedNpi struct.
"""
make_events(x::TimedNpi, effect_on::Function, effect_off::Function) = begin
    npi_times_on = x.resparams.time_on
    npi_times_off = x.resparams.time_off

    cb_npi_on = DiscreteCallback(make_time_condition(npi_times_on), effect_on)
    cb_npi_off = DiscreteCallback(make_time_condition(npi_times_off), effect_off)

    return CallbackSet(cb_npi_on, cb_npi_off)
end

get_times(x::TimedResponseData) = begin
    return [x.time_on; x.time_off]
end

get_times(x::TimedNpi) = begin
    return get_times(x.resparams)
end

"""
    make_events(x::ReactiveNpi)

Make a CallbackSet from a ReactiveNpi struct.
"""
make_events(x::ReactiveNpi, effect_on::Function, effect_off::Function) = begin
    npi_idx_on = x.resparams.id_state_on
    npi_idx_off = x.resparams.id_state_off
    npi_value_on = x.resparams.value_state_on
    npi_value_off = x.resparams.value_state_off

    cb_npi_on = ContinuousCallback(
        make_state_condition(npi_value_on, npi_idx_on, "up"), effect_on
    )
    cb_npi_off = ContinuousCallback(
        make_state_condition(npi_value_off, npi_idx_off, "down"), effect_off
    )

    return CallbackSet(cb_npi_on, cb_npi_off)
end

get_coef(x::Npi) = x.params.coef

end
