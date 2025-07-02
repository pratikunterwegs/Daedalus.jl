
module Events

using ..Constants
using DifferentialEquations

export make_state_condition, make_time_condition, start_vax!, reduce_beta!

"""
    make_cond(threshold)::Function

Factory function for conditions. Makes a function that checks whether a root is
    found at the sum of a state index.
"""
function make_state_condition(threshold, index)::Function
    function fn_cond(u, t, integrator)
        X = @view u[:, index, :]
        state_sum = sum(X)

        state_sum - threshold
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
    reduce_beta!(integrator)

An event function that reduces beta by a fixed value.
"""
function reduce_beta!(integrator)
    integrator.p.beta *= 0.8
end

"""
    restore_beta!(integrator)

An event function that restores beta by a fixed value.
"""
function restore_beta!(integrator)
    integrator.p.beta *= 1.6
end

end
