
module Events

using ..Constants

export make_cond, cond_vax, start_vax!, reduce_beta!

"""
    make_cond(threshold)::Function

Factory function for conditions.
"""
function make_cond(threshold)::Function
    function fn_cond(u, t, integrator)
        H = @view u[:, iH, :]
        total_hosp = sum(H)

        total_hosp - threshold
    end

    return fn_cond
end

"""
    condition(u, t, integrator)

A condition function that triggers an event at a specific time.
"""
function cond_vax(u, t, integrator) # Event when condition(u,t,integrator) == 0
    # trigger vaccination at t_vax, access by position
    t == integrator.p[15]
end

"""
    affect!(integrator)

An event function.
"""
function start_vax!(integrator)
    integrator.p[16] = true # initial value is 0.0
end

function reduce_beta!(integrator)
    integrator.p[3] *= 0.2
end

end