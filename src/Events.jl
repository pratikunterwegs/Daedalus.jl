
module Events

using ..Constants
using ..DaedalusStructs
using ..Helpers

using DiffEqCallbacks: SavedValues, SavingCallback, PresetTimeCallback
using LinearAlgebra: eigen
using OrdinaryDiffEq

export make_state_condition,
    make_events, make_param_changer, make_param_reset,
    get_coef, make_save_events, make_rt_logger

"""
    make_cond(threshold)::Function

Factory function for conditions. Makes a function that checks whether a root is
    found at the sum of a state index. Passing `crossing = "up"` checks for an
    increasing root (-1 to +1), while `crossing = "down"` checks for a
    decreasing root (+1 to -1).
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
    make_param_changer(param_name, func, coef)

A function factory to generate effect functions.
"""
function make_param_changer(param_name, func, coef)
    function effect!(integrator)
        println("time = ", integrator.t, "; changing param ", param_name, "!")
        new_value = func(getproperty(integrator.p, Symbol(param_name)), coef)
        setproperty!(integrator.p, Symbol(param_name * "_now"), new_value)
    end
end

function make_param_reset(param_name)
    function effect!(integrator)
        println("time = ", integrator.t, "; resetting param ", "!")
        original_value = getproperty(integrator.p, Symbol(param_name))
        setproperty!(integrator.p, Symbol(param_name * "_now"), original_value)
    end
end

"""
    make_events(x::Npi)

Make a CallbackSet from a Npi struct.
"""
make_events(x::Npi, effect_on::Function, effect_off::Function, savepoints) = begin
    idx_on = DaedalusStructs.get_indices(x.params.comp_on)
    idx_off = DaedalusStructs.get_indices(x.params.comp_off)

    value_on = get_comp_threshold(x.params.comp_on)
    value_off = get_comp_threshold(x.params.comp_off)

    # pass appropriate fn `>` or `<` for comparison
    # checking against U can be replaced by checking x.saved_values
    function affect_on!(integrator)
        _u = similar(integrator.u)
        get_du!(_u, integrator)
        growing = any(_u[idx_on] .> 0.0)

        U = @view integrator.u[idx_on]
        if sum(U) > value_on && growing && !x.ison # only supporting gt for now
            effect_on(integrator)
            x.ison = true
        else
            nothing
        end
    end

    cb_on = PresetTimeCallback(savepoints, affect_on!)

    function affect_off!(integrator)
        U = @view integrator.u[idx_off]
        if sum(U) < value_off && x.ison # only supporting lt for now
            effect_off(integrator)
            x.ison = false
        else
            nothing
        end
    end

    cb_off = PresetTimeCallback(savepoints, affect_off!)

    return CallbackSet(cb_on, cb_off)
end

get_coef(x::Npi) = begin
    return x.coefs.coef
end

"""
    make_save_events(x::Npi, savepoints)

Make a CallbackSet of SavingCallbacks from an Npi struct. `savepoints` is
    expected to be a float range.
"""
make_save_events(x::Npi, savepoints) = begin
    # the saving callbacks save to x.saved_values
    savingcb = SavingCallback(
        (u, t, integrator) -> begin
            # handle comp on
            idx_on = DaedalusStructs.get_indices(x.params.comp_on)
            sum_on = sum(u[idx_on])

            # handle comp off, usually Rt
            idx_off = DaedalusStructs.get_indices(x.params.comp_off)
            sum_off = sum(u[idx_off])

            return (sum_on, sum_off)
        end,
        x.saved_values, saveat=savepoints
    )

    return savingcb
end

"""
    make_rt_logger(savepoints)

Make a PresetTimeCallback to update Rt.
"""
make_rt_logger(savepoints) = begin
    # make a PresetTimeCallback that updates u at integerish times
    iRt = Constants.get_indices("Rt")
    function affect!(integrator)

        U = @view integrator.u[1:integrator.p.size]
        p_susc = sum_by_age(U, iS) ./ integrator.p.demography
        ngm_susc = integrator.p.ngm .* p_susc
        rt = maximum(eigen(ngm_susc).values)

        integrator.u[iRt] = rt
    end

    pstcb_rt = PresetTimeCallback(savepoints, affect!)

    return pstcb_rt
end

end
