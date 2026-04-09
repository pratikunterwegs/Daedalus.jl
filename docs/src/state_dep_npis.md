```@meta
CurrentModule = Daedalus
```

# State-dependent interventions

This module shows how to use the `Npi` types to specify rule-based interventions during epidemics.
The idea is to be able to simulate policies that come into effect based on the state of the epidemic.

This is in contrast with time-dependent NPIs which require policy timespans to be known in advance.
State-dependent NPIs may be useful when there is clarity on the system state that would trigger a response (e.g. '80% of emergency hospital capacity is exhausted').

## Reactive NPI struct

`Npi` is a flexible framework for reactive (state-dependent) interventions. Each `ParamEffect` specifies:
- Which parameter to modify (e.g., `:beta`, `:omega`)
- How to modify it (transformation function)
- When to activate and deactivate (independent trigger conditions per effect)

## Basic example: Single parameter with prevalence trigger

This example shows how to implement an NPI that reduces transmission by 20% once hospital occupancy crosses 20000 beds.
The NPI is removed once hospital occupancy falls below 10000 beds.

This is the simplest kind of reactive event, as it is conditioned on the state of the system `x`.
For a compartmental epidemic model this equates to an event that is reactive to prevalence.

```@example reactive_npi
using Daedalus
using Plots

# Create an effect: reduce beta by 20% when hospitalizations exceed 20000
trigger_on = Daedalus.DaedalusStructs.ReactiveTrigger(20000.0, "H")
trigger_off = Daedalus.DaedalusStructs.ReactiveTrigger(10000.0, "H")
effect = Daedalus.DaedalusStructs.ParamEffect(
    :beta,
    x -> x .* 0.8,         # reduce by 20%
    x -> x ./ 0.8,         # reset: divide by 0.8
    trigger_on,             # Activate when H > 20000
    trigger_off             # Deactivate when H < 10000
)

# Create NPI from the effect
npi = Daedalus.DaedalusStructs.Npi([effect])

data = daedalus("Australia", "sars-cov-2 delta", npi=npi, time_end=600.0);
data_default = daedalus("Australia", "sars-cov-2 delta", time_end=600.0);

# plot hospitalizations over time
hosp_default = Daedalus.Outputs.get_values(data_default, "H", 1)
hosp_npi = Daedalus.Outputs.get_values(data, "H", 1)

rt = Daedalus.Outputs.get_values(data, "Rt", 1)
rt_default = Daedalus.Outputs.get_values(data_default, "Rt", 1)

times = Daedalus.Outputs.get_times(data)

plot(times, hosp_npi, label="Hosp w/ NPI")
plot!(times, hosp_default, label="Hosp w/o NPI")
xlabel!("Time (days)")
ylabel!("Hosp occupancy")
```

Plot $R_t$ to see the effect on transmission reduction.

```@example reactive_npi
plot(times, rt, label="Rt w/ NPI")
plot!(times, rt_default, label="Rt w/o NPI")
xlabel!("Time (days)")
ylabel!("Rt")
```

## Advanced example: Multiple parameters with independent triggers

Each effect can respond to different state variables independently:

```@example multi_param_npi
using Daedalus
using Plots

# Create multiple effects with different triggers
trigger_beta_on = Daedalus.DaedalusStructs.ReactiveTrigger(5000.0, "H")
trigger_beta_off = Daedalus.DaedalusStructs.ReactiveTrigger(1.0, "Rt")
effect_beta = Daedalus.DaedalusStructs.ParamEffect(
    :beta,
    x -> x .* 0.8,         # reduce by 20%
    x -> x ./ 0.8,         # reset: divide by 0.8
    trigger_beta_on,        # Activate on hospitalizations
    trigger_beta_off        # Deactivate when Rt < 1.0
)

hosp_capacity = 20e3
trigger_omega_on = Daedalus.DaedalusStructs.ReactiveTrigger(hosp_capacity, "H")
trigger_omega_off = Daedalus.DaedalusStructs.ReactiveTrigger(hosp_capacity / 2.0, "H")
effect_omega = Daedalus.DaedalusStructs.ParamEffect(
    :omega,
    x -> x .* 1.2,         # increase by 20%
    x -> x ./ 1.2,         # reset: divide by 1.2
    trigger_omega_on,       # Activate at high occupancy
    trigger_omega_off       # Deactivate at medium occupancy
)

# Create NPI with both effects
npi = Daedalus.DaedalusStructs.Npi([effect_beta, effect_omega])

result = daedalus("Australia", "sars-cov-2 delta", npi=npi, time_end=600.0);
```
