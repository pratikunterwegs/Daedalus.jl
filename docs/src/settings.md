```@meta
CurrentModule = Daedalus
```

# Multiple contact settings

`CountryData.contact_matrix` accepts either a single `Matrix{Float64}` or a `Vector{Matrix{Float64}}`.
Passing a vector lets the model treat distinct social contact regimes.
For example, community mixing and workplace mixing — as separate matrices that are combined inside the ODE force-of-infection calculation.

## Creating a multi-setting CountryData struct

Retrieve a country struct, then replace its `contact_matrix` field with a
vector of matrices:

```@example multi_settings
using Daedalus

cd = Daedalus.DataLoader.get_country("Australia")
cm = deepcopy(cd.contact_matrix) # baseline 4×4 contact matrix
cm_work = cm .* 2.0 # a second setting with doubled contacts

cd_multi = deepcopy(cd)
cd_multi.contact_matrix = [cm, cm_work]
```

`Data.get_settings` returns the number of active contact settings:

```@example multi_settings
Daedalus.Data.get_settings(cd)
```

```@example multi_settings
Daedalus.Data.get_settings(cd_multi)
```

## Running the model with multiple settings

Pass the modified struct directly to [`daedalus`](@ref) along with infection parameters:

```@example multi_settings
result_multi = daedalus(cd_multi, "sars-cov-2 delta", time_end = 300.0)

times  = Daedalus.Outputs.get_times(result_multi)
deaths = Daedalus.Outputs.get_values(result_multi, "D", 1)
println("Total deaths at end: ", round(Int, last(deaths)))
```

## Calculating β with multiple settings

`daedalus` computes `β` from the **sum of all contact matrices** via [`Data.total_contacts`](@ref):

```julia
# NOTE: example, this code is not run
contacts_unscaled = total_contacts(prepare_contacts(cd; scaled = false))
beta = get_beta(contacts_unscaled, r0, ...)
```

The effective reproduction number is therefore calibrated against the *total* contact rate across all settings.
As a consequence, splitting contacts across multiple equal settings does not change the epidemic.
Adding a second identical matrix doubles the total contacts, which halves `β`.

```@example multi_settings
cd_double = deepcopy(cd)
cd_double.contact_matrix = [cm, cm]

infection = Daedalus.DataLoader.get_pathogen("sars-cov-2 delta");

# get beta used for each country
beta_single = Daedalus.Helpers.get_beta(
    Daedalus.Data.total_contacts(
    Daedalus.Data.prepare_contacts(cd; scaled=false)
), infection)

beta_double = Daedalus.Helpers.get_beta(
    Daedalus.Data.total_contacts(
    Daedalus.Data.prepare_contacts(cd_double; scaled=false)
), infection)
```

Print ``\beta`` values to check that `beta_double` is half of `beta_single`.

```@example multi_settings
println("β (1 setting): $beta_single")
println("β (2 settings): $beta_double")
```

## How settings appear in the ODE

The matrices are stacked into a 3D array of shape `(N_TOTAL_GROUPS, N_TOTAL_GROUPS, K)` where `K` is the number of settings.
In each ODE evaluation, a weight vector of length `K` contracts the third dimension via [`Helpers.weighted_slice_sum!`](@ref):

```@example multi_settings
c3d = Daedalus.Data.contacts3d(cd_multi)
size(c3d) # (49, 49, 2) for two settings
```

Current weights are `ones(K)` (all settings fully active).
This design is intended to support future NPI scaling, where individual settings — for example, workplace contacts during an economic closure — can be reduced independently without affecting community contacts.
```
