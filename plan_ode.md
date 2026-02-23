# Plan: Implement Workplace-Specific Contacts in Daedalus.jl

## Background

The R/C++ `{daedalus}` package models three distinct routes of transmission in its ODE system:

1. **Community contacts** — a 49×49 contact matrix over all strata (4 age groups + 45 economic sectors)
2. **Within-sector workplace contacts** — a 45-element vector `cm_work` of contacts between workers in the same sector
3. **Consumer-to-worker contacts** — a 45×4 matrix `cm_cons_work` of contacts from community members to workers in each sector

These are combined in the force of infection (FOI) as:

```
FOI_comm  = β × cm × I                         (all groups)
FOI_work  = β × cm_work × I_workers × ω²      (workers only, ω = openness)
FOI_cw    = β × cm_cons_work × I_age × ω      (workers only)
```

New infections for community members use only `FOI_comm`. New infections for workers use `FOI_comm + FOI_work + FOI_cw`.

The Julia package currently **blends all three routes** into a single 49×49 contact matrix in `prepare_contacts()` (see `Data.jl:53-70`). This produces correct steady-state behaviour but prevents:
- Separate scaling of workplace transmission by NPI openness coefficients
- Accurate FOI decomposition matching the reference implementation
- Future sector-specific closure modelling

---

## Current State

### Data.jl
- `prepare_contacts()` builds a single 49×49 community matrix that adds consumer-worker and within-sector diagonal contacts inline before scaling by demography.
- `worker_contacts()` returns a 45-element `SVector` of within-sector contacts (already scaled by workforce). This is passed as `cw` to `Params` but is **not used** in the ODE.
- `consumer_worker_contacts()` returns a 45×4 `SMatrix` of consumer-worker contacts (scaled by demography). Not exposed as a parameter.

### DaedalusStructs.jl
- `Params` has a `cw::StaticVector` field that is unused in the ODE.
- No `cm_work`, `cm_cons_work`, or `openness` fields.

### Ode.jl
- Single FOI: `community_infectious = sum(Is .+ Ia * p.epsilon, dims=2)` then `foi = p.beta_now * p.contacts * community_infectious`
- No separation of workplace or consumer-worker routes.

### Model.jl / Events.jl
- NPIs currently modify `beta_now` (a scalar multiplier on total transmission).
- No `openness` concept or per-route NPI scaling.

---

## Proposed Changes

### 1. `Data.jl` — Separate contact matrices

**Refactor `prepare_contacts()`** to produce only the community contact matrix (without adding workplace or consumer-worker contacts). The new name `prepare_community_contacts()` makes the role explicit. Keep the old name as a deprecated alias or remove it.

```julia
# New: community-only 49×49 matrix
function prepare_community_contacts(cm=australia_contacts(); scaled=true)
    cm_x = ones(N_TOTAL_GROUPS, N_TOTAL_GROUPS) .* cm[i_WORKING_AGE, i_WORKING_AGE]
    cm_x[i_AGE_GROUPS, i_AGE_GROUPS] = cm
    cm_x[i_AGE_GROUPS, i_ECON_GROUPS] .= cm[:, i_WORKING_AGE]
    cm_x[i_ECON_GROUPS, i_AGE_GROUPS] .= reshape(cm[i_WORKING_AGE, :], 1, N_AGE_GROUPS)
    # NOTE: Do NOT add consumer_worker_contacts or worker_contacts here
    if scaled
        cm_x *= Diagonal(1 ./ prepare_demog())
    end
    return SMatrix{N_TOTAL_GROUPS, N_TOTAL_GROUPS}(cm_x)
end
```

**`worker_contacts()`** — already correct; returns a 45-element scaled vector. No change needed.

**`consumer_worker_contacts()`** — already correct; returns a 45×4 scaled matrix. No change needed. Export it.

**Update exports** to include `prepare_community_contacts` and `consumer_worker_contacts`.

> **Compatibility note**: `prepare_contacts()` used in `Helpers.jl` for beta/NGM calculation uses the unscaled matrix (`scaled=false`). The NGM should be computed from community contacts only, since workplace routes are relative to community beta. No change needed there, but verify after implementation.

---

### 2. `DaedalusStructs.jl` — Update `Params`

Replace the unused `cw` field with properly typed `cm_work` and `cm_cons_work`, and add an `openness` scalar.

```julia
mutable struct Params
    contacts::StaticArray         # community contacts (49×49), scaled by demography
    cm_work::StaticVector         # within-sector contacts (45-element), scaled by workers
    cm_cons_work::StaticMatrix    # consumer-worker contacts (45×4), scaled by demography
    ngm::StaticArray
    demography::StaticVector
    beta::Float64
    beta_now::Float64
    sigma::Float64
    p_sigma::Float64
    epsilon::Float64
    rho::Float64
    eta::Vector{Float64}
    omega::Vector{Float64}
    omega_now::Vector{Float64}
    gamma_Ia::Float64
    gamma_Is::Float64
    gamma_H::Vector{Float64}
    nu::Float64
    psi::Float64
    openness::Float64             # NPI scalar on workplace contacts [0,1]; 1 = fully open
    size::Int
end
```

The `SMatrix` type parameter for `cm_cons_work` will be `SMatrix{N_ECON_GROUPS, N_AGE_GROUPS, Float64}`.

---

### 3. `Ode.jl` — Three-component force of infection

Replace the single FOI with three separate components following the C++ implementation.

**Core ODE changes** (replacing lines 40-44 of current `Ode.jl`):

```julia
# Sum over vaccine strata to get total infectious per group
infectious = sum(Is .+ Ia * p.epsilon, dims=2)  # (N_TOTAL_GROUPS, 1)

# --- FOI Component 1: Community contacts (all groups) ---
foi_comm = p.beta_now * p.contacts * infectious  # (N_TOTAL_GROUPS, 1)
# Broadcast over vaccine strata
foi_comm_bcast = repeat(foi_comm, 1, N_VACCINE_STRATA)  # (N_TOTAL_GROUPS, N_VAX_STRATA)

# New infections from community route
new_I = S .* foi_comm_bcast  # (N_TOTAL_GROUPS, N_VAX_STRATA)

# --- FOI Component 2: Workplace within-sector contacts (workers only) ---
infectious_workers = infectious[(N_AGE_GROUPS + 1):end, :]  # (N_ECON_GROUPS, 1)
# cm_work is (N_ECON_GROUPS,) vector, element-wise multiply
foi_work = p.beta_now .* p.cm_work .* infectious_workers[:,1]  # (N_ECON_GROUPS,)
# Apply openness² (sector closure effect)
foi_work_scaled = foi_work .* (p.openness^2)

# --- FOI Component 3: Consumer-worker contacts (workers only) ---
infectious_age = infectious[1:N_AGE_GROUPS, 1]  # (N_AGE_GROUPS,)
foi_cw = p.beta_now .* (p.cm_cons_work * infectious_age)  # (N_ECON_GROUPS,)
# Apply openness¹
foi_cw_scaled = foi_cw .* p.openness

# --- Add workplace + consumer-worker infections for workers ---
# foi_work_scaled and foi_cw_scaled are (N_ECON_GROUPS,) vectors
# broadcast over vaccine strata and add to new_I for worker rows
worker_extra_foi = (foi_work_scaled .+ foi_cw_scaled)  # (N_ECON_GROUPS,)
worker_extra_foi_bcast = repeat(worker_extra_foi, 1, N_VACCINE_STRATA)
worker_rows = (N_AGE_GROUPS + 1):N_TOTAL_GROUPS
S_workers = @view S[worker_rows, :]
new_I[worker_rows, :] .+= S_workers .* worker_extra_foi_bcast
```

The rest of the ODE (dS, dE, dIs, dIa, dH, dR, dD) remains unchanged, as they depend only on `new_I`.

---

### 4. `Model.jl` — Pass new parameters

Update `daedalus()` to:
- Call `prepare_community_contacts()` instead of `prepare_contacts()` for the runtime contact matrix
- Compute `cm_work` via `worker_contacts()` (already called but passed as unused `cw`)
- Compute `cm_cons_work` via `consumer_worker_contacts()`
- Pass `openness=1.0` as default (fully open, no NPI)
- Construct `Params` with the new fields

**Function signature** — add `openness::Float64=1.0` as a parameter or leave it fully managed by NPIs.

**Params construction** change:
```julia
# old
parameters = Params(contacts, ngm, demog, cw, beta, beta, ...)
# new
parameters = Params(
    prepare_community_contacts(),
    worker_contacts(),
    consumer_worker_contacts(),
    ngm, demog,
    beta, beta,
    ...,
    openness=1.0,
    size
)
```

---

### 5. `Events.jl` — NPI openness handling

Currently the `Npi` and `TimedNpi` callbacks set `beta_now = beta * coef`. With the new structure, two options exist:

**Option A (minimal change)**: Keep modifying `beta_now`. Workplace and community contacts are then scaled equally. This is less faithful to the R package but requires no Events changes.

**Option B (faithful)**: Modify `openness` in the callbacks instead of (or in addition to) `beta_now`. The NPI coefficient becomes the openness scalar, and:
- Community transmission: unchanged (full `beta_now`)
- Workplace: scaled by `openness²`
- Consumer-worker: scaled by `openness`

Option B requires updating the callback effect functions in `make_param_changer`, `make_timed_npi_callbacks` to target `openness` instead of `beta_now`. This better matches the R package semantics where NPIs represent sector openness (e.g., 0.0 = full closure of workplace contacts, 1.0 = fully open).

**Recommendation**: Implement Option B. The `Params` struct already supports arbitrary field mutation via `setproperty!`, so the callbacks only need to target `:openness`.

---

## Files to Modify

| File | Change |
|------|--------|
| `src/Data.jl` | Add `prepare_community_contacts()`, export `consumer_worker_contacts()` |
| `src/DaedalusStructs.jl` | Replace `cw` with `cm_work`, `cm_cons_work`; add `openness` to `Params` |
| `src/Ode.jl` | Implement three-component FOI |
| `src/Model.jl` | Pass new parameters to `Params`; call new data functions |
| `src/Events.jl` | (Option B) Update NPI callbacks to set `openness` |

---

## Invariants and Constraints

- **State vector layout**: No change. `u` is still `(N_TOTAL_GROUPS, N_COMPARTMENTS, N_VACCINE_STRATA)`.
- **Worker rows**: `(N_AGE_GROUPS+1):N_TOTAL_GROUPS` = rows 5:49 of the state.
- **`cm_work` dimensions**: `SVector{N_ECON_GROUPS, Float64}` (45 elements), already scaled by workforce.
- **`cm_cons_work` dimensions**: `SMatrix{N_ECON_GROUPS, N_AGE_GROUPS, Float64}` (45×4), already scaled by age-group demography.
- **`openness` range**: `[0.0, 1.0]`. Default `1.0` = no NPI effect on workplace contacts.
- **NGM / beta calculation**: Use community contacts only (`prepare_community_contacts(scaled=false)`) to maintain R₀ calibration. This ensures beta calibration is unchanged by the new workplace routes being explicit rather than blended.

---

## Open Questions

1. **NGM validity**: Is the community-only contact matrix the correct input for R₀/beta calibration, or should the NGM include the expected contribution from workplace contacts? The R package passes a separate NGM pre-computed from community contacts, suggesting community-only is correct.

2. **Backward compatibility of `prepare_contacts()`**: Should the old function be kept as a deprecated alias returning the blended matrix? Users calling it directly would get different results from the ODE which now uses `prepare_community_contacts()`.

3. **Susceptibility scaling**: The C++ implementation multiplies `new_I` by a `susc` matrix to reduce vaccinated infections. The Julia package currently does not implement this step. Consider adding this as a follow-up.

4. **NPI openness vs. beta**: If Option B is adopted for Events, the semantics of existing `Npi` and `TimedNpi` `coefs` change from "fraction of beta" to "openness level". This should be clearly documented.
