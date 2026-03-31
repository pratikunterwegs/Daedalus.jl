
module Helpers

export get_beta, get_ngm, sum_by_age, dominant_eigenvalue, weighted_slice_sum!

using ..Constants
import ..DataLoader as DL

using LinearAlgebra

"""
    get_beta(cm, r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is)::Union{Float64, Vector{Float64}}

Calculate transmission rate(s) `beta` given a contact matrix and epidemiological parameters.

Uses the next-generation matrix (NGM) method to compute beta values that achieve
the target basic reproduction number(s) R₀.

# Arguments
- `cm::Matrix{Float64}`: Contact matrix (n_groups × n_groups)
- `r0::Union{Float64, Vector{Float64}}`: Target basic reproduction number(s)
- `sigma::Float64`: Rate of progression from exposed to infectious
- `p_sigma::Float64`: Proportion of exposed who develop symptoms
- `epsilon::Float64`: Reduction factor for asymptomatic transmission
- `gamma_Ia::Float64`: Recovery rate for asymptomatic infectious
- `gamma_Is::Float64`: Recovery rate for symptomatic infectious

# Returns
- If `r0` is a scalar: returns a single `Float64` value for beta
- If `r0` is a vector: returns a `Vector{Float64}` of beta values
"""
function get_beta(cm::Matrix{Float64}, r0::Union{Float64, Vector{Float64}},
        sigma::Float64, p_sigma::Float64, epsilon::Float64,
        gamma_Ia::Float64, gamma_Is::Float64)::Union{Float64, Vector{Float64}}
    n_groups = size(cm)[1]
    sigma_1 = sigma * (1.0 - p_sigma)
    sigma_2 = sigma * p_sigma

    foi_a = epsilon * cm
    foi_s = cm

    f_mat = hcat(
        zeros(n_groups, n_groups),
        foi_a, foi_s
    )

    vvec = [repeat([sigma], n_groups);
            repeat([gamma_Ia], n_groups);
            repeat([gamma_Is], n_groups)]

    v_mat = Matrix(Diagonal(vvec))
    v_mat[(1:n_groups) .+ n_groups, 1:n_groups] = Matrix(Diagonal(repeat(
        [-sigma_1], n_groups)))
    v_mat[(1:n_groups) .+ n_groups * 2, 1:n_groups] = Matrix(Diagonal(repeat(
        [-sigma_2], n_groups)))

    v_inv = inv(v_mat)

    ngm = f_mat * v_inv[:, 1:n_groups]

    r0a = maximum(real(eigen(ngm).values))

    beta = r0 ./ r0a

    return beta
end

function get_beta(cm::Matrix{Float64}, infdata::DL.InfectionData)
    r0 = infdata.r0
    sigma = infdata.sigma
    p_sigma = infdata.p_sigma
    epsilon = infdata.epsilon
    gamma_Ia = infdata.gamma_Ia
    gamma_Is = infdata.gamma_Is

    get_beta(cm, r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is)
end

"""
    get_ngm(cm, beta, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is)::Union{Matrix{Float64}, Vector{Matrix{Float64}}}

Calculate the next-generation matrix (NGM) given a contact matrix and transmission rate(s).

The NGM is used to compute the effective reproduction number Rt accounting for
current susceptibility levels in the population.

# Arguments
- `cm::Matrix{Float64}`: Contact matrix (n_groups × n_groups)
- `beta::Union{Float64, Vector{Float64}}`: Transmission rate(s)
- `sigma::Float64`: Rate of progression from exposed to infectious
- `p_sigma::Float64`: Proportion of exposed who develop symptoms
- `epsilon::Float64`: Reduction factor for asymptomatic transmission
- `gamma_Ia::Float64`: Recovery rate for asymptomatic infectious
- `gamma_Is::Float64`: Recovery rate for symptomatic infectious

# Returns
- If `beta` is a scalar: returns a single `Matrix{Float64}`
- If `beta` is a vector: returns a `Vector{Matrix{Float64}}` of NGMs
"""
function get_ngm(cm::Matrix{Float64}, beta::Float64,
        sigma::Float64, p_sigma::Float64, epsilon::Float64,
        gamma_Ia::Float64, gamma_Is::Float64)::Matrix{Float64}
    n_groups = size(cm)[1]
    sigma_1 = sigma * (1.0 - p_sigma)
    sigma_2 = sigma * p_sigma

    foi_a = epsilon * beta * cm
    foi_s = beta * cm

    f_mat = hcat(
        zeros(n_groups, n_groups),
        foi_a, foi_s
    )

    vvec = [repeat([sigma], n_groups);
            repeat([gamma_Ia], n_groups);
            repeat([gamma_Is], n_groups)]

    v_mat = Matrix(Diagonal(vvec))
    v_mat[(1:n_groups) .+ n_groups, 1:n_groups] = Matrix(Diagonal(repeat(
        [-sigma_1], n_groups)))
    v_mat[(1:n_groups) .+ n_groups * 2, 1:n_groups] = Matrix(Diagonal(repeat(
        [-sigma_2], n_groups)))

    v_inv = inv(v_mat)

    ngm = f_mat * v_inv[:, 1:n_groups]

    return ngm
end

"""
    get_ngm(cm::Matrix, beta::Vector, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is)::Vector{Matrix{Float64}}

Vector-dispatch version: compute one NGM for each transmission rate in `beta`.
"""
function get_ngm(cm::Matrix{Float64}, beta::Vector{Float64},
        sigma::Float64, p_sigma::Float64, epsilon::Float64,
        gamma_Ia::Float64, gamma_Is::Float64)::Vector{Matrix{Float64}}
    return [get_ngm(cm, beta_i, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is)
            for
            beta_i in beta]
end

"""
    sum_by_age(state, index)::Vector{Float64}

Sum state array over economic groups into age groups for a given compartment.

Takes the full state array and sums all economic sector rows into the working-age
group, producing a 4-element vector (one per age group).

# Arguments
- `state`: 3D state array (N_TOTAL_GROUPS × N_COMPARTMENTS × N_VACCINE_STRATA)
- `index`: Compartment index (integer)

# Returns
A `Vector{Float64}` of length 4 (one value per age group) with economic sectors
summed into the working-age group.
"""
function sum_by_age(state, index)::Vector{Float64}
    # subset indices
    state_ = @view state[:, index, :]
    state_ = sum(state_, dims = 2)

    # sum working age and add to index
    s_ = state_[1:N_AGE_GROUPS]
    s_[i_WORKING_AGE] += sum(last(state_, N_ECON_GROUPS))

    return s_
end

"""
    dominant_eigenvalue(A; v_init=nothing, max_iter=100, tol=1e-6)

Compute the dominant (largest magnitude) eigenvalue using power iteration.
This is significantly faster than computing all eigenvalues when only the
largest eigenvalue is needed.

# Arguments
- `A::AbstractMatrix`: The matrix whose dominant eigenvalue is to be computed
- `v_init::Union{Nothing,AbstractVector}`: Optional initial vector for iteration.
  If not provided, a random vector is used.
- `max_iter::Int`: Maximum number of iterations (default: 100)
- `tol::Float64`: Convergence tolerance (default: 1e-6)

# Returns
- The dominant eigenvalue (largest in magnitude)

# Algorithm
Uses the power iteration method: repeatedly multiplies a vector by the matrix
and normalizes. The Rayleigh quotient converges to the dominant eigenvalue.

# Performance
For n×n matrices, this method is O(kn²) where k is the number of iterations
(typically k << n), compared to O(n³) for full eigendecomposition.

# Example
```julia
A = rand(49, 49)
λ_dom = dominant_eigenvalue(A)
```
"""
function dominant_eigenvalue(A::AbstractMatrix;
        v_init::Union{Nothing, AbstractVector} = nothing,
        max_iter::Int = 100,
        tol::Float64 = 1e-6)
    n = size(A, 1)

    # Initialize with provided vector or random vector
    if isnothing(v_init)
        v = randn(n)
    else
        v = copy(v_init)
    end

    # Normalize initial vector
    v_norm = norm(v)
    if v_norm > 0
        v = v / v_norm
    else
        # Handle zero initial vector
        v = randn(n)
        v = v / norm(v)
    end

    λ = zero(eltype(A))

    for i in 1:max_iter
        # Apply matrix
        v_new = A * v

        # Compute Rayleigh quotient (eigenvalue estimate)
        λ_new = dot(v, v_new)

        # Normalize new vector
        v_new_norm = norm(v_new)
        if v_new_norm > 0
            v_new = v_new / v_new_norm
        else
            # Matrix annihilated the vector, eigenvalue is zero
            return zero(eltype(A))
        end

        # Check convergence
        if i > 1 && abs(λ_new - λ) < tol * max(abs(λ_new), 1.0)
            return λ_new
        end

        λ = λ_new
        v = v_new
    end

    # If we reach here, convergence was not achieved, but return best estimate
    return λ
end

"""
    weighted_slice_sum!(X, v, result)

Compute the weighted sum of 2D slices of a 3D array in-place, writing the
result into `result`.

Mathematically: `result[m, k] = Σᵢ X[m, k, i] * v[i]`

This is a tensor contraction over the third dimension of `X`, equivalent to a
matrix-vector product after reshaping. `X` is reinterpreted (without copying)
as a `(M*K) × N` matrix whose columns are the flattened slices `X[:,:,i]`;
BLAS `gemv` then computes the weighted sum in a single call.

# Arguments
- `X::Array{T,3}`: Three-dimensional input array of size `(M, K, N)`
- `v::Vector{T}`: Weight vector of length `N`; element `v[i]` scales slice `X[:,:,i]`
- `result::Array{T,2}`: Pre-allocated output array of size `(M, K)`; overwritten in-place

# Notes
- `reshape` and `vec` are O(1) view operations — no data is copied.
- Efficient for large `M`, `K`. For small `N` (e.g. `N_VACCINE_STRATA = 2`),
  the BLAS dispatch overhead may dominate; a plain `@inbounds @simd` loop can
  be faster in that regime.
"""
function weighted_slice_sum!(X::Array{T, 3}, v::Vector{T}, result::Array{T, 2}) where {T}
    M, K, N = size(X)
    mul!(vec(result), reshape(X, M * K, N), v)
end

end
