
module Helpers

export get_beta

using ..Constants

using LinearAlgebra

function get_beta(cm, r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is)
    n_groups = size(cm)[1]
    sigma_1 = sigma * (1.0 - p_sigma)
    sigma_2 = sigma * p_sigma

    foi_a = epsilon * cm
    foi_s = cm

    f_mat = hcat(
        zeros(n_groups, n_groups),
        foi_a, foi_s
    )

    vvec = [
        repeat([sigma], n_groups);
        repeat([gamma_Ia], n_groups);
        repeat([gamma_Is], n_groups)
    ]

    v_mat = Matrix(Diagonal(vvec))
    v_mat[(1:n_groups).+n_groups, 1:n_groups] =
        Matrix(Diagonal(repeat([-sigma_1], n_groups)))
    v_mat[(1:n_groups).+n_groups*2, 1:n_groups] =
        Matrix(Diagonal(repeat([-sigma_2], n_groups)))

    v_inv = inv(v_mat)

    ngm = f_mat * v_inv[:, 1:n_groups]

    r0a = maximum(real(eigen(ngm).values))

    beta = r0 / r0a

    return beta
end

end