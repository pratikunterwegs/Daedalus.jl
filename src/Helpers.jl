
module Helpers

export get_beta

using ..Constants

using LinearAlgebra

function get_beta(cm, r0, sigma, p_sigma, epsilon, gamma_Ia, gamma_Is)
    sigma_1 = sigma * (1.0 - p_sigma)
    sigma_2 = sigma * p_sigma

    foi_a = epsilon * cm
    foi_s = cm

    f_mat = hcat(
        zeros(N_AGE_GROUPS, N_AGE_GROUPS),
        foi_a, foi_s
    )

    vvec = [
        repeat([sigma], N_AGE_GROUPS);
        repeat([gamma_Ia], N_AGE_GROUPS);
        repeat([gamma_Is], N_AGE_GROUPS)
    ]

    v_mat = Matrix(Diagonal(vvec))
    v_mat[(1:N_AGE_GROUPS).+N_AGE_GROUPS, 1:N_AGE_GROUPS] =
        Matrix(Diagonal(repeat([-sigma_1], N_AGE_GROUPS)))
    v_mat[(1:N_AGE_GROUPS).+N_AGE_GROUPS*2, 1:N_AGE_GROUPS] =
        Matrix(Diagonal(repeat([-sigma_2], N_AGE_GROUPS)))

    v_inv = inv(v_mat)

    ngm = f_mat * v_inv[:, 1:N_AGE_GROUPS]

    r0a = maximum(real(eigen(ngm).values))

    beta = r0 / r0a

    return beta
end

end