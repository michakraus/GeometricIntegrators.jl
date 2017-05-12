
abstract AbstractNonlinearFunctionParametersVPRK{DT,TT,ΑT,FT,D,S} <: NonlinearFunctionParameters{DT}
abstract AbstractIntegratorVPRK{DT,TT} <: Integrator{DT,TT}

immutable NonlinearFunctionCacheVPRK{ST}
    Q::Matrix{ST}
    V::Matrix{ST}
    P::Matrix{ST}
    F::Matrix{ST}

    v::Array{ST,1}
    y::Array{ST,1}
    z::Array{ST,1}

    function NonlinearFunctionCacheVPRK(D,S)
        # create internal stage vectors
        Q = zeros(ST,D,S)
        V = zeros(ST,D,S)
        P = zeros(ST,D,S)
        F = zeros(ST,D,S)

        # create update vectors
        v = zeros(ST,D)
        y = zeros(ST,D)
        z = zeros(ST,D)

        new(Q, V, P, F, v, y, z)
    end
end

immutable NonlinearFunctionCacheVPRKprojection{ST}
    q̅::Vector{ST}
    p̅::Vector{ST}
    λ::Vector{ST}

    U::Array{ST,2}
    G::Array{ST,2}

    u::Array{ST,1}
    g::Array{ST,1}

    function NonlinearFunctionCacheVPRKprojection(D)
        # create projected solution vectors
        q̅ = zeros(ST,D)
        p̅ = zeros(ST,D)
        λ = zeros(ST,D)

        # create projection vectors
        U = zeros(ST,D,2)
        G = zeros(ST,D,2)

        # create update vectors
        u = zeros(ST,D)
        g = zeros(ST,D)

        new(q̅, p̅, λ, U, G, u, g)
    end
end


function compute_stages_vprk!(x, Q, V, P, F, params)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages_vprk!(x, q̅, p̅, λ, Q, V, U, P, F, G, params)
    # copy x to V
    compute_stages_v_vprk!(x, V, params)

    # compute U, G and p̅
    compute_projection_vprk!(x, q̅, p̅, λ, V, U, G, params)

    # compute Q
    compute_stages_q_vprk!(Q, V, U, params)

    # compute P and F
    compute_stages_p_vprk!(Q, V, P, F, params)
end

function compute_stages_v_vprk!{ST,DT,TT,AT,FT,D,S}(x::Vector{ST}, V::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S})
    @assert D == size(V,1)
    @assert S == size(V,2)

    # copy x to V
    for i in 1:S
        for k in 1:D
            V[k,i] = x[D*(i-1)+k]
        end
    end
end

function compute_stages_q_vprk!{ST,DT,TT,AT,FT,D,S}(Q::Matrix{ST}, V::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S})
    @assert D == size(Q,1) == size(V,1)
    @assert S == size(Q,2) == size(V,2)

    local y::ST

    # compute Q
    for i in 1:S
        for k in 1:D
            y = 0
            for j in 1:S
                y += params.t_q.a[i,j] * V[k,j]
            end
            Q[k,i] = params.q[k] + params.Δt * y
        end
    end
end


function compute_stages_q_vprk!{ST,DT,TT,AT,FT,D,S}(Q::Matrix{ST}, V::Matrix{ST}, U::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S})
    @assert D == size(Q,1) == size(V,1) == size(U,1)
    @assert S == size(Q,2) == size(V,2)

    local y::ST

    # compute Q
    for i in 1:S
        for k in 1:D
            y = params.R[1] * U[k,1]
            for j in 1:S
                y += params.t_q.a[i,j] * V[k,j]
            end
            Q[k,i] = params.q[k] + params.Δt * y
        end
    end
end


@generated function compute_stages_p_vprk!{ST,DT,TT,AT,FT,D,S}(Q::Matrix{ST}, V::Matrix{ST}, P::Matrix{ST}, F::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S})
    # create temporary vectors
    tQ = zeros(ST,D)
    tV = zeros(ST,D)
    tP = zeros(ST,D)
    tF = zeros(ST,D)

    compute_stages_vprk = quote
        @assert D == size(Q,1) == size(V,1) == size(P,1) == size(F,1)
        @assert S == size(Q,2) == size(V,2) == size(P,2) == size(F,2)

        local tᵢ::TT

        # compute P=α(Q) and F=f(Q)
        for i in 1:S
            tᵢ = params.t + params.Δt * params.t_q.c[i]
            simd_copy_xy_first!($tQ, Q, i)
            simd_copy_xy_first!($tV, V, i)
            params.α(tᵢ, $tQ, $tV, $tP)
            params.f(tᵢ, $tQ, $tV, $tF)
            simd_copy_yx_first!($tP, P, i)
            simd_copy_yx_first!($tF, F, i)
        end
    end

    return compute_stages_vprk
end


@generated function compute_rhs_vprk!{ST,DT,TT,AT,FT,D,S}(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S})
    compute_stages_vprk = quote
        local z::ST

        # compute b = - [(P-AF)]
        for i in 1:S
            for k in 1:D
                z = 0
                for j in 1:S
                    z += params.t_p.a[i,j] * F[k,j]
                end
                b[D*(i-1)+k] = - (P[k,i] - params.p[k]) + params.Δt * z
            end
        end
    end

    return compute_stages_vprk
end


@generated function compute_rhs_vprk!{ST,DT,TT,AT,FT,D,S}(b::Vector{ST}, P::Matrix{ST}, F::Matrix{ST}, G::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S})
    compute_stages_vprk = quote
        local z::ST

        # compute b = - [(P-AF)]
        for i in 1:S
            for k in 1:D
                z = params.R[1] * G[k,1]
                for j in 1:S
                    z += params.t_p.a[i,j] * F[k,j]
                end
                b[D*(i-1)+k] = - (P[k,i] - params.p[k]) + params.Δt * z
            end
        end
    end

    return compute_stages_vprk
end


function compute_rhs_vprk_projection_q!{ST,DT,TT,AT,FT,D,S}(b::Vector{ST}, q̅::Vector{ST}, V::Matrix{ST}, U::Matrix{ST}, offset::Int, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S})
    local y::ST

    for k in 1:D
        y = params.R[1] * U[k,1] + params.R[2] * U[k,2]
        for j in 1:S
            y += params.t_q.b[j] * V[k,j]
        end
        b[offset+k] = - (q̅[k] - params.q[k]) + params.Δt * y
    end
end


function compute_rhs_vprk_projection_p!{ST,DT,TT,AT,FT,D,S}(b::Vector{ST}, p̅::Vector{ST}, F::Matrix{ST}, G::Matrix{ST}, offset::Int, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S})
    local z::ST

    for k in 1:D
        z = params.R[1] * G[k,1] + params.R[2] * G[k,2]
        for j in 1:S
            z += params.t_p.b[j] * F[k,j]
        end
        b[offset+k] = - (p̅[k] - params.p[k]) + params.Δt * z
    end
end


@generated function compute_rhs_vprk_correction!{ST,DT,TT,AT,FT,D,S}(b::Vector{ST}, V::Matrix{ST}, params::AbstractNonlinearFunctionParametersVPRK{DT,TT,AT,FT,D,S})
    μ = zeros(ST,D)

    compute_stages_vprk = quote
        local sl::Int = div(S+1, 2)

        if length(params.d_v) > 0
            # compute μ
            for k in 1:D
                $μ[k] = params.t_p.b[sl] / params.d_v[sl] * b[D*(sl-1)+k]
            end

            # replace equation for Pₗ with constraint on V
            for k in 1:D
                b[D*(sl-1)+k] = 0
                for i in 1:S
                    b[D*(sl-1)+k] += V[k,i] * params.d_v[i]
                end
            end

            # modify P₁, ..., Pₛ except for Pₗ
            for i in 1:S
                if i ≠ sl
                    z = params.d_v[i] / params.t_p.b[i]
                    for k in 1:D
                        b[D*(i-1)+k] -= z * $μ[k]
                    end
                end
            end
        end
    end

    return compute_stages_vprk
end


function update_solution!{DT,TT}(int::AbstractIntegratorVPRK{DT,TT}, cache::NonlinearFunctionCacheVPRK{DT})
    simd_mult!(cache.y, cache.V, int.tableau.q.b)
    simd_mult!(cache.z, cache.F, int.tableau.p.b)
    simd_axpy!(int.Δt, cache.y, int.q, int.qₑᵣᵣ)
    simd_axpy!(int.Δt, cache.z, int.p, int.pₑᵣᵣ)
end


function project_solution!{DT,TT}(int::AbstractIntegratorVPRK{DT,TT}, cache::NonlinearFunctionCacheVPRKprojection{DT}, R::Vector{TT})
    simd_mult!(cache.u, cache.U, R)
    simd_mult!(cache.g, cache.G, R)
    simd_axpy!(int.Δt, cache.u, int.q, int.qₑᵣᵣ)
    simd_axpy!(int.Δt, cache.g, int.p, int.pₑᵣᵣ)
end

function cut_periodic_solution!(int::AbstractIntegratorVPRK)
    # take care of periodic solutions
    for k in 1:int.equation.d
        if int.equation.periodicity[k] ≠ 0
            int.q[k] = mod(int.q[k], int.equation.periodicity[k])
        end
    end
end