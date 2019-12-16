
"Parameters for right-hand side function of Specialised Partitioned Additive Runge-Kutta methods for Variational systems."
mutable struct AbstractParametersVSPARK{IT,DT,TT,D,S,R,P,FT,PT,UT,GT,ϕT,tabType} <: Parameters{DT,TT}
    f_f::FT
    f_p::PT
    f_u::UT
    f_g::GT
    f_ϕ::ϕT

    Δt::TT

    tab::tabType

    t::TT
    q::Vector{DT}
    p::Vector{DT}
    λ::Vector{DT}

    function AbstractParametersVSPARK{IT,DT,D,S,R,P}(f_f::FT, f_p::PT, f_u::UT, f_g::GT, f_ϕ::ϕT, Δt::TT, tab::tabType) where {IT,DT,TT,D,S,R,P,FT,PT,UT,GT,ϕT,tabType}
        # create solution vectors
        q = zeros(DT,D)
        p = zeros(DT,D)
        λ = zeros(DT,D)

        new{IT,DT,TT,D,S,R,P,FT,PT,UT,GT,ϕT,tabType}(f_f, f_p, f_u, f_g, f_ϕ, Δt, tab, zero(TT), q, p, λ)
    end
end


function update_params!(params::AbstractParametersVSPARK, cache::IntegratorCacheSPARK)
    # set time for nonlinear solver and copy previous solution
    params.t  = cache.t
    params.q .= cache.q
    params.p .= cache.p
    params.λ .= cache.λ
end


equation(int::AbstractIntegratorVSPARK) = int.equation
timestep(int::AbstractIntegratorVSPARK) = int.params.Δt
tableau(int::AbstractIntegratorVSPARK) = int.tableau
nstages(int::AbstractIntegratorVSPARK) = int.tableau.s
pstages(int::AbstractIntegratorVSPARK) = int.tableau.r


function initial_guess!(int::AbstractIntegratorVSPARK, cache::IntegratorCacheSPARK)
    for i in 1:nstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q.c[i], tableau(int).p.c[i])

        for k in 1:ndims(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
            int.solver.x[3*(ndims(int)*(i-1)+k-1)+3] = cache.ṽ[k]
        end
    end

    for i in 1:pstages(int)
        evaluate!(int.iguess, cache.q, cache.p, cache.v, cache.f,
                              cache.q̅, cache.p̅, cache.v̅, cache.f̅,
                              cache.q̃, cache.p̃, cache.ṽ, cache.f̃,
                              tableau(int).q̃.c[i], tableau(int).p̃.c[i])

        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+1] = (cache.q̃[k] - cache.q[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+2] = (cache.p̃[k] - cache.p[k])/timestep(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(ndims(int)*(i-1)+k-1)+3] = 0
        end
    end

    if isdefined(tableau(int), :λ) && tableau(int).λ.c[1] == 0
        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*(k-1)+3] = cache.λ[k]
        end
    end

    if length(tableau(int).d) > 0
        for k in 1:ndims(int)
            int.solver.x[3*ndims(int)*nstages(int)+3*ndims(int)*pstages(int)+k] = 0
        end
    end
end