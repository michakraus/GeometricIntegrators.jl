"Holds the tableau of an Specialised Partitioned Additive Runge-Kutta method."
struct AbstractTableauSPARK{IT, DT <: Number, S, R, Ρ} <: AbstractTableau{DT}
    name::Symbol
    o::Int
    s::Int
    r::Int
    ρ::Int

    q::CoefficientsARK{DT,S,R}
    p::CoefficientsARK{DT,S,R}

    q̃::CoefficientsPRK{DT,S,R}
    p̃::CoefficientsPRK{DT,S,R}

    λ::CoefficientsMRK{DT,R}

    ω::Matrix{DT}
    δ::Matrix{DT}
    d::Vector{DT}

    function AbstractTableauSPARK{IT,DT}(name::Symbol, o::Int, s::Int, r::Int, ρ::Int, q, p, q̃, p̃, λ, ω, δ, d=DT[]) where {IT, DT}
        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"
        @assert ρ ≥ 0 && ρ ≤ r+1

        @assert s==q.s==p.s==q̃.s==p̃.s
        @assert r==q.r==p.r==q̃.r==p̃.r==λ.r
        @assert size(ω, 1)==r-ρ
        # @assert size(ω, 2)==r+1
        @assert size(δ, 1)==ρ

        @assert length(d)==0 || length(d)==s

        new{IT,DT,s,r,ρ}(name, o, s, r, ρ, q, p, q̃, p̃, λ, ω, δ, d)
    end

    function AbstractTableauSPARK{IT}(name::Symbol, order::Int,
                             a_q::AbstractMatrix{DT}, a_p::AbstractMatrix{DT},
                             α_q::AbstractMatrix{DT}, α_p::AbstractMatrix{DT},
                             a_q̃::AbstractMatrix{DT}, a_p̃::AbstractMatrix{DT},
                             α_q̃::AbstractMatrix{DT}, α_p̃::AbstractMatrix{DT},
                             b_q::AbstractVector{DT}, b_p::AbstractVector{DT},
                             β_q::AbstractVector{DT}, β_p::AbstractVector{DT},
                             c_q::AbstractVector{DT}, c_p::AbstractVector{DT},
                             c_λ::AbstractVector{DT}, d_λ::AbstractVector{DT},
                             ω::AbstractMatrix{DT}, δ::AbstractMatrix{DT}, d::AbstractVector{DT}=DT[]) where {IT, DT <: Number}

        s = length(c_q)
        r = length(c_λ)
        ρ = size(δ, 1)

        @assert s > 0 "Number of stages s must be > 0"
        @assert r > 0 "Number of stages r must be > 0"

        @assert s==size(a_q,1)==size(a_q,2)==length(b_q)==length(c_q)
        @assert s==size(a_p,1)==size(a_p,2)==length(b_p)==length(c_p)
        @assert s==size(α_q,1)==size(α_p,1)
        @assert r==size(α_q,2)==size(α_p,2)
        @assert r==length(c_λ)==length(d_λ)
        @assert r==size(a_q̃,1)==size(a_p̃,1)
        @assert s==size(a_q̃,2)==size(a_p̃,2)
        @assert r==size(α_q̃,1)==size(α_q̃,2)==length(β_q)
        @assert r==size(α_p̃,1)==size(α_p̃,2)==length(β_p)

        q = CoefficientsARK{DT}(name, order, s, r, a_q, b_q, c_q, α_q, β_q)
        p = CoefficientsARK{DT}(name, order, s, r, a_p, b_p, c_p, α_p, β_p)
        q̃ = CoefficientsPRK{DT}(name, order, s, r, a_q̃, c_λ, α_q̃)
        p̃ = CoefficientsPRK{DT}(name, order, s, r, a_p̃, c_λ, α_p̃)
        λ = CoefficientsMRK{DT}(name, r, d_λ, c_λ)

        AbstractTableauSPARK{IT,DT}(name, order, s, r, ρ, q, p, q̃, p̃, λ, ω, δ, d)
    end

    function AbstractTableauSPARK{IT}(name::Symbol, order::Int,
                             a_q::AbstractMatrix{DT}, a_p::AbstractMatrix{DT},
                             α_q::AbstractMatrix{DT}, α_p::AbstractMatrix{DT},
                             a_q̃::AbstractMatrix{DT}, a_p̃::AbstractMatrix{DT},
                             α_q̃::AbstractMatrix{DT}, α_p̃::AbstractMatrix{DT},
                             b_q::AbstractVector{DT}, b_p::AbstractVector{DT},
                             β_q::AbstractVector{DT}, β_p::AbstractVector{DT},
                             c_q::AbstractVector{DT}, c_p::AbstractVector{DT},
                             c_λ::AbstractVector{DT}, d_λ::AbstractVector{DT},
                             d::AbstractVector{DT}=DT[]) where {IT, DT <: Number}

        R = length(c_λ)
        AbstractTableauSPARK{IT}(name, order,
                                 a_q, a_p, α_q, α_p,
                                 a_q̃, a_p̃, α_q̃, α_p̃,
                                 b_q, b_p, β_q, β_p,
                                 c_q, c_p, c_λ, d_λ,
                                 hcat(Array(Diagonal(ones(DT,R))), zeros(DT,R)),
                                 zeros(DT,0,R), d)
    end
end

# TODO function readAbstractTableauSPARKFromFile(dir::AbstractString, name::AbstractString)


"Holds the tableau of an Specialised Partitioned Additive Runge-Kutta method."
const TableauSPARK = AbstractTableauSPARK{:spark}
