
function getTableauVSPARKMidpointProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}
    @assert q.s == p.s
    s = q.s

    g = getCoefficientsGLRK(1)
    α̃ = g.a
    β = g.b
    γ = g.c

    α = SMatrix{s, 1, T}(0.5 .* ones(T, s, 1))

    _q_ã = zeros(T, 1, s)
    for i in 1:s
        _q_ã[1,i] = q.b[i] / β[1] * ( β[1] - α[i,1] )
    end
    q_ã = SMatrix{1, s, T}(_q_ã)

    _p_ã = zeros(T, 1, s)
    for i in 1:s
        _p_ã[1,i] = p.b[i] / β[1] * ( β[1] - α[i,1] )
    end
    p_ã = SMatrix{1, s, T}(_p_ã)

    β = @SVector [α̃[1,1] * (1 + R∞)]

    d = ones(T, 1)

    ω  = zeros(T, 1, 2)
    ω .= [0.0  1.0]
    δ  = zeros(T, 0, 1)

    return TableauVSPARKprimary(name, min(q.o, p.o),
                        q.a, p.a, α, α,
                        q_ã, p_ã, α̃, α̃,
                        q.b, p.b, β, β,
                        q.c, p.c, γ, d,
                        ω, δ)
end

"Tableau for Gauss-Legendre method with s stages and midpoint projection."
function getTableauVSPARKGLRKpMidpoint(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKMidpointProjection(Symbol("vsparkglrk", s, "pMidpoint"), glrk, glrk; R∞=(-1)^s)
end



function getTableauVSPARKSymmetricProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

    @assert q.s == p.s

    o = min(q.o, p.o)

    a_q = q.a
    a_p = p.a

    _α_q = zeros(T, q.s, 2)
    _α_q[:,1] .= 0.5
    α_q = SMatrix{q.s, 2, T}(_α_q)

    _α_p = zeros(T, p.s, 2)
    _α_p[:,1] .= 0.5
    α_p = SMatrix{p.s, 2, T}(_α_p)

    a_q̃ = SMatrix(transpose(hcat(zero(q.b), q.b)))
    a_p̃ = SMatrix(transpose(hcat(zero(p.b), p.b)))

    α_q̃ = @SMatrix [0.0  0.0
                    0.5  R∞*0.5]
    α_p̃ = @SMatrix [0.0  0.0
                    0.5  R∞*0.5]

    b_q = q.b
    b_p = p.b
    β_q = @SVector [0.5, R∞*0.5]
    β_p = @SVector [0.5, R∞*0.5]

    c_q = q.c
    c_p = p.c
    c_λ = @SVector [ 0.0, 1.0]
    d_λ = @SVector [ 0.5, 0.5]

    ω_λ  = zeros(T, 1, 3)
    ω_λ .= [0.5 R∞*0.5 0.0]

    δ_λ  = zeros(T, 1, 2)
    δ_λ .= [-1.0 +1.0]


    if length(d) == 0
        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVSPARKprimary(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end


"Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB2pSymmetric()
    getTableauVSPARKSymmetricProjection(:LobIIIAIIIB2pSymmetric, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), get_lobatto_d_vector(2); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB3pSymmetric()
    getTableauVSPARKSymmetricProjection(:LobIIIAIIIB3pSymmetric, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), get_lobatto_d_vector(3); R∞=+1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with four stages and symmetric projection."
function getTableauVSPARKLobIIIAIIIB4pSymmetric()
    getTableauVSPARKSymmetricProjection(:LobIIIAIIIB4pSymmetric, getCoefficientsLobIIIA4(), getCoefficientsLobIIIB4(), get_lobatto_d_vector(4); R∞=-1)
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauVSPARKGLRKpSymmetric(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKSymmetricProjection(Symbol("vpglrk", s, "pSymmetric"), glrk, glrk; R∞=(-1)^s)
end


function getTableauVSPARKSymplecticProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=+1) where {T}

    la = getCoefficientsLobIIIA2()
    lb = getCoefficientsLobIIIB2()

    @assert q.s == p.s
    @assert q.b == p.b
    @assert q.c == p.c

    @assert la.s == lb.s
    @assert la.b == lb.b
#    @assert la.c == lb.c

    s = q.s
    σ = la.s
    ρ = 0

    a_q = q.a
    a_p = p.a
    b_q = q.b
    b_p = p.b

    # β_q = la.b
    # β_p = lb.b
    β_q = @SVector [0.5, R∞*0.5]
    β_p = @SVector [0.5, R∞*0.5]

    _α_q = zeros(T, s, 2)
    _α_q[:,1] .= 0.5
    α_q = SMatrix{s, 2, T}(_α_q)

    _α_p = zeros(T, s, 2)
    _α_p[:,1] .= 0.5
    α_p = SMatrix{s, 2, T}(_α_q)

#    α_q = zeros(T, s, σ)
#    α_p = zeros(T, s, σ)
#    for i in 1:s
#        for j in 1:σ
#            α_q[i,j] = #q.b[i] / β[1] * ( β[1] - α[i,1] )
#            α_p[i,j] = #p.b[i] / β[1] * ( β[1] - α[i,1] )
#        end
#    end

    _a_q̃ = zeros(T, σ, s)
    _a_p̃ = zeros(T, σ, s)
    for i in 1:σ
        for j in 1:s
            _a_q̃[i,j] = b_q[j] / β_p[i] * (β_p[i] - α_p[j,i])
            _a_p̃[i,j] = b_p[j] / β_q[i] * (β_q[i] - α_q[j,i])
        end
    end
    a_q̃ = SMatrix{σ, s, T}(_a_q̃)
    a_p̃ = SMatrix{σ, s, T}(_a_p̃)

    α_q̃ = la.a
    α_p̃ = lb.a

    c_q = q.c
    c_p = p.c
    c_λ = la.c
    d_λ = @SVector [0.5, 0.5]


    ω_λ = @SMatrix [0.5 0.5 0.0
                    0.0 0.0 1.0]

    δ_λ = @SMatrix zeros(T, ρ, σ)


    if length(d) == 0
        return TableauVSPARKprimary(name, min(q.o, p.o),
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ)
    else
        @assert length(d) == s

        return TableauVSPARKprimary(name, min(q.o, p.o),
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            ω_λ, δ_λ, d)
    end

end

function getTableauVSPARKGLRKpSymplectic(s)
    glrk = getCoefficientsGLRK(s)
    getTableauVSPARKSymplecticProjection(Symbol("VSPARK", s, "pSymplectic"), glrk, glrk; R∞=(-1)^s)
end
