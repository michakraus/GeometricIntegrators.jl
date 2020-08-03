
function getTableauSymplecticProjection(name, q::CoefficientsRK{T}, p::CoefficientsRK{T}, d=[]; R∞=1) where {T}

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
    c_λ = @SVector [0.0, 1.0]
    d_λ = @SVector [0.0, 1.0]


    if length(d) == 0
        return TableauVPARK(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ)
    else
        @assert length(d) == q.s == p.s

        return TableauVPARK(name, o,
                            a_q, a_p, α_q, α_p,
                            a_q̃, a_p̃, α_q̃, α_p̃,
                            b_q, b_p, β_q, β_p,
                            c_q, c_p, c_λ, d_λ,
                            d)
    end

end


"Tableau for Gauss-Lobatto IIIA-IIIB method with two stages and symplectic projection."
function getTableauLobIIIAIIIB2pSymplectic()
    getTableauSymplecticProjection(:LobIIIAIIIB2pSymplectic, getCoefficientsLobIIIA2(), getCoefficientsLobIIIB2(), get_lobatto_d_vector(2); R∞=-1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with three stages and symplectic projection."
function getTableauLobIIIAIIIB3pSymplectic()
    getTableauSymplecticProjection(:LobIIIAIIIB3pSymplectic, getCoefficientsLobIIIA3(), getCoefficientsLobIIIB3(), get_lobatto_d_vector(3); R∞=+1)
end

"Tableau for Gauss-Lobatto IIIA-IIIB method with four stages and symplectic projection."
function getTableauLobIIIAIIIB4pSymplectic()
    getTableauSymplecticProjection(:LobIIIAIIIB4pSymplectic, getCoefficientsLobIIIA4(), getCoefficientsLobIIIB4(), get_lobatto_d_vector(4); R∞=-1)
end

"Tableau for Gauss-Legendre method with s stages and symplectic projection."
function getTableauGLRKpSymplectic(s)
    glrk = getCoefficientsGLRK(s)
    getTableauSymplecticProjection(Symbol("vpglrk", s, "pSymplectic"), glrk, glrk; R∞=(-1)^s)
end
