"""
The `real_roots` field contains the real roots. The complex roots are stored in
`complex_roots`, with only one half of each complex conjugate root pair stored.
To be precise, of the two roots in the pair, only the root with the positive
imaginary part is stored.
"""
struct ComplexConjugateRootTheoremRoots{F<:AbstractFloat}
    real_roots::Vector{F}
    complex_roots::Vector{Complex{F}}
end

Base.:(==)(l::ComplexConjugateRootTheoremRoots, r::ComplexConjugateRootTheoremRoots) =
    (l.real_roots == r.real_roots) && (l.complex_roots == r.complex_roots)

Base.isequal(l::ComplexConjugateRootTheoremRoots, r::ComplexConjugateRootTheoremRoots) =
    isequal(l.real_roots == r.real_roots) && isequal(l.complex_roots == r.complex_roots)

Base.hash(r::ComplexConjugateRootTheoremRoots, h::UInt) = hash(r.real_roots, hash(r.complex_roots, h))

"""
Indicates an internal error.
"""
struct ComplexConjugateException <: Exception
    count_real::Int
    count_cc_negimag::Int
    count_cc_posimag::Int
end

"""
Roots of a real polynomial with real coefficients.
"""
function real_polynomial_roots(coefs::AbstractVector{F}) where {F<:AbstractFloat}
    rts = roots(coefs)::AbstractVector{Complex{F}}
    real_roots = F[]
    complex_roots = Complex{F}[]

    # Count of complex roots with negative imaginary part, used as a sanity check.
    cn = 0

    for r ∈ rts
        if isreal(r)
            push!(real_roots, real(r))
        elseif false < imag(r)
            push!(complex_roots, r)
        else
            cn += true
        end
    end

    cp = length(complex_roots)
    (cn == cp) || throw(ComplexConjugateException(length(real_roots), cn, cp))

    ComplexConjugateRootTheoremRoots(real_roots, complex_roots)
end
