# Alternative Hamiltonians and Lindblad dissipators for simulation_GCF.jl.
#
# Simple use:
#   include("simulation_GCF.jl")
#   include("Hs_and_Ls.jl")
#   gates = Time_propagator(H_oper_XXZ, Lind_Dissipator_dephasing, sites,
#                           delta, scheme;
#                           eps, Da, eta, mu, gamma, count_u)
#
# You can also copy one function below over H_oper or Lind_Dissipator in
# simulation_GCF.jl and remove its model suffix. Every public function keeps the
# keyword interface expected by Time_propagator:
#
#   H_oper_MODEL(sites; eps, Da, count_u)
#   Lind_Dissipator_MODEL(sites; eta, mu, gamma)
#
# Scope: nearest-neighbour spin-1/2 Hamiltonians and on-site jump operators.
# Long-range Hamiltonians and multi-site jumps need changes to local_generators.
# All formulas use Pauli matrices, not spin operators S=sigma/2.
# The choices are a compact set of canonical arXiv condensed-matter models,
# not a quantitative bibliometric ranking of every model on arXiv.
#
# the choices below:
#   XXZ and boundary driving
#   Transverse-field Ising chain
#   Anisotropic XY chain
#   Boundary-driven XX dynamics
#   Bulk dephasing and boundary drive
#
# The optional validation scripts are archived separately from this minimal package.

using ITensors, LinearAlgebra

# This file may be included after simulation_GCF.jl or directly by itself.
if !isdefined(@__MODULE__, :LOCAL_OPS)
    include(joinpath(@__DIR__, "mpdo_basics.jl"))
end

# =============================================================================
# Hamiltonians
# =============================================================================

"""
    H_oper_XXZ(sites; eps, Da, count_u)

XXZ Hamiltonian used by `simulation_GCF.jl`:

    H = eps * sum(sigma_x*sigma_x + sigma_y*sigma_y + Da*sigma_z*sigma_z).

The total counting field `count_u` is distributed uniformly over the `L-1`
bonds. This function is the suffixed, copy-ready version of `H_oper`.
"""
function H_oper_XXZ(sites; eps, Da, count_u)
    L = length(sites)
    L >= 2 || error("Need at least two sites")

    tilt = count_u / (2 * (L - 1))
    terms = []
    for j in 1:L-1
        a, b = sites[j], sites[j+1]
        ket = eps * (
            2exp(-tilt) * op("K_p", a)*op("K_m", b) +
            2exp(tilt) * op("K_m", a)*op("K_p", b) +
            Da * op("K_z", a)*op("K_z", b))
        bra = eps * (
            2exp(-tilt) * op("B_p", a)*op("B_m", b) +
            2exp(tilt) * op("B_m", a)*op("B_p", b) +
            Da * op("B_z", a)*op("B_z", b))
        push!(terms, (ket=ket, bra=bra))
    end
    return terms
end



"""
    H_oper_TFIM(sites; eps, Da, count_u)

Transverse-field Ising model in Pauli convention:

    H = -eps * sum(sigma_z[j]*sigma_z[j+1]) - Da * sum(sigma_x[j]).

Parameter map: `eps = J` and `Da = h`. The on-site field is shared between
neighbouring bond generators; edge sites receive the full field exactly once.

Set `count_u=0`: no current-counting deformation is defined here. This model
breaks z-magnetization conservation, so create states with `QNcharge=false`.

"""
function H_oper_TFIM(sites; eps, Da, count_u)
    L = length(sites)
    L >= 2 || error("Need at least two sites")
    iszero(count_u) || error("H_oper_TFIM has no counting-field tilt; use count_u=0")
    hasqns(first(sites)) && error("TFIM requires QNcharge=false because sigma_x breaks U(1)")

    terms = []
    for j in 1:L-1
        a, b = sites[j], sites[j+1]
        left_weight = j == 1 ? 1.0 : 0.5
        right_weight = j == L-1 ? 1.0 : 0.5

        ket = -eps * op("K_z", a)*op("K_z", b)
        ket -= Da * (left_weight * op("K_x", a)*op("Iden", b) +
                     right_weight * op("Iden", a)*op("K_x", b))
        bra = -eps * op("B_z", a)*op("B_z", b)
        bra -= Da * (left_weight * op("B_x", a)*op("Iden", b) +
                     right_weight * op("Iden", a)*op("B_x", b))
        push!(terms, (ket=ket, bra=bra))
    end
    return terms
end


"""
    H_oper_XY(sites; eps, Da, count_u)

Anisotropic XY chain without a transverse field:

    H = eps * sum((1+Da)*sigma_x*sigma_x + (1-Da)*sigma_y*sigma_y).

Parameter map: `eps = J` and `Da = gamma_XY`. At `Da=0` this is the XX chain;
at `Da=1` it is `2eps*sum(sigma_x*sigma_x)`.

Set `count_u=0`: no tilt is defined for the pair-creation terms. For `Da != 0`
the Hamiltonian conserves parity rather than z magnetization, so use
`QNcharge=false`.

"""
function H_oper_XY(sites; eps, Da, count_u)
    L = length(sites)
    L >= 2 || error("Need at least two sites")
    iszero(count_u) || error("H_oper_XY has no counting-field tilt; use count_u=0")
    hasqns(first(sites)) && error("Anisotropic XY requires QNcharge=false")

    terms = []
    for j in 1:L-1
        a, b = sites[j], sites[j+1]
        # (1+g)XX+(1-g)YY = 2(PM+MP)+2g(PP+MM).
        ket = 2eps * (
            op("K_p", a)*op("K_m", b) + op("K_m", a)*op("K_p", b) +
            Da * (op("K_p", a)*op("K_p", b) + op("K_m", a)*op("K_m", b)))
        bra = 2eps * (
            op("B_p", a)*op("B_m", b) + op("B_m", a)*op("B_p", b) +
            Da * (op("B_p", a)*op("B_p", b) + op("B_m", a)*op("B_m", b)))
        push!(terms, (ket=ket, bra=bra))
    end
    return terms
end


# =============================================================================
# Lindblad dissipators
# =============================================================================

"""
    Lind_Dissipator_none(sites; eta, mu, gamma)

No dissipation: returns a zero local generator at every site. The three keywords
are accepted only to preserve the simulation interface.

"""
function Lind_Dissipator_none(sites; eta, mu, gamma)
    terms = ITensor[]
    for s in sites
        # A zero 4x4 matrix is the one-site generator D = 0.
        D = zeros(ComplexF64, 4, 4)
        push!(terms, ITensor(D, s', dag(s)))
    end
    return terms
end


"""
    Lind_Dissipator_dephasing(sites; eta, mu, gamma)

Independent bulk dephasing on every site with jump `sqrt(eta)*sigma_z`.
`mu` and `gamma` are unused. For one site the Liouvillian spectrum is
`[0,0,-2eta,-2eta]`: populations are fixed and coherences decay.


"""
function Lind_Dissipator_dephasing(sites; eta, mu, gamma)
    eta >= 0 || error("eta must be nonnegative")

    terms = ITensor[]
    for s in sites
        jumps = [sqrt(eta) * sigma_z]
        D = zeros(ComplexF64, 4, 4)
        for Lb in jumps
            Q = Lb' * Lb
            D += kron(Lb, conj(Lb)) -
                 (kron(Q, I2) + kron(I2, transpose(Q))) / 2
        end
        push!(terms, ITensor(D, s', dag(s)))
    end
    return terms
end


"""
    Lind_Dissipator_loss(sites; eta, mu, gamma)

Local amplitude damping on every site with jump `sqrt(eta)*sigma_minus`.
`mu` and `gamma` are unused. The unique one-site steady state is spin down.

"""
function Lind_Dissipator_loss(sites; eta, mu, gamma)
    eta >= 0 || error("eta must be nonnegative")

    terms = ITensor[]
    for s in sites
        jumps = [sqrt(eta) * sigma_minus]
        D = zeros(ComplexF64, 4, 4)
        for Lb in jumps
            Q = Lb' * Lb
            D += kron(Lb, conj(Lb)) -
                 (kron(Q, I2) + kron(I2, transpose(Q))) / 2
        end
        push!(terms, ITensor(D, s', dag(s)))
    end
    return terms
end


"""
    Lind_Dissipator_polarizing(sites; eta, mu, gamma)

Local gain and loss on every site:

    L_plus  = sqrt(eta*(1+mu))*sigma_plus
    L_minus = sqrt(eta*(1-mu))*sigma_minus.

The one-site stationary polarization is `<sigma_z>=mu`. `gamma` is unused.
Require `eta>=0` and `abs(mu)<=1`.

"""
function Lind_Dissipator_polarizing(sites; eta, mu, gamma)
    eta >= 0 || error("eta must be nonnegative")
    abs(mu) <= 1 || error("Require abs(mu) <= 1")

    terms = ITensor[]
    for s in sites
        jumps = [sqrt(eta*(1+mu)) * sigma_plus,
                 sqrt(eta*(1-mu)) * sigma_minus]
        D = zeros(ComplexF64, 4, 4)
        for Lb in jumps
            Q = Lb' * Lb
            D += kron(Lb, conj(Lb)) -
                 (kron(Q, I2) + kron(I2, transpose(Q))) / 2
        end
        push!(terms, ITensor(D, s', dag(s)))
    end
    return terms
end


"""
    Lind_Dissipator_boundary(sites; eta, mu, gamma)

Oppositely polarized reservoirs acting only at the first and last sites. This is
the boundary part of the current simulation; `eta` is unused. The left jumps are
`sqrt(gamma*(1-mu))*sigma_plus` and `sqrt(gamma*(1+mu))*sigma_minus`; the right
rates are reversed.

"""
function Lind_Dissipator_boundary(sites; eta, mu, gamma)
    gamma >= 0 || error("gamma must be nonnegative")
    abs(mu) <= 1 || error("Require abs(mu) <= 1")

    terms = ITensor[]
    for (j, s) in enumerate(sites)
        # Bulk sites have no jumps in this boundary-only model.
        jumps = Matrix{ComplexF64}[]
        if j == 1
            push!(jumps, sqrt(gamma*(1-mu)) * sigma_plus,
                         sqrt(gamma*(1+mu)) * sigma_minus)
        elseif j == length(sites)
            push!(jumps, sqrt(gamma*(1+mu)) * sigma_plus,
                         sqrt(gamma*(1-mu)) * sigma_minus)
        end

        D = zeros(ComplexF64, 4, 4)
        for Lb in jumps
            Q = Lb' * Lb
            D += kron(Lb, conj(Lb)) -
                 (kron(Q, I2) + kron(I2, transpose(Q))) / 2
        end
        push!(terms, ITensor(D, s', dag(s)))
    end
    return terms
end


"""
    Lind_Dissipator_dephasing_boundary(sites; eta, mu, gamma)

Bulk dephasing plus oppositely polarized boundary reservoirs. This is identical
to `Lind_Dissipator` in simulation_GCF.jl and is the standard choice for the
current tilted-XXZ example.

"""
function Lind_Dissipator_dephasing_boundary(sites; eta, mu, gamma)
    eta >= 0 || error("eta must be nonnegative")
    gamma >= 0 || error("gamma must be nonnegative")
    abs(mu) <= 1 || error("Require abs(mu) <= 1")

    terms = ITensor[]
    for (j, s) in enumerate(sites)
        # Every site dephases; the two edge sites also couple to reservoirs.
        jumps = [sqrt(eta) * sigma_z]
        if j == 1
            push!(jumps, sqrt(gamma*(1-mu)) * sigma_plus,
                         sqrt(gamma*(1+mu)) * sigma_minus)
        elseif j == length(sites)
            push!(jumps, sqrt(gamma*(1+mu)) * sigma_plus,
                         sqrt(gamma*(1-mu)) * sigma_minus)
        end

        D = zeros(ComplexF64, 4, 4)
        for Lb in jumps
            Q = Lb' * Lb
            # Row-vectorized form of
            # Lb*rho*Lb' - (Lb'*Lb*rho + rho*Lb'*Lb)/2.
            D += kron(Lb, conj(Lb)) -
                 (kron(Q, I2) + kron(I2, transpose(Q))) / 2
        end
        push!(terms, ITensor(D, s', dag(s)))
    end
    return terms
end
