# Row-vectorized local density matrix: (rho_up,up, rho_up,dn, rho_dn,up, rho_dn,dn).
# A rho B is represented by A ⊗ transpose(B). K operators act on the ket, B operators on the bra.
const I2 = Matrix{Float64}(I, 2, 2)
const sigma_x = [0.0 1.0; 1.0 0.0]
const sigma_y = [0.0 -im; im 0.0]
const sigma_z = [1.0 0.0; 0.0 -1.0]
# sigma_plus = (sigma_x + i sigma_y)/2; these are Pauli ladder operators.
const sigma_plus = [0.0 1.0; 0.0 0.0]
const sigma_minus = transpose(sigma_plus)

# Keep the original site name and charge convention for snapshot compatibility.
function ITensors.space(::SiteType"mpdo=1/2"; conserve_qns=false)
    conserve_qns || return 4
    return [QN("Qz", 0)=>1, QN("Qz", 1)=>1, QN("Qz", -1)=>1, QN("Qz", 0)=>1]
end

# Ket: K_a = sigma_a ⊗ I. Bra: B_a = I ⊗ sigma_a.
# With row vectorization B_a represents RIGHT multiplication by transpose(sigma_a):
# B_p means rho*sigma_minus, B_m means rho*sigma_plus, B_y means -rho*sigma_y.
# Keeping this convention preserves the tilted generator used by the script.
const LOCAL_OPS = Dict(
    "Iden" => kron(I2, I2),
    "K_x" => kron(sigma_x, I2), "K_y" => kron(sigma_y, I2),
    "K_z" => kron(sigma_z, I2), "K_p" => kron(sigma_plus, I2),
    "K_m" => kron(sigma_minus, I2),
    "B_x" => kron(I2, sigma_x), "B_y" => kron(I2, sigma_y),
    "B_z" => kron(I2, sigma_z), "B_p" => kron(I2, sigma_plus),
    "B_m" => kron(I2, sigma_minus),
)
# Register these names so op("K_z", site), etc. build the corresponding ITensor.
# x/y combine opposite charges: use dense sites (conserve_qns=false) for them.
# The simulation instead uses p/m, which have definite charges and support QNs.
for name in keys(LOCAL_OPS)
    @eval function ITensors.op!(tensor::ITensor, ::OpName{$(QuoteNode(Symbol(name)))},
                                ::SiteType"mpdo=1/2", site::Index)
        if hasqns(site) && $name in ("K_x", "K_y", "B_x", "B_y")
            error("x/y have mixed charge; use dense sites or the p/m ladder operators")
        end
        matrix = LOCAL_OPS[$name]
        for a in 1:4, b in 1:4
            iszero(matrix[a,b]) || (tensor[site'=>a, site=>b] = matrix[a,b])
        end
    end
end

# Unnormalized identity density matrix: trace = 2^L. Divide each local tensor
# by 2 to obtain the infinite-temperature state I/2^L with trace one.
function ident(sites)
    psi = MPS(sites)
    for j in eachindex(sites)
        psi[j] = ITensor([1.0, 0.0, 0.0, 1.0], sites[j])
    end
    return psi
end

# Apply all insertions before taking the trace; repeated sites are allowed.
# These are trace-normalized observables, including for a tilted state.
function measure_mpdo_long(psi, names, positions)
    phi = deepcopy(psi)
    for (name, j) in zip(names, positions)
        phi[j] = noprime(op(name, siteind(psi, j)) * phi[j])
    end
    identity = ident(siteinds(psi))
    return inner(identity, phi) / inner(identity, psi)
end
measure_mpdo_single(psi, name, j) = measure_mpdo_long(psi, [name], [j])
