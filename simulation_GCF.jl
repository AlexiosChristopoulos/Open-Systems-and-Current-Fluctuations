# Run with --project=. so imports use this folder's Project.toml.
using ITensors, ITensorMPS, LinearAlgebra, HDF5
include(joinpath(@__DIR__, "mpdo_basics.jl"))

# New example states conserve the ket-minus-bra charge. Loaded MPS indices keep
# their original charge structure; gates are always built using those indices.
const QNcharge = true

"""
    H_oper(sites; eps, Da, count_u)

Return the ket and bra Hamiltonian terms, one pair per nearest-neighbour bond.
At count_u=0, H = eps Σ(σxσx + σyσy + Da σzσz)=eps Σ(2 (σpσm + σmσp) + Da σzσz) , with σz=diag(1,-1).
The counting field count_u is distributed over the bulk such that:
the ket and bra hopping terms carry exp(∓count_u/[2(L-1)]).
The coherent part of the Liouvillian is -im*(ket-bra) and the  Hamiltonian is not Hermitian at finite tilt.
"""
function H_oper(sites; eps, Da, count_u)
    L = length(sites)
    L >= 2 || error("Need at least two sites")
    
    # Share the total counting field u over the L-1 bonds. K and B denote
    # multiplication on the ket and bra factors of the vectorized density matrix.
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
    Lind_Dissipator(sites; eta, mu, gamma)

Define physical jump matrices J on each site and return their dissipators.
D[L]ρ = LρL† - (L†Lρ + ρL†L)/2;
Bulk: sqrt(eta)*σz. Left: sqrt(gamma*(1∓mu))*σ±.
Right: sqrt(gamma*(1±mu))*σ±. This is the averaged deterministic v7 generator.

"""
function Lind_Dissipator(sites; eta, mu, gamma)
    
    terms = ITensor[]
    for (j, s) in enumerate(sites)
        jumps = [sqrt(eta) * sigma_z]# this is the dephasing term and jumps is a list filled with all of the Lindblad operators in the model
        if j == 1
            push!(jumps, sqrt(gamma*(1-mu))*sigma_plus, sqrt(gamma*(1+mu))*sigma_minus) #this is for the boundary driving term at site 1
        elseif j == length(sites)
            push!(jumps, sqrt(gamma*(1+mu))*sigma_plus, sqrt(gamma*(1-mu))*sigma_minus) # this is the boundary driving term at site N
        end
        D = zeros(ComplexF64, 4, 4)
        for Lb in jumps
            Q = Lb' * Lb #this implements the Lb†Lb
            D += kron(Lb, conj(Lb)) - (kron(Q, I2) + kron(I2, transpose(Q)))/2 #This is the D[Lb] dissipator in the vectorised picture
        end
        push!(terms, ITensor(D, s', dag(s)))
    end
    return terms
end

# Assign each one-site dissipator exactly once. The last bond carries both
# D[L-1] and D[L]. Disjoint odd/even bonds then commute within each layer.
function local_generators(H_function, Lind_function, sites; eps, Da, eta, mu, gamma, count_u)
    H = H_function(sites; eps, Da, count_u)
    D = Lind_function(sites; eta, mu, gamma)
    terms = ITensor[]
    for j in 1:length(sites)-1
        term = -im*(H[j].ket - H[j].bra) + D[j]*op("Iden", sites[j+1])
        if j == length(sites)-1
            term += op("Iden", sites[j])*D[j+1]
        end
        push!(terms, term)
    end
    return terms
end



"""
    Time_propagator(H_oper, Lind_Dissipator, sites, Deltat, scheme; parameters...)

Return gates in application order. TEBD = first order; TEBD2 = Strang second
order; TEBD4 = five symmetric second-order steps (Suzuki fourth order).
"""
function Time_propagator(H_function, Lind_function, sites, Deltat, scheme;
                         eps, Da, eta, mu, gamma, count_u)
    Deltat > 0 || error("Deltat must be positive")
    scheme in ("TEBD", "TEBD2", "TEBD4") || error("Use TEBD, TEBD2 or TEBD4")
    terms = local_generators(H_function, Lind_function, sites; eps, Da, eta, mu, gamma, count_u)
    # Gates on disjoint bonds commute, so each layer can be applied in order.
    odd = 1:2:length(terms)
    even = 2:2:length(terms)
    # The explicit element type also keeps an empty even layer valid for L=2.
    layer(bonds, dt) = ITensor[exp(dt*terms[j]) for j in bonds]
    second_order(dt) = [layer(odd, dt/2); layer(even, dt); layer(odd, dt/2)]
    scheme == "TEBD" && return [layer(odd, Deltat); layer(even, Deltat)]
    scheme == "TEBD2" && return second_order(Deltat)
    p = 1 / (4 - 4^(1/3))  # Suzuki coefficient; the middle substep is negative
    
    gates = ITensor[]
    for weight in (p, p, 1-4p, p, p)  # A TEBD4 is consisted by  five TEBD2  substeps of time step p*Deltat, p*Deltat , (1-4p)*Deltat, p*Deltat,p*Deltat
        append!(gates, second_order(weight*Deltat))
    end
    return gates
end

# Choose one observable by name. All insertions are Pauli operators.
# The current and connected zz correlation are the supported observables.
function measurements(rho, observable::AbstractString; eps)
    L = length(rho)
    if observable == "current"
        # Positive current flows from i to i+1. For Pauli magnetization:
        # d<sigma_z(i)>/dt = <j_(i-1) - j_i> (plus boundary dissipators).
        # j_i = 2eps*(sigma_x sigma_y - sigma_y sigma_x)
        #     = 4im*eps*(sigma_plus sigma_minus - sigma_minus sigma_plus).
        # The eps factor is essential, including its sign. The current of Sz
        # or occupation (1+sigma_z)/2 is HALF this Pauli magnetization current.
        # At nonzero tilt this is an observable in the tilted state, not by
        # itself a derivative of the dominant tilted eigenvalue.
        return [4eps*real(im*(measure_mpdo_long(rho, ["K_p","K_m"], [i,i+1]) -
                         measure_mpdo_long(rho, ["K_m","K_p"], [i,i+1]))) for i in 1:L-1]
    elseif observable == "zz_correlation"
        z = [real(measure_mpdo_single(rho, "K_z", i)) for i in 1:L]
        return [real(measure_mpdo_long(rho, ["K_z","K_z"], [i,j])) - z[i]*z[j]
                for i in 1:L, j in 1:L]
    else
        error("Unknown observable: $observable")
    end
end

const OBSERVABLES = ["current", "zz_correlation"]

# Apply one time step and keep the MPS norm manageable.
# log_scale records the normalization removed from rho, so Log_trace is correct.
function evolution_step(rho, gates, chi, cf, log_scale)
    rho = apply(gates, rho; cutoff=cf, maxdim=chi)
    # rescale_factor is the norm removed from the MPS after this time step.
    # Dividing by it keeps the stored tensors numerically well scaled. Its
    # logarithm is accumulated in log_scale so the physical trace can later be
    # reconstructed without losing the exponential growth or decay information.
    rescale_factor = norm(rho)
    isfinite(rescale_factor) && rescale_factor > 0 ||
        error("Invalid state norm; check the time step")
    rho[1] /= rescale_factor
    return rho, log_scale + log(rescale_factor)
end

# Always start from a file containing t and rho. Steps means ADDITIONAL steps.
function Time_Evolution(initial_file; directory, chi, chi_low, cf,
                        eps, eta, Da, mu, gamma, count_u, delta, Steps,
                        Steps_low_chi, save_step, scheme, observables=OBSERVABLES)
    all(name -> name in OBSERVABLES, observables) || error("Use observables from $OBSERVABLES")
    observables = unique(observables)
    Steps >= 0 && Steps_low_chi >= 0 && save_step > 0 || error("Invalid step counts")
    chi > 0 && chi_low > 0 && cf >= 0 || error("Invalid truncation settings")
    # HDF5 stores scalars/arrays as datasets and the MPS as a structured group.
    # read(..., MPS) uses ITensor's HDF5 reader to reconstruct the tensor indices.
    initial = h5open(initial_file, "r") do f
        Dict("t" => read(f,"t"), "rho" => read(f,"rho",MPS),
             "log_scale" => (haskey(f,"log_scale") ? read(f,"log_scale") : 0.0))
    end
    t0 = initial["t"]                  # starting time, not necessarily zero
    rho = initial["rho"]              # saved vectorized density matrix (MPS)
    log_scale = get(initial, "log_scale", 0.0)  # removed log normalization
    sites = siteinds(rho)  # Gates must use the actual indices of the saved state.
    L = length(rho)
    identity = ident(sites)
    gates = Time_propagator(H_oper, Lind_Dissipator, sites, delta, scheme;
                           eps, Da, eta, mu, gamma, count_u)

    # Save both dictionaries directly in the directory chosen by the user.
    # Reusing the same directory overwrites these two files.
    mkpath(directory)
    datafile = joinpath(directory, "data.h5")
    snapshot = joinpath(directory, "rho.h5")

    # Measure the initial state, every save_step steps, and the final state.
    saved_steps = collect(0:save_step:Steps)
    last(saved_steps) != Steps && push!(saved_steps, Steps)
    times = t0 .+ delta .* saved_steps
    count = length(times)
    # Allocate only the requested measurements; time is the first array dimension.
    shapes = Dict("current" => (L-1,), "zz_correlation" => (L,L))
    measured = Dict(name => zeros(count, shapes[name]...) for name in observables)
    logs = zeros(count)
    checkpoint_every = max(1, div(Steps, 15))
    k = 0  # number of measurements completed, not the number of evolution steps

    for step in 0:Steps
        # Step zero is the loaded state. The first Steps_low_chi applications
        # use chi_low; all subsequent applications use the requested chi.
        if step > 0
            bond_dimension = step <= Steps_low_chi ? chi_low : chi
            rho, log_scale = evolution_step(rho, gates, bond_dimension, cf, log_scale)
            # Collect unused allocations more often during the larger-chi stage.
            gc_interval = step <= Steps_low_chi ? 25 : 10
            if step % gc_interval == 0
                GC.gc(false)
            end
        end
        t = t0 + step * delta
        if step % save_step == 0 || step == Steps
            k += 1
            for name in observables
                selectdim(measured[name], 1, k) .= measurements(rho, name; eps)
            end
            logs[k] = log(abs(inner(identity, rho))) + log_scale
        end
        # Measurements accumulate in memory; checkpoints rewrite the two files
        # roughly 15 times per run. Always save the final state, even for Steps=0.
        if step == Steps || (step > 0 && step % checkpoint_every == 0)
            # Save only completed measurements, plus the current state and time.
            # Two ordinary dictionaries: measurements and the current state.
            # Slicing through k excludes future rows that have not been measured.
            parameters = Dict("L"=>L, "eps"=>eps, "eta"=>eta, "Da"=>Da,
                "mu"=>mu, "gamma"=>gamma, "count_u"=>count_u, "delta"=>delta,
                "chi"=>chi, "chi_low"=>chi_low, "cf"=>cf, "Steps"=>Steps,
                "Steps_low_chi"=>Steps_low_chi, "save_step"=>save_step, "scheme"=>scheme)
            data = Dict{String,Any}(
                "t"=>times[1:k], "saved_steps"=>saved_steps[1:k],
                "Log_trace"=>logs[1:k], "observables"=>observables,
                "parameters"=>parameters, "initial_file"=>abspath(initial_file))
            for name in observables
                data[name] = Array(selectdim(measured[name], 1, 1:k))
            end
            state = Dict("t"=>t, "rho"=>rho, "log_scale"=>log_scale,
                         "parameters"=>parameters)
            # Save the dictionaries directly: ordinary entries become datasets,
            # and the nested parameter dictionary becomes the HDF5 group parameters/.
            # ITensor supplies write(f, "rho", rho) for the MPS stored in state.
            for (filename, contents) in ((datafile, data), (snapshot, state))
                h5open(filename, "w") do f
                    for (key, value) in contents
                        if key == "parameters"
                            group = create_group(f, key)
                            for (parameter, setting) in value
                                write(group, parameter, setting)
                            end
                            close(group)
                        else
                            write(f, key, value)
                        end
                    end
                end
            end
            println("Checkpoint: step=$step/$Steps, t=$t, measurements=$k, maxlinkdim=$(maxlinkdim(rho))")
            flush(stdout)  # show progress promptly in the cluster log
        end
    end
    return (; rho, times, logs, log_scale, datafile, snapshot)
end



## This is the function that can be used when you want to externally input the parameters of the simulation. Otherwise call directly the Time_Evolution(...) function as showed in run_example.jl  
function main(args=ARGS)
    if isempty(args) || args[1] in ("--help", "-h")
        println("julia --project=. simulation_GCF.jl initial.h5 directory chi chi_low cf eps eta Da mu gamma count_u delta Steps Steps_low_chi save_step current,zz_correlation [TEBD2]")
        println("The input file supplies t and rho; Steps additional steps are applied.")
        return
    end
    length(args) in (16,17) || error("Expected 16 arguments and an optional TEBD scheme; use --help")

    # Read the arguments in the same order as the help line above.
    initial_file = args[1]
    directory = args[2]
    chi = parse(Int, args[3])
    chi_low = parse(Int, args[4])
    cf = parse(Float64, args[5])
    eps = parse(Float64, args[6])
    eta = parse(Float64, args[7])
    Da = parse(Float64, args[8])
    mu = parse(Float64, args[9])
    gamma = parse(Float64, args[10])
    count_u = parse(Float64, args[11])
    delta = parse(Float64, args[12])
    Steps = parse(Int, args[13])
    Steps_low_chi = parse(Int, args[14])
    save_step = parse(Int, args[15])
    observables = strip.(split(args[16], ','))
    all(name -> name in OBSERVABLES, observables) || error("Use observables from $OBSERVABLES")
    scheme = length(args) == 17 ? args[17] : "TEBD2"
    println("Warm-up: 5 sites, 10 steps, save_step=2, Steps_low_chi=5")
    # Julia compiles methods on first use. Running this small example first
    # compiles evolution, measurements and saving before timing the main run.
    # It is a compilation warm-up, not equilibration of the user's state.
    # A five-site run cannot reproduce the memory cost of a large system.
    # Use chi_low for five steps, then chi for five steps. Temporary files are removed.
    mktempdir() do warm_directory
        sites = siteinds("mpdo=1/2", 5; conserve_qns=QNcharge)
        rho = ident(sites)
        for j in eachindex(sites)
            rho[j] /= 2  # infinite temperature, trace one
        end
        warm_file = joinpath(warm_directory, "initial.h5")
        h5open(warm_file, "w") do f
            write(f, "t", 0.0)
            write(f, "rho", rho)
            write(f, "log_scale", 0.0)
        end
        Time_Evolution(warm_file; directory=warm_directory, chi, chi_low, cf,
            eps, eta, Da, mu, gamma, count_u, delta, Steps=10,
            Steps_low_chi=5, save_step=2, scheme, observables)
    end
    println("Evolving the saved state for $Steps additional steps using $scheme")
    result = @time Time_Evolution(initial_file; directory, chi, chi_low, cf,
        eps, eta, Da, mu, gamma, count_u, delta, Steps, Steps_low_chi, save_step, scheme, observables)
    println("Saved measurements: ", result.datafile)
    println("Saved final state: ", result.snapshot)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
