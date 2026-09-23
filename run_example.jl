# Run: julia --project=. run_example.jl /path/to/initial_N6.h5 [output_directory]
# The file must contain t (initial time) and rho (an MPS).

# FIRST, CREATE AN INFINITE-TEMPERATURE INPUT (run this block in a Julia REPL):
# From this folder start: julia --project=.
# Copy the following lines without their leading "# " characters; edit L and path.
# include("simulation_GCF.jl")
# L = 6
# sites = siteinds("mpdo=1/2", L; conserve_qns=QNcharge)
# rho = ident(sites)
# for j in eachindex(sites)
#     rho[j] /= 2                  # local I/2, so global rho = I/2^L and Tr(rho)=1
# end
# h5open("initial_N6.h5", "w") do f
#     write(f, "t", 0.0)          # physical starting time
#     write(f, "rho", rho)       # ITensor writes the MPS as an HDF5 group
#     write(f, "log_scale", 0.0) # no normalization has yet been removed
# end
# Then exit the REPL and run: julia --project=. run_example.jl initial_N6.h5
# The creation block is commented so a simulation never overwrites your input.

include(joinpath(@__DIR__, "simulation_GCF.jl"))
length(ARGS) in (1,2) || error("Supply initial.h5 and optionally an output directory")
initial_file = ARGS[1]
# Omit the output argument to save data.h5 and rho.h5 directly in this folder.
# Supply a different directory when you want to keep previous outputs.
directory = length(ARGS) == 2 ? ARGS[2] : @__DIR__
chi_low = 16          # bond-dimension limit for the first Steps_low_chi steps
chi = 32          # bond-dimension limit afterwards
cf = 1e-10           # truncation cutoff
eps, eta, Da = -1.0, 0.2, 0.5  # exchange, dephasing rate, XXZ anisotropy
mu, gamma, count_u = 0.3, 1.0, -0.5  # polarization, boundary rate, counting field u
delta = 0.01         # time step
Steps = 100           # additional steps; final time = saved time + Steps*delta
Steps_low_chi = 50
save_step = 5         # measure every five steps, plus initial and final states
scheme = "TEBD2"

observables = ["current", "zz_correlation"]  # remove unwanted entries

# Julia compiles functions when they are first used. This small run compiles
# evolution, selected measurements and saving before the actual simulation.
# It starts from its own state, so it does not equilibrate or modify your input.
# Its memory cost does not predict the cost of your full system.
# Use chi_low for five steps, then chi for five steps. Temporary files are removed.
mktempdir() do warm_directory
    sites = siteinds("mpdo=1/2", 5; conserve_qns=QNcharge)
    rho = ident(sites)
    for j in eachindex(sites)
        rho[j] /= 2  # identity divided by 2 at each site gives total trace one
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
@time Time_Evolution(initial_file; directory, chi, chi_low, cf, eps, eta, Da,
               mu, gamma, count_u, delta, Steps, Steps_low_chi, save_step, scheme, observables)
