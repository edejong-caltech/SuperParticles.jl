"Test case with a single exponential distribution"

using DifferentialEquations
#using Plots

using SuperParticles.KernelFunctions
using SuperParticles.ParticleDistributions
using SuperParticles.MultiParticleSources
using SuperParticles.MultiParticleSources: update_S_coalescence_matrix!,
    update_Q_coalescence_matrix!, update_R_coalescence_matrix!,
    q_integrand_inner, q_integrand_outer
using JET: @test_opt

FT = Float64
tol = 1e-4

function rhs!(ddist_moments, dist_moments, p, t)
    # update the ParticleDistributions
    for i=1:p.Ndist
        update_dist_from_moments!(p.pdists[i], dist_moments[i,:])
    end
    # update the information
    update_coal_ints!(p.Nmom, p.kernel_func, p.pdists, p.coal_data)
    @show p.coal_data.coal_ints
    ddist_moments .= p.coal_data.coal_ints
    @show ddist_moments
end

function main()
    T_end = 1.0
    coalescence_coeff = 1e-4
    kernel = LinearKernelFunction(coalescence_coeff)

    # Initial condition 
    Ndist = 1
    particle_number = [10.0]
    mass_scale = [100.0]
    Nmom = 2

    # Initialize distributions
    pdists = map(1:Ndist) do i
        ExponentialParticleDistribution(particle_number[i], mass_scale[i])
    end

    dist_moments = zeros(FT, Ndist, Nmom)
    for i in 1:Ndist
        dist_moments[i,:] = get_moments(pdists[i])
    end

    # Set up ODE information
    ddist_moments = zeros(FT, Ndist, Nmom)
    coal_data = initialize_coalescence_data(Ndist, Nmom)
    p = (Ndist=Ndist, Nmom=Nmom, pdists=pdists, kernel_func=kernel, coal_data=coal_data)

    tspan = (0.0, T_end)
    prob = ODEProblem(rhs!, dist_moments, tspan, p; progress=true)
    sol = solve(prob, Tsit5(), reltol=tol, abstol=tol)

end

@time main()