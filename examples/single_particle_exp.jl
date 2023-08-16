"Test case with a single exponential distribution"

#using DifferentialEquations
#using Plots

using SuperParticles.KernelFunctions
using SuperParticles.ParticleDistributions
using SuperParticles.MultiParticleSources
using SuperParticles.MultiParticleSources: update_S_coalescence_matrix!,
    update_Q_coalescence_matrix!, update_R_coalescence_matrix!,
    q_integrand_inner, q_integrand_outer
using JET: @test_opt
FT = Float64

# function rhs!(ddist_moments, dist_moments, p, t)
#     # update the ParticleDistributions
#     for i=1:p.Ndist
#         update_dist_from_moments!(p.pdists[i], dist_moments[i])
#     end
#     #get_coalescence_integral_moment_qrs!(1.0, kernel_func, pdists, coal_data.Q, coal_data.R, coal_data.S)
#     update_coal_ints!(p.Nmom, p.kernel_func, p.pdists, p.coal_data)
#     #ddist_moments = p.coal_data.coal_ints
# end

function main()
    # T_end = 1.0
    # coalescence_coeff = 1.0
    kernel = LinearKernelFunction(1.0)

    # Initial condition 
    Ndist = 1
    particle_number = [10.0]
    mass_scale = [100.0]
    Nmom = 2

    # Initialize distributions
    pdists = map(1:Ndist) do i
        ExponentialParticleDistribution(particle_number[i], mass_scale[i])
    end
    @show pdists

    dist_moments = map(1:Ndist) do i 
        get_moments(pdists[i])
    end
    @show dist_moments

    cd = initialize_coalescence_data(length(pdists), 2)
    moment_order = 0.0
    x = 1.0
    y = 20.0
    j = 1
    k = 1
    #@test_opt q_integrand_outer(x, j, k, kernel, pdists, 0.0)
    #@test_opt update_Q_coalescence_matrix!(moment_order, kernel, pdists, cd.Q)
    #update_S_coalescence_matrix!(moment_order, kernel, pdists, cd.S)
    # get_coalescence_integral_moment_qrs!(moment_order, kernel, pdists, cd.Q, cd.R, cd.S)
    update_coal_ints!(2, kernel, pdists, cd)
    # # Set up ODE information
    # ddist_moments = zeros(FT, (Ndist, Nmom))
    # coal_data = initialize_coalescence_data(Ndist, Nmom)
    # p = (Ndist=Ndist, Nmom=Nmom, pdists=pdists, kernel_func=kernel_func, coal_data=coal_data)
    # @show ddist_moments

    # # rhs!(ddist_moments, dist_moments, p, 0.0)
    # update_coal_ints!(p.Nmom, p.kernel_func, p.pdists, p.coal_data)
    # @show p.coal_data.coal_ints

end

#@time main()