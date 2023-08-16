"Test case with a single exponential distribution"

#using DifferentialEquations
#using Plots

using SuperParticles.KernelFunctions
using SuperParticles.ParticleDistributions
using SuperParticles.MultiParticleSources

FT = Float64

T_end = 1.0
coalescence_coeff = 1.0
kernel_func = ConstantKernelFunction(coalescence_coeff)

# Initial condition 
Ndist = 1
particle_number = [10.0]
k = [2.0]
mass_scale = [1.0]
Nmom = 3 #2

# Initialize distributions
pdists = map(1:Ndist) do i
    GammaParticleDistribution(particle_number[i], k[i], mass_scale[i])
    #ExponentialParticleDistribution(particle_number[i], mass_scale[i])
end
dist_moments = map(1:Ndist) do i 
    get_moments(pdists[i])
end

# Set up ODE information
ddist_moments = zeros(FT, (Ndist, Nmom))
coal_data = initialize_coalescence_data(Ndist, Nmom)
p = (Ndist=Ndist, Nmom=Nmom, pdists=pdists, kernel_func=kernel_func, coal_data=coal_data)
@show ddist_moments

function rhs!(ddist_moments, dist_moments, p, t)
    # update the ParticleDistributions
    for i=1:p.Ndist
        update_dist_from_moments!(p.pdists[i], dist_moments[i])
    end

    update_coal_ints!(p.Nmom, p.kernel_func, p.pdists, p.coal_data)
    ddist_moments = p.coal_data.coal_ints
end

#rhs!(ddist_moments, dist_moments, p, 0.0)
update_coal_ints!(p.Nmom, p.kernel_func, p.pdists, p.coal_data)
@show ddist_moments

