using SuperParticles.ParticleDistributions
using SuperParticles.KernelFunctions
using SuperParticles.MultiParticleSources: weighting_fn, q_integrand_inner,
    q_integrand_outer, r_integrand_inner, r_integrand_outer,
    s_integrand1, s_integrand2, s_integrand_inner,
    update_R_coalescence_matrix!
using JET: @test_opt

rtol = 1e-4

# weighting function
dist1 = GammaParticleDistribution(10.0, 3.0, 100.0)
dist2 = GammaParticleDistribution(1.0, 5.0, 10.0)
pdists = [dist1, dist2]
@test_opt weighting_fn(10.0, 1, pdists)
@test_opt weighting_fn(8.0, 2, pdists)

# q_integrands
x = 50.0
y = 20.0
j = 1
k = 2
kernel = LinearKernelFunction(1.0)
@test_opt q_integrand_inner(x, y, j, k, kernel, pdists)
@test_opt q_integrand_outer(x, j, k, kernel, pdists, 0.0)
@test_opt q_integrand_outer(x, j, k, kernel, pdists, 1.0)
@test_opt q_integrand_outer(x, j, k, kernel, pdists, 1.5)

# r_integrands
@test_opt r_integrand_inner(x, y, j, k, kernel, pdists)
@test_opt r_integrand_outer(x, j, k, kernel, pdists, 0.0)
@test_opt r_integrand_outer(x, j, k, kernel, pdists, 1.0)

# s_integrands
@test_opt s_integrand_inner(x, k, kernel, pdists, 1.0)
@test_opt s_integrand_inner(x, k, kernel, pdists, 1.0)
@test_opt s_integrand_inner(x, k, kernel, pdists, 1.0)

# overall Q R S fill matrices 
# n = 1
pdists = [dist1]
Q = zeros(1)
R = zeros(1)
S = zeros(1)
moment_order = 0.0
#@test_opt update_R_coalescence_matrix!(moment_order, kernel, pdists, R)
#@test_opt get_coalescence_integral_moment_qrs!(moment_order, kernel, pdists, Q, R, S)
