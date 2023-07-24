using SuperParticles.ParticleDistributions
using SuperParticles.KernelFunctions
using SuperParticles.MultiParticleSources: weighting_fn, q_integrand_inner,
    q_integrand_outer
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
@test_opt q_integrand_outer(x, j, k, kernel, pdists, 1.0)

