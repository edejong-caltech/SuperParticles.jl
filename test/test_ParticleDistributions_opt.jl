using SuperParticles.ParticleDistributions
using JET: @test_opt

rtol = 1e-3

# Initialization
@test_opt GammaParticleDistribution(1.0, 2.0, 3.0)
@test_opt ExponentialParticleDistribution(1.0, 2.0)
# Evaluations
dist = GammaParticleDistribution(1.0, 2.0, 3.0)
@test_opt dist(1.0)
dist = ExponentialParticleDistribution(1.0, 2.0)
@test_opt dist(1.0)
