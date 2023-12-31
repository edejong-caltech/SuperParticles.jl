using SuperParticles.ParticleDistributions
using JET: @test_opt

rtol = 1e-3

# Initialization
@test_opt GammaParticleDistribution(1.0, 2.0, 3.0)
@test_opt ExponentialParticleDistribution(1.0, 2.0)
# Evaluations
dist1 = GammaParticleDistribution(1.0, 2.0, 3.0)
@test_opt dist1(1.0)
dist2 = ExponentialParticleDistribution(1.0, 2.0)
@test_opt dist2(1.0)

# moments <-> params
@test_opt get_moments(dist1)
moments1 = [10.0, 50.0, 300.0]
@test_opt update_dist_from_moments!(dist1, moments1)
@test_opt get_moments(dist2)
moments2 = [10.0, 50.0]
@test_opt update_dist_from_moments!(dist2, moments2)