
using SuperParticles.ParticleDistributions

rtol = 1e-3

# Gamma distribution
# Initialization
dist = GammaParticleDistribution(1.0, 2.0, 3.0)
@test (dist.n, dist.θ, dist.k) == (FT(1.0), FT(3.0), FT(2.0))
@test_throws Exception GammaParticleDistribution(-1.0, 2.0, 3.0)
@test_throws Exception GammaParticleDistribution(1.0, -2.0, 3.0)
@test_throws Exception GammaParticleDistribution(1.0, 2.0, -3.0)

# Evaluate
@test dist(1.0) == 0.07961459006375433
@test max_mass(dist) == 48.42640687119285