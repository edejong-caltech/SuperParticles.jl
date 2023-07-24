using SuperParticles.ParticleDistributions
using SuperParticles.KernelFunctions
using SuperParticles.MultiParticleSources: weighting_fn, q_integrand_inner,
    q_integrand_outer, r_integrand

rtol = 1e-4

# weighting function
dist1 = GammaParticleDistribution(10.0, 3.0, 100.0)
pdists = [dist1]
@test weighting_fn(10.0, 1, pdists) == 1.0
@test_throws AssertionError weighting_fn(10.0, 2, pdists)

dist2 = GammaParticleDistribution(1.0, 5.0, 10.0)
pdists = [dist1, dist2]
@test weighting_fn(10.0, 1, pdists) == 0.02866906313141952
@test weighting_fn(10.0, 2, pdists) == 1.0
@test abs(weighting_fn(rtol, 1, pdists) - 1.0) <= rtol

# q_integrands
kernel = LinearKernelFunction(1.0)
dist3 = GammaParticleDistribution(2.0, 6.0, 50.0)
pdists = [dist1, dist2, dist3]
x = 50.0
y = 20.0
for j in 1:3
    for k in 1:3
        if j == k
            @test_throws AssertionError q_integrand_inner(x, y, j, k, kernel, pdists)
        else
            @test q_integrand_inner(x, y, j, k, kernel, pdists) > 0.0
            @test_throws AssertionError q_integrand_inner(y, x, j, k, kernel, pdists)
            for moment_order in 0:2
                @test q_integrand_outer(x, j, k, kernel, pdists, FT(moment_order)) > 0.0
                @test q_integrand_outer(y, j, k, kernel, pdists, FT(moment_order)) > 0.0
            end
        end
    end
end

