using Test

FT = Float64

@testset "SuperParticles" begin
    @testset "Correctness" begin
	    @testset "ParticleDistributions" begin
	        include("test_ParticleDistributions_correctness.jl")
          include("test_KernelFunctions_correctness.jl")
          include("test_MultiParticleSources_correctness.jl")
        end
    end

    @testset "Type stability" begin
        @testset "ParticleDistributions" begin
          include("test_ParticleDistributions_opt.jl")
          include("test_KernelFunctions_opt.jl")
          include("test_MultiParticleSources_opt.jl")
    end
  end
end
