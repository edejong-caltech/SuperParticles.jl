using Test

FT = Float64

@testset "SuperParticles" begin
    @testset "Correctness" begin
	    @testset "ParticleDistributions" begin
	        include("test_ParticleDistributions_correctness.jl")
        end
    end

    @testset "Type stability" begin
        @testset "ParticleDistributions" begin
          include("test_ParticleDistributions_opt.jl")
    end
  end
end
