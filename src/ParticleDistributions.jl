"""
  particle mass distribution module

Particle mass distribution functions for microphysical process modeling:
  - computation of densities of distributions
  - sampling from densities
"""
module ParticleDistributions

using Distributions: Distribution, Gamma, Exponential, pdf
using DocStringExtensions
using SpecialFunctions: gamma
using Random: rand

# particle mass distributions available for microphysics
#export ParticleDistribution
export GammaParticleDistribution
export ExponentialParticleDistribution

# methods that query particle mass distributions

"""
  ParticleDistribution{FT}

A particle mass distribution function, which can be initialized
for various subtypes of assumed shapes in the microphysics parameterization.
"""
abstract type ParticleDistribution{FT} end


"""
  GammaParticleDistribution{FT} <: ParticleDistribution{FT}

Represents particle mass distribution function of gamma shape.

# Constructors
  GammaParticleDistribution(n::Real, θ::Real, k::Real)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct GammaParticleDistribution{FT, D <: Distribution} <: ParticleDistribution{FT}
  "normalization constant (e.g., droplet number concentration)"
  n::FT
  "shape parameter"
  k::FT
  "scale parameter"
  θ::FT
  "underlying probability distribution"
  dist::D

  function GammaParticleDistribution(n::FT, k::FT, θ::FT) where {FT<:Real}
    if n < 0
        throw(DomainError(n, "n must be nonnegative"))
    end
    dist = Gamma(k, θ)
    return new{FT, typeof(dist)}(n, k, θ, dist)
  end
end

"""
  ExponentialParticleDistribution{FT} <: ParticleDistribution{FT}

Represents particle mass distribution function of exponential shape.

# Constructors
  ExponentialParticleDistribution(n::Real, θ::Real)

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ExponentialParticleDistribution{FT, D <: Distribution} <: ParticleDistribution{FT}
  "normalization constant (e.g., droplet number concentration)"
  n::FT
  "scale parameter"
  θ::FT
  "underlying probability distribution"
  dist::D

  function ExponentialParticleDistribution(n::FT, θ::FT) where {FT<:Real}
    if n < 0
      throw(DomainError(n, "n must be nonnegative"))
    end
    dist = Exponential(θ)
    return new{FT, typeof(dist)}(n, θ, dist)
  end
end

"""
  (pdist::ParticleDistribution{FT}(x::FT)

  - `x` - is an array of points to evaluate the density of `pdist` at
Returns the particle mass density evaluated at `x`.
"""
function (pdist::ParticleDistribution{FT})(x::FT) where {FT<:Real}
  return pdist.n * pdf(pdist.dist, x)
end



end