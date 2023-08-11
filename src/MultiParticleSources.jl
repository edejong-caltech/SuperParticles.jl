"""
  Microphysical source functions involving interactions between particles
  
  Includes terms such as collisional coalescence
"""
module MultiParticleSources

using QuadGK
using HCubature
using StaticArrays
# using SpecialFunctions: gamma

export get_coalescence_integral_moment_qrs

FT = Float64

"""
get_coalescence_integral_moment_qrs(x::Array{FT}, kernel::KernelFunction{FT}, pdist::ParticleDistribution{FT}, n_samples::Int)

Returns the collision-coalescence integral at points `x`.
Q: source term to particle k via collisions with smaller particles j
"""
# function get_coalescence_integral_moment_qrs!(
#   moment_order, kernel, pdists, Q, R, S)

#   Ndist = length(pdists)
  
#     for j in 1:Ndist
#     max_mass = ParticleDistributions.max_mass(pdists[j])
#     s1 = x -> s_integrand1(x, j, kernel, pdists, moment_order)
#     s2 = x -> s_integrand2(x, j, kernel, pdists, moment_order)
#     S[j,1] = quadgk(s1, 0.0, max_mass; rtol=1e-8)[1]
#     S[j,2] = quadgk(s2, 0.0, max_mass; rtol=1e-8)[1]
#         for k in 1:Ndist
#             max_mass = max(ParticleDistributions.max_mass(pdists[j]), ParticleDistributions.max_mass(pdists[k]))
#             R[j,k] = hcubature(xy -> r_integrand(xy[1], xy[2], j, k, kernel, pdists, moment_order), [0.0, 0.0], [max_mass, max_mass]; rtol=1e-8, maxevals=1000)[1]
#             if j < k
#                 Q[j,k] = quadgk(x -> q_integrand_outer(x, j, k, kernel, pdists, moment_order), 0.0, max_mass; rtol=1e-8)[1]
#             else
#                 Q[j,k] = 0.0
#             end
#         end
#     end
# end

function update_R_coalescence_matrix!(
    moment_order, kernel, pdists, R
)
    Ndist = length(pdists)
    for j in 1:Ndist
        for k in 1:Ndist
            R[j,k] = quadgk(x -> r_integrand_outer(x, j, k, kernel, pdists, moment_order), 0.0, Inf)[1]
            #max_mass = max(ParticleDistributions.max_mass(pdists[j]), ParticleDistributions.max_mass(pdists[k]))
            #R[j,k] = hcubature(xy -> SA[r_integrand(xy[1], xy[2], j, k, kernel, pdists, moment_order)], (0.0, 0.0), (Inf, Inf); rtol=1e-8, maxevals=1000)[1]
        end
    end
end

function weighting_fn(x, k, pdists)
    denom = 0.0
    num = 0.0
    if k > length(pdists)
        throw(AssertionError("k out of range"))
    end
    for j=1:length(pdists)
      denom += pdists[j](x) / pdists[j].n
      if j<= k
        num += pdists[j](x) / pdists[j].n
      end
    end
    if denom == 0.0
      return 0.0
    else
      return num / denom
    end
end

function q_integrand_inner(x, y, j, k, kernel, pdists)
    if j==k
       throw(AssertionError("q_integrand called on j==k, should call s instead"))
    elseif y > x
        throw(AssertionError("x <= y required in Q integrals"))
    end
    integrand = 0.5 * kernel(x - y, y) * (pdists[j](x-y) * pdists[k](y) + pdists[k](x-y) * pdists[j](y))
    return integrand
end

function q_integrand_outer(x, j, k, kernel, pdists, moment_order)
    if j==k
        throw(AssertionError("q_integrand called on j==k, should call s instead"))
    end
    outer = x.^moment_order * quadgk(yy -> q_integrand_inner(x, yy, j, k, kernel, pdists), 0.0, x)[1]
    return outer
end

function r_integrand_inner(x, y, j, k, kernel, pdists)
    integrand = kernel(x, y) * pdists[j](x) * pdists[k](y)
    return integrand
end

function r_integrand_outer(x, j, k, kernel, pdists, moment_order)
    outer = x.^moment_order * quadgk(yy -> r_integrand_inner(x, yy, j, k, kernel, pdists), 0.0, Inf)[1]
    return outer
end

function s_integrand_inner(x, k, kernel, pdists, moment_order)
    integrand_inner = y -> 0.5 * kernel(x - y, y) * pdists[k](x-y) * pdists[k](y)
    integrand_outer = x.^moment_order * quadgk(yy -> integrand_inner(yy), 0.0, x)[1]
    return integrand_outer
  end
  
function s_integrand1(x, k, kernel, pdists, moment_order)
    integrandj = weighting_fn(x, k, pdists) * s_integrand_inner(x, k, kernel, pdists, moment_order)
    return integrandj
end
  
function s_integrand2(x, k, kernel, pdists, moment_order)
    integrandk = (1 - weighting_fn(x, k, pdists)) * s_integrand_inner(x, k, kernel, pdists, moment_order)
    return integrandk
end

end # module
