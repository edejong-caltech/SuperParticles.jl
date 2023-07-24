"""
  Microphysical source functions involving interactions between particles
  
  Includes terms such as collisional coalescence
"""
module MultiParticleSources

using QuadGK
# using HCubature
# using SpecialFunctions: gamma

export get_coalescence_integral_moment_qrs

FT = Float64

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

function r_integrand(x, y, j, k, kernel, pdists, moment_order)
    integrand = x.^moment_order * kernel(x, y) * pdists[j](x) * pdists[k](y)
    return integrand
end


end # module
