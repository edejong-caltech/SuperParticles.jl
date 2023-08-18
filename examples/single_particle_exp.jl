"Test case with a single exponential distribution"

using DifferentialEquations
using Plots

using SuperParticles.KernelFunctions
using SuperParticles.ParticleDistributions
using SuperParticles.MultiParticleSources
using SuperParticles.MultiParticleSources: update_S_coalescence_matrix!,
    update_Q_coalescence_matrix!, update_R_coalescence_matrix!,
    q_integrand_inner, q_integrand_outer
using JET: @test_opt

FT = Float64
tol = 1e-4

function rhs!(ddist_moments, dist_moments, p, t)
    # update the ParticleDistributions
    for i=1:p.Ndist
        update_dist_from_moments!(p.pdists[i], dist_moments[i,:])
    end
    # update the information
    update_coal_ints!(p.Nmom, p.kernel_func, p.pdists, p.coal_data)
    ddist_moments .= p.coal_data.coal_ints
end

function plot_moments!(sol, p; plt_title="test_moments")
    Plots.gr()
    time = sol.t
    moment_plots = []
    for i in 1:p.Nmom
        ploti = plot(time,
                vcat(sol.u...)[:,i],
                linewidth=2,
                xaxis="time [s]",
                yaxis="M"*string(i-1),
                label="M_{"*string(i-1)*",1}",
                xlims=(0, maximum(time)),
                ylims=(0, 1.2*maximum(vcat(sol.u...)[:,i]))
            )
        push!(moment_plots, 
            ploti
        )
    end
    plot(moment_plots...)
    savefig("examples/"*plt_title*".png")
end

function plot_spectra!(sol, p; plt_title="test_spectra", logxrange=(0, 8))
    Plots.gr()
    x = 10 .^ (collect(range(logxrange[1], logxrange[2], 100)))
    r = (x * 3 / 4 / π) .^ (1/3)

    # initial distribution
    for i=1:p.Ndist
        update_dist_from_moments!(p.pdists[i], vcat(sol.u...)[1,:])
    end
    pinit = plot(r,
        3*x.*p.pdists[1].(x),
        linewidth=2,
        xaxis=:log,
        yaxis="dV / d(ln r)",
        xlabel="r",
        xlim=(minimum(r), maximum(r)),
        label="Initial, pdist 1"
    )

    # final distribution
    for i=1:p.Ndist
        update_dist_from_moments!(p.pdists[i], vcat(sol.u...)[end,:])
    end
    pfin = plot(r,
        3*x.*p.pdists[1].(x),
        linewidth=2,
        xaxis=:log,
        yaxis="dV / d(ln r)",
        xlabel="r",
        xlim=(minimum(r), maximum(r)),
        label="Final, pdist 1"
    )

    plot(pinit, pfin)
    savefig("examples/"*plt_title*".png")
end

function main()
    T_end = 1.0
    coalescence_coeff = 1e-3
    kernel = LinearKernelFunction(coalescence_coeff)

    # Initial condition 
    Ndist = 1
    particle_number = [100.0]
    mass_scale = [100.0]
    Nmom = 2

    # Initialize distributions
    pdists = map(1:Ndist) do i
        ExponentialParticleDistribution(particle_number[i], mass_scale[i])
    end

    dist_moments = zeros(FT, Ndist, Nmom)
    for i in 1:Ndist
        dist_moments[i,:] = get_moments(pdists[i])
    end

    # Set up ODE information
    coal_data = initialize_coalescence_data(Ndist, Nmom)
    p = (Ndist=Ndist, Nmom=Nmom, pdists=pdists, kernel_func=kernel, coal_data=coal_data)

    tspan = (0.0, T_end)
    prob = ODEProblem(rhs!, dist_moments, tspan, p; progress=true)
    sol = solve(prob, Tsit5(), reltol=tol, abstol=tol)
    @show sol.u
    plot_moments!(sol, p; plt_title="single_particle_exp_moments")
    plot_spectra!(sol, p; plt_title="single_particle_exp_spectra")
end

@time main()