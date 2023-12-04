using Pkg
function installed()
	deps = Pkg.dependencies()
	installs = Dict{String, VersionNumber}()
	for (uuid, dep) in deps
		dep.is_direct_dep || continue
		dep.version === nothing && continue
		installs[dep.name] = dep.version
	end
	return installs
end
# Check if packages are installed, else install them
Packages = ["Plots", "Test", "LaTeXStrings", "CSV", "DataFrames"]
installed_Packages = keys(installed())
for Package in Packages
	if !(Package in installed_Packages)
		try
			eval(Meta.parse("using $Package"))
		catch
			println("Package $Package was not found. Installation started")
			Pkg.add(Package)
			eval(Meta.parse("using $Package"))
		end
	else
		eval(Meta.parse("using $Package"))
	end
end

include("src/system_structs.jl")
include("src/solver/solver.jl")
include("src/solver/solver_z.jl")
include("src/utils.jl")
include("test/test_solver.jl")

function test()
	η0 = 0.2
	L = 3
	M = 30 + L
	v_ext = vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L))
	sys1 = make_RodLat(L, M, η0, v_ext)
	println("starting system:")
	println(sys1)
	
	
	α = 0.05
	steps = 2
	
	# sys4 = step_ρ(sys1, α)
	# for i in 1:steps-1
	#   sys4 = step_ρ(sys4, α)
	# end

	sys4, steps = solve_RodLat(sys1, 1e-9, 1e-1)
	
	# sys2 = sys1
	# sys3 = sys1
	# sys4 = sys1
	
	plot(eachindex(sys1.ρ[1:L*10]) |> collect, [sys1.ρ[1:L*10], sys4.ρ[1:L*10]], label=["initial state" "$steps steps"], legend=:bottom)
	vline!([L], label="Wall")
end


function test_each_function()
	L  = 3
	M  = 50
	η0 = 0.9

	v_ext = vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L))
	sys1 = make_RodLat(L, M, η0, v_ext)

	plt = plot(eachindex(sys1.ρ), log10.(sys1.ρ), label="ρ_start")
	n0 = n_0(sys1)
	n1 = n_1(sys1, n0)
	plot!(plt, eachindex(n0), [log10.(n0), log10.(n1)], label=["n0" "n1"])
	# plot!(plt, eachindex(n0), log10.(μ1(n1, sys1.L)), label="μ1")
	# plot!(plt, eachindex(n0), log10.(μ2(n0, sys1.L)), label="μ2")
	plot!(plt, eachindex(n0), log10.(exp.(-1 .* μ1(n1, sys1.L) .- μ2(n0, sys1.L))), label="exp(-μ1-μ2)")
	ρ_new = sys1.ρ0*μ_ex_ρ0(sys1) .* exp.(-1 .* μ1(n1, sys1.L) .- μ2(n0, sys1.L)) .* sys1.v_ext
	plot!(plt, eachindex(n0), log10.(ρ_new), label="ρ_new")
	plot!(plt, eachindex(n0), log10.(ones(size(n0)) .* μ_ex_ρ0(sys1)), label="μ_ex_ρ0")
	display(plt)
	# picard iteration
	# ρ_new = (1-α) .* sys1.ρ .+ α .* ρ_new


end


function test_z(L=3, M=400, ηs=0.1:0.1:0.9, v_ext=vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L)))

	# Parameter des Systems:
	# L     = 3
	# M     = 1000 + L
	# ηs    = 0.1:0.1:0.9
	# v_ext = vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L))
	
	# Parameter für Genauigkeit bzw. für Plots
	plot_range = L:(10*L)
	α_minmax   = (0.01, 0.1)
	ϵ          = 5e-12

	# Generiere initiale Systeme
	systems = []
	for i in eachindex(ηs)
		push!(systems, make_RodLat(L, M, ηs[i], v_ext))
	end

	comp_times = zeros(Float64, size(ηs))
	# Berechne Systeme
	for i in eachindex(systems)
		starttime = time()
		println("Stats for Sys", i)
		# systems[i] = solve_RodLat(systems[i], α_minmax .* (1 - ηs[i]), ϵ)
		systems[i] =      ρ_steps(systems[i], α_minmax .* (1 - ηs[i]), ϵ)
		println("\n")
		comp_times[i] = time() - starttime
	end
	display(scatter(ηs, comp_times, xlabel="η0", ylabel="computation time"))
	
	
	# Compute surface tension and save it
	γs = γ.(systems)
	γs_ana = γ_analytical.(ηs ./ L, L)
	γ_df = DataFrame(density=ηs, surface_tension=γs)
	CSV.write("results/surface_tension_L_$(L).csv", γ_df)

	# plot surface tension
	γ_plot = scatter(ηs, γs, xlabel=L"\eta_0", ylabel=L"\gamma_0", label="numerical", title="Surface tension for different densities (L = $(L))", dpi=300)
	plot!(γ_plot, ηs, γs_ana, label="analytical")
	savefig(γ_plot, "results/surface_tension_L_$(L).png")


	# Compute gibbs_adsorption and save it
	
	Γs = Γ.(systems)
	Γs_ana = Γ_analytical.(ηs ./ L, L)
	Γ_df = DataFrame(density=ηs, gibbs_adsorption=Γs)
	CSV.write("results/gibbs_adsorption_L_$(L).csv", Γ_df)

	# plot gibbs adsorption
	Γ_plot = scatter(ηs, Γs, xlabel=L"\eta_0", ylabel=L"\Gamma", label="numerical", title="gibbs adsorption for different densities (L = $(L))", dpi=300)
	plot!(Γ_plot, ηs, Γs_ana, label="analytical")
	savefig(Γ_plot, "results/gibbs_adsorption_L_$(L).png")


	# Fertige Plots an
	colors = palette([:blue, :red, :orange], size(systems, 1)) |> collect |> permutedims
	if size(ηs, 1) <= 3
		colors = repeat([:auto], size(ηs, 1))
	end
	plot((eachindex(systems[1].ρ) |> collect)[plot_range], (x -> x.ρ[plot_range]).(systems), label = permutedims(L"\eta_0 = " .* string.(0.1:0.1:0.9)), xlabel = "s", ylabel = L"\rho", title = "L = " * string(L), legend = :topright, dpi = 300, foreground_color_legend = nothing, linecolor=colors)
	savefig("results/L_$(L)_Dichteprofil.png")

end

# test_z(3)
test_z(10)



# TODO-List
# - Gamma analytisch mit plotten ✓
# - Gibbs analytisch mit plotten ✓
# - integrieren von gamma um surface tension zu bestimmen. ⨯