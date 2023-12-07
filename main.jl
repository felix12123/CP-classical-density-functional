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
Packages = ["Plots", "LaTeXStrings", "CSV", "DataFrames", "Statistics", "Measurements", "LsqFit", "Latexify"]
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
include("src/solver/basic_functions.jl")
include("src/solver/solver.jl")
include("src/solver/analyze_system.jl")


function start_evaluation(L=3, M=250*L, ηs=0.1:0.1:0.9 |> collect, v_ext=vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L)))
	# Parameter für Genauigkeit bzw. für Plots
	α_minmax   = (0.01, 0.5)
	ϵ          = 1e-12
	println("ϵ = $ϵ")
	println("α_minmax = $α_minmax")
	println("M = $M")


	# Generiere initiale Systeme
	systems = Vector{RodLat}(undef, size(ηs))
	for i in eachindex(ηs)
		systems[i] = make_RodLat(L, M, ηs[i], v_ext)
	end

	# we want to compare the computation time for each density
	comp_times = zeros(Float64, size(ηs))

	# Berechne Systeme
	for i in eachindex(systems)
		starttime = time()
		println("Stats for Sys", i, " η0 = $(ηs[i])")
		systems[i] = solve_RodLat(systems[i], α_minmax .* (1 - ηs[i]), ϵ, 300_000)
		comp_times[i] = time() - starttime
		println()
	end

	# save figure of the computation time
	savefig(scatter(ηs, comp_times, xlabel="η0", ylabel="computation time (s)"), "results/computation_timesL_$L")
	
	# create further visualisations
	analyze_results(systems, ηs)
end

function analyze_results(systems::Vector{RodLat}, ηs::Vector{Float64})
	# choose region to plot
	L = systems[1].L
	plot_range = L:(10*L)


	# Compute surface tension and save it
	γs = γ.(systems)
	γs_ana = γ_analytical.(ηs ./ L, L)
	γ_df = DataFrame(density=ηs, surface_tension=γs)
	CSV.write("results/surface_tension_L_$(L).csv", γ_df)
	println("Error of numerical results of surface tension: ", mean(abs.((γs .- γs_ana) ./ γs_ana)))

	# plot surface tension
	γ_plot = scatter(ηs, γs, xlabel=L"\eta_0", ylabel=L"\gamma_0", label="numerical", title="Surface tension for different densities (L = $(L))", dpi=300)
	plot!(γ_plot, ηs, γs_ana, label="analytical")
	savefig(γ_plot, "results/surface_tension_L_$(L).png")


	# Compute gibbs_adsorption and save it
	Γs = Γ.(systems)
	Γs_ana = Γ_analytical.(ηs ./ L, L)
	Γ_df = DataFrame(density=ηs, gibbs_adsorption=Γs)
	CSV.write("results/gibbs_adsorption_L_$(L).csv", Γ_df)
	println("Error of numerical results of gibbs adsorption: ", mean(abs.((Γs .- Γs_ana) ./ Γs_ana)))

	# plot gibbs adsorption
	Γ_plot = scatter(ηs, Γs, xlabel=L"\eta_0", ylabel=L"\Gamma", label="numerical", title="gibbs adsorption for different densities (L = $(L))", dpi=300)
	plot!(Γ_plot, ηs, Γs_ana, label="analytical")
	savefig(Γ_plot, "results/gibbs_adsorption_L_$(L).png")


	# Fertige Plots an
	colors = palette([:blue, :red, :orange], size(systems, 1)) |> collect |> permutedims
	if size(ηs, 1) <= 3
		colors = repeat([:auto], size(ηs, 1))
	end
	plot((eachindex(systems[1].ρ) |> collect)[plot_range], (x -> x.ρ[plot_range]).(systems), label = permutedims(L"\eta_0 = " .* string.(ηs)), xlabel = "s", ylabel = L"\rho", title = "L = " * string(L), legend = :topright, dpi = 300, size=(600, 400) .* 1.1, foreground_color_legend = nothing, linecolor=colors)
	savefig("results/L_$(L)_Dichteprofil.png")
end


function examine_comp_time(L=3, M=250*L, ηs=0.4:0.01:0.95 |> collect, N=25)
	v_ext=vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L))
	# Parameter für Genauigkeit bzw. für Plots
	α_minmax   = (0.01, 0.5)
	ϵ          = 1e-5
	println("ϵ = $ϵ")
	println("α_minmax = $α_minmax")
	println("M = $M")


	# Generiere initiale Systeme
	systems = Vector{RodLat}(undef, size(ηs))
	for i in eachindex(ηs)
		systems[i] = make_RodLat(L, M, ηs[i], v_ext)
	end

	# we want to compare the computation time for each density
	comp_times = zeros(Float64, size(ηs))

	# Berechne Systeme
	for i in 1:N
		for i in eachindex(systems)
			starttime = time()
			solve_RodLat(systems[i], α_minmax .* (1 - ηs[i]), ϵ, 300_000, print_res=false)
			comp_times[i] += (time() - starttime) / N
		end
	end
	model(x, p) = p[1] ./ (p[2].- x) .^ p[3] .+ p[4]
	p0 = [5.0, 1.0, 3.0, 0.0]
	lower = [-Inf, -Inf, 1.0, -Inf]
	# weights = comp_times ./ sum(comp_times)
	fit = LsqFit.curve_fit(model, ηs, comp_times, p0, lower=lower)
	params = fit.param .± stderror(fit)
	fitdata = model(ηs, Measurements.value.(params))
	# save figure of the computation time
	plt1 = scatter(ηs, comp_times, label="Results", xlabel="η0", ylabel="computation time (s)", dpi=300, legend=:outerbottom)
	p1 = Measurements.value.(params)
	p1 = params
	p1 = round.(p1, digits=5)
	plot!(ηs, fitdata, label="Fit = $(p1[1])/(($(p1[2]) - x) ^$(p1[3])) + $(p1[4])")
	savefig(plt1, "results/computation_timesL_$L")
	println("Computation time scales to 1/x^$(params[3])")
end

# start_evaluation(3)
# start_evaluation(10)

examine_comp_time()
