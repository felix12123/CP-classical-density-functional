function prefactor(sys::RodLat)::Float64
	return sys.ρ0 * exp(sys.L * log(1 + (sys.ρ0/ (1 - sys.L * sys.ρ0))) - log(1 + sys.ρ0 - sys.L * sys.ρ0))
end

function ρ_step!(sys::RodLat, α::Float64)
	ρs_new_i = prefactor(sys) .* exp.(-1 .* μ_ex(sys)) .* sys.v_ext
	ϵ = sum((ρs_new_i .- sys.ρ) .^ 2) 										# calculate difference to previous state
	sys.ρ .= (1 - α) .* sys.ρ .+ α .* ρs_new_i 						# picard iteration
	sys.ρ0 = sys.ρ[div(end, 2)]
	sys.η0 = sys.ρ0 * sys.L
	return ϵ
end

function solve_RodLat(sys::RodLat, α_min_max::NTuple{2, Float64}, ϵ_crit::Float64, max_iter::Int=200_000; print_res=true)
	sys     = deepcopy(sys)
	counter = 0
	ϵ_old   = Inf
	α       = α_min_max[1]
	ϵ_new   = Inf
	while !(ϵ_new < ϵ_crit || counter >= max_iter)		# interrutp program after so many iterations to prevent endless loop, remove for more accurate calculations
		counter += 1

		# Update schritt für das system
		ϵ_new = ρ_step!(sys, α)

		# Variables α:
		if ϵ_new < ϵ_old
			α = min(α_min_max[2], 1.1 * α)
		else
			α = max(α_min_max[1], α/5)
		end
		ϵ_old = ϵ_new
	end

	if print_res
		println("Konvergenz-Summe = ", ϵ_new)
		println("Iterationen bis zur Konvergenz: ", counter)
	end
	
	# Gebe System aus
	return sys
end