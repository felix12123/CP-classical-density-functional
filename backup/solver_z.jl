# Berechnung der n_s(0) bzw. n_(1)
function n(start::Integer, stop::Integer, ρs::Vector{Float64})::Float64
	
	len_ρs = size(ρs, 1)

	sum = 0
	for i in max(1, start):min(len_ρs, stop)
		sum += ρs[i]
	end
	return sum

end


Φ_OD_deriv(x) = -log(1-x)

function μ1(n1::Vector{Float64}, L::Int)::Vector{Float64}
	return [sum(Φ_OD_deriv.(n1[s:min(end, s+L-1)])) for s in eachindex(n1)]
end
function μ1(sys::RodLat)::Vector{Float64} # can probably be deleted
	n1 = n_1(sys)
	L = sys.L
	[sum(Φ_OD_deriv.(n1[s:min(end, s+L-1)])) for s in eachindex(n1)]
end


function μ2(n0::Vector{Float64}, L::Int)::Vector{Float64}
	[sum(Φ_OD_deriv.(n0[min(end, s+1):min(end, s+L-1)])) for s in eachindex(n0)]
end
function μ2(sys::RodLat)::Vector{Float64} # can probably be deleted
	n0 = n_0(sys)
	L = sys.L
	[sum(Φ_OD_deriv.(n0[min(end, s+1):max(end, s+L-1)])) for s in eachindex(n0)]
end

function μ_ex(sys::RodLat)
	n0 = n_0(sys)
	n1 = n_1(sys, n0)

	return μ1(n1, sys.L) - μ2(n0, sys.L)
end

# Berechne µ_ex
function μ_ex(s::Integer, L::Integer, ρs::Vector{Float64})::Float64

	res = 0

   # Berechne µs1
   for i in s:(s+L-1)
	res += -log(1-n((i-L+1), i, ρs))
   end

   # Berechne µs2
   for i in (s+1):(s+L-1)
	res -= -log(1-n((i-L+1), (i-1), ρs))
   end
   
   return res
   
end



# Berechne den Vorfaktor für ρs
function prefactor(L::Integer, ρ0::Float64)::Float64
	return ρ0 * exp(L * log(1 + (ρ0/ (1 - L * ρ0))) - log(1 + ρ0 - L * ρ0))
end
function prefactor(sys::RodLat)::Float64
	return sys.ρ0 * exp(sys.L * log(1 + (sys.ρ0/ (1 - sys.L * sys.ρ0))) - log(1 + sys.ρ0 - sys.L * sys.ρ0))
end



function ρ_step(sys::RodLat, α::Float64)
	
	# Ziehe Systemdaten
	L     = sys.L
	ρ0    = sys.ρ0
	v_ext = sys.v_ext
	ρs    = sys.ρ

	# Ziehe Parameter für Berechnung der neuen Dichte-Werte nach Picard & mixing
	len_ρs            = size(ρs, 1)
	ρs_i1             = deepcopy(ρs)
	convergence_array = ones(Float64, len_ρs)

	# Berechne für jeden index s in ρs den neuen Wert 
	pref = prefactor(sys)
	for s in 1:len_ρs
		ρs_new_i = pref * exp(-1 * μ_ex(s, L, ρs)) * v_ext[s]
		convergence_array[s] = (ρs_new_i - ρs[s])^2
		ρs_i1[s] = (1 - α) * ρs[s] + α * ρs_new_i
	end

	# Update System
	sys.ρ = ρs_i1
	sys.ρ0 = sys.ρ[div(end, 2)]
	sys.η0 = sys.ρ0 * sys.L

	# Gebe System aus
	return sys

end


function μ_ρ0(sys::RodLat)::Float64
	# log(sys.ρ0) + μ_ex(div(size(sys.ρ, 1), 2), sys.L, sys.ρ)
	L, ρ0 = sys.L, sys.ρ0
	L * log(1 + (ρ0/ (1 - L * ρ0))) - log(1 + ρ0 - L * ρ0) + log(ρ0)
end
function μ_ρ0_analytical(ρ0::Float64, L::Int=3)::Float64
	L * log(1 + (ρ0/ (1 - L * ρ0))) - log(1 + ρ0 - L * ρ0) + log(ρ0)
end

function ω(sys::RodLat)::Vector{Float64}
	(sys.ρ .* (log.(sys.ρ) .- 1) .+ Φ_OD.(n_1(sys)) .- Φ_OD.(n_0(sys)) - μ_ρ0(sys) .* sys.ρ)
end

function p(sys::RodLat)::Vector{Float64}
	log.((1 .- sys.ρ .* (sys.L - 1)) ./ (1 .- sys.ρ .* sys.L))
end
function p(ρ::Float64, L::Int)::Float64
	log((1-ρ*(L-1)) / (1-ρ*L))
end

function γ(sys)
	0.5 * (ω(sys) .+ p(sys.ρ0, sys.L))[sys.L+1:end-sys.L] |> sum
end

function γ_analytical(ρ0::Float64, L::Int)
	0.5 * ((1-L) * log(1 + ρ0 / (1 - L * ρ0)) - log(1 + ρ0*(1 - L)))
end

function Γ_analytical(ρ0::Float64, L::Int, Δρ::Float64=1e-3)
	-(γ_analytical(ρ0 + Δρ, L) - γ_analytical(ρ0 - Δρ, L)) / (μ_ρ0_analytical(ρ0 + Δρ, L) - μ_ρ0_analytical(ρ0 - Δρ, L))
end


function Γ(sys::RodLat)
	0.5 * sum((sys.ρ .- sys.ρ0)[sys.L+1:end-sys.L])
end

function ρ_step!(sys::RodLat, α::Float64)
	ρs_new_i = prefactor(sys.L, sys.ρ0) .* exp.(-1 .* μ_ex(sys)) .* sys.v_ext
	ϵ = sum((ρs_new_i .- sys.ρ) .^ 2)
	sys.ρ .= (1 - α) .* sys.ρ .+ α .* ρs_new_i
	sys.ρ0 = sys.ρ[div(end, 2)]
	sys.η0 = sys.ρ0 * sys.L
	return ϵ
end

function solve_RodLat(sys::RodLat, α_min_max::NTuple{2, Float64}, ϵ_crit::Float64)
	sys = deepcopy(sys)
	counter = 0
	# prev_ρ = deepcopy(sys.ρ)
	ϵ_old  = Inf
	α      = α_min_max[1]

	while true
		counter += 1
		# prev_ρ = sys.ρ |> deepcopy

		# Update schritt für das system
		ϵ_new = ρ_step!(sys, α)

		# Variables α:
		if ϵ_new < ϵ_old
			α = min(α_min_max[2], 1.1 * α)
			# println("α got bigger")
		else
			α = max(α_min_max[1], α/5)
			# println("α got smaller")
		end
		ϵ_old = ϵ_new
		
		# Update System
		if ϵ_new < ϵ_crit || counter > 75_000
			println("Konvergenz-Summe = ", ϵ_new)
			println("Iterationen bis zur Konvergenz: ", counter)
			break
		end
	end
	
	# Gebe System aus
	return sys
end


function ρ_steps(sys::RodLat, α_min_max::NTuple{2, Float64}, ϵ_crit::Float64)
	
	convergent = false
	convergence_array = ones(Float64, size(sys.ρ, 1))
	counter = 1
	
	while !convergent
		counter += 1
		# Ziehe wichtige Daten
		ρ0     = sys.ρ0
		v_ext  = sys.v_ext
		ρs     = sys.ρ
		len_ρs = size(ρs, 1)
		L      = sys.L
		ρs_i1  = deepcopy(ρs)
		ϵ_old  = Inf
		α      = α_min_max[1]
		
		# Berechne für jeden index s in ρs den neuen Wert 
		for s in 1:len_ρs
			ρs_new_i = prefactor(L, ρ0) * exp(-1 * μ_ex(s, L, ρs)) * v_ext[s]
			ρs_i1[s] = (1 - α) * ρs[s] + α * ρs_new_i
			convergence_array[s] = (ρs_new_i - ρs[s])^2
		end
		ϵ_new = sum(convergence_array)
		
		# Variables α:
		if ϵ_new < ϵ_old
			α = min(α_min_max[2], 1.1 * α)
		else
			α = max(α_min_max[1], α/5)
			println("α got smaller")
		end
		ϵ_old = ϵ_new
		
		# Update System
		if ϵ_new < ϵ_crit || counter > 75_000
			convergent = true
			println("Konvergenz-Summe = ", ϵ_new)
			println("Iterationen bis zur Konvergenz: ", counter)
			# println("γ = ", γ(sys))
		end
		sys.ρ = ρs_i1
		sys.ρ0 = sys.ρ[div(end, 2)]
		sys.η0 = sys.ρ0 * sys.L
	end

	# Gebe System aus
	return sys

end