# Berechnung der n_s(0) bzw. n_(1)
function n(start::Integer, stop::Integer, ρs::Vector{Float64})::Float64
	
	len_ρs = size(ρs, 1)

	sum = 0
	for i in max(1, start):min(len_ρs, stop)
		sum += ρs[i]
	end
	return sum

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

function n_0(sys::RodLat)::Vector{Float64}
	L, ρ = sys.L, sys.ρ
	N = size(ρ, 1)

	[sum(ρ[max(1, s-L+1):min(N, s-1)]) for s in eachindex(ρ)]
end

# you can input n_0 if it was already calculated, to save a little bit of time.
function n_1(sys::RodLat, n_0::Vector{Float64}=[NaN])::Vector{Float64}
	n1 = zeros(Float64, size(sys.ρ))
	if isequal(n_0, [NaN])
		n1 = [sum(sys.ρ[max(1, s-sys.L+1):s]) for s in eachindex(sys.ρ)]
	else
		n1 = n_0 .+ sys.ρ
	end

	if maximum(n1) > 1
		display(plot([n1, sys.ρ]))
		error("n_1 has elements larger than 1: $(maximum(n1))")
	elseif minimum(n1) < 0
		error("n_1 has elements smaller than 0: $(minimum(n1))")
	end
	return n1
end

function Φ_OD(x::Float64)::Float64
	x + (1 - x)*log(1 - x)
end

function μ_ρ0(sys::RodLat)::Float64
	log(sys.ρ0) + μ_ex(div(size(sys.ρ, 1), 2), sys.L, sys.ρ)
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

function γ_from_int(ρ0, L)
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
		else
			α = max(α_min_max[1], α/5)
		end
		ϵ_old = ϵ_new
		
		# Update System
		if ϵ_new < ϵ_crit || counter > 300_000
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
			convergence_array[s] = (ρs_i1[s] - ρs[s])^2
		end

		# Variables α:
		if sum(convergence_array) < ϵ_old
			α = min(α_min_max[2], 1.1 * α)
		else
			α = max(α_min_max[1], α/5)
		end
		ϵ_old = sum(convergence_array)
		
		# Update System
		if sum(convergence_array) < ϵ_crit || counter > 300_000
			convergent = true
			println("Konvergenz-Summe = ", sum(convergence_array))
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