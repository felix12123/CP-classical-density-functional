Φ_OD_deriv(x) = -log(max(0.0, 1-x))
function Φ_OD(x::Float64)::Float64
	x + (1 - x)*log(1 - x)
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
		return min.(0.999999, n1)
		# display(plot([n1, sys.ρ]))
		# error("n_1 has elements larger than 1: $(maximum(n1))")
	elseif minimum(n1) < 0
		error("n_1 has elements smaller than 0: $(minimum(n1))")
	end
	return n1
end

function μ1(n1::Vector{Float64}, L::Int)::Vector{Float64}
	return [sum(Φ_OD_deriv.(n1[s:min(end, s+L-1)])) for s in eachindex(n1)]
end

function μ2(n0::Vector{Float64}, L::Int)::Vector{Float64}
	[sum(Φ_OD_deriv.(n0[min(end, s+1):min(end, s+L-1)])) for s in eachindex(n0)]
end

function μ_ex(sys::RodLat)
	n0 = n_0(sys)
	n1 = n_1(sys, n0)

	return μ1(n1, sys.L) - μ2(n0, sys.L)
end

function μ_ρ0(ρ0::Float64, L::Int=3)::Float64
	L * log(1 + (ρ0/ (1 - L * ρ0))) - log(1 + ρ0 - L * ρ0) + log(ρ0)
end
function μ_ρ0(sys::RodLat)::Float64
	# log(sys.ρ0) + μ_ex(div(size(sys.ρ, 1), 2), sys.L, sys.ρ)
	L, ρ0 = sys.L, sys.ρ0
	L * log(1 + (ρ0/ (1 - L * ρ0))) - log(1 + ρ0 - L * ρ0) + log(ρ0)
end