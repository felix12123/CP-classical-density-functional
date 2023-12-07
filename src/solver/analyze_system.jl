function p(sys::RodLat)::Vector{Float64}
	log.((1 .- sys.ρ .* (sys.L - 1)) ./ (1 .- sys.ρ .* sys.L))
end
function p(ρ::Float64, L::Int)::Float64
	log((1-ρ*(L-1)) / (1-ρ*L))
end

function Γ(sys::RodLat)
	0.5 * sum((sys.ρ .- sys.ρ0)[sys.L+1:end-sys.L])
end

function γ_analytical(ρ0::Float64, L::Int)
	0.5 * ((1-L) * log(1 + ρ0 / (1 - L * ρ0)) - log(1 + ρ0*(1 - L)))
end

function Γ_analytical(ρ0::Float64, L::Int, Δρ::Float64=1e-3)
	-(γ_analytical(ρ0 + Δρ, L) - γ_analytical(ρ0 - Δρ, L)) / (μ_ρ0(ρ0 + Δρ, L) - μ_ρ0(ρ0 - Δρ, L))
end

function γ(sys)
	0.5 * (ω(sys) .+ p(sys.ρ0, sys.L))[sys.L+1:end-sys.L] |> sum
end

function ω(sys::RodLat)::Vector{Float64}
	(sys.ρ .* (log.(sys.ρ) .- 1) .+ Φ_OD.(n_1(sys)) .- Φ_OD.(n_0(sys)) - μ_ρ0(sys) .* sys.ρ)
end
