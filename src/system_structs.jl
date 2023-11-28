mutable struct RodLat
  ρ::Vector{Float64}
  L::Int
  M::Int
  ρ0::Float64
  v_ext::Vector{Float64}
  β::Float64
  η0::Float64
  
  function RodLat(ρ::Vector{Float64}, L::Int, M::Int, η0::Float64, v_ext::Vector{Float64})
    ρ0 = η0/L
    if M < 2*L #FIXME correct condition?
      error("lattice to small: M=$M; L=$L")
    elseif ρ0 > 1 || ρ0 < 0
      error("ρ0 not valid: $ρ0. Must be ∈ [0,1]")
    end
    new(ρ, L, M, ρ0, v_ext::Vector{Float64}, 1.0, η0)
  end
end

function make_RodLat(L::Int, M::Int, η0::Float64, v_ext::Vector{Float64})
  RodLat(η0/L .* v_ext, L, M, η0, v_ext)
end