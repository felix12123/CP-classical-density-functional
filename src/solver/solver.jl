function solve_RodLat(sys::RodLat, steps::Int)

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

Φ_OD_deriv(x) = -log(1-x)

function μ1(n1::Vector{Float64}, L::Int)::Vector{Float64}
  return [sum(Φ_OD_deriv.(n1[s:min(end, s+L-1)])) for s in eachindex(n1)]
  # catch
  #   for s in eachindex(n1)
  #     println("\nArgument of log:")
  #     display(n1[s:min(end, s+L-1)])
  #     sum(Φ_OD_deriv.(n1[s:min(end, s+L-1)]))
  #   end
  # end
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
  N = size(sys.ρ)
  n0 = n_0(sys)
  n1 = n_1(sys, n0)

  return μ1(n1, sys.L) - μ2(n0, sys.L)
end

# TODO there is (probably) a faster analytical solution to this
# Qestion: should the return value be a vector, or a value in the middle?
function μ_ex_ρ0(sys::RodLat)::Float64
  sys0 = make_RodLat(sys.L, sys.M, sys.ρ0, sys.v_ext)
  μ_ex(sys0)[div(end, 2)]
end

function step_ρ(sys::RodLat, α::Float64)
  n0 = n_0(sys)
  n1 = n_1(sys, n0)
  ρ_new = sys.ρ0*μ_ex_ρ0(sys) .* exp.(-1 .* μ1(n1, sys.L) .- μ2(n0, sys.L)) .* sys.v_ext
  ρ_new = (1-α) .* sys.ρ .+ α .* ρ_new
  sys = deepcopy(sys)
  sys.ρ .= ρ_new

  # FIXME do we have to update ρ0 every time?
  # sys.ρ0 = sys.ρ[div(end, 2)]
  # sys.η0 = sys.ρ0 * sys.L
  return sys
end