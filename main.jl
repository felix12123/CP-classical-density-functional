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
Packages = ["Plots", "Test"]
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
include("src/utils.jl")
include("test/test_solver.jl")

function test()
  η0 = 0.995
  L = 4
  M = 350
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
# test_each_function()
test()