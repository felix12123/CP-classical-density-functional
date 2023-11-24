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
include("src/solver/solver_z.jl")
include("src/utils.jl")
include("test/test_solver.jl")

function test()
  η0 = 0.1
  L = 3
  M = 30 + L
  v_ext = vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L))
  sys1 = make_RodLat(3, 20, η0, v_ext)
  println(sys1)
  α = 1.0
  sys2 = step_ρ(sys1, α)
  println(sys2)
  sys3 = step_ρ(sys2, α)
  sys4 = step_ρ(sys3, α)
  for i in 1:100
    sys4 = step_ρ(sys4, α)
  end
  
  # sys2 = sys1
  # sys3 = sys1
  # sys4 = sys1

  plot(eachindex(sys1.ρ) |> collect, [sys1.ρ, sys2.ρ, sys3.ρ, sys4.ρ], label=permutedims(string.(0:3|>collect) .* " steps"), legend=:bottom)
  vline!([L, M], label="Wall")
end

function test_z()

  # Lege Parameter und externes Potential fest
  η0 = 0.1
  L  = 3
  M  = 30 + L
  α  = 1.0
  v_ext = vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L))

  # Generiere System zum testen
  it_0 = make_RodLat(L, M, η0, v_ext)
  #println(it_0)

  # Testsystem nach 1 Iteration
  it_1 = ρ_step(it_0, α)
  it_2 = ρ_step(deepcopy(it_1), α)
  it_3 = ρ_step(deepcopy(it_2), α)
  it_4 = ρ_step(deepcopy(it_3), α)

  plot(eachindex(it_1.ρ) |> collect, [it_1.ρ, it_2.ρ, it_3.ρ, it_4.ρ], label=permutedims(string.(0:3|>collect) .* " steps"), legend=:bottom)
  vline!([L, M], label="Wall")
end

test_z()

#test()