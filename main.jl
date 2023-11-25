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
Packages = ["Plots", "Test", "LaTeXStrings"]
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
  η0 = 0.2
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
  η1 = 0.1
  η2 = 0.2
  η3 = 0.3
  η4 = 0.4
  η5 = 0.5
  η6 = 0.6
  η7 = 0.7
  η8 = 0.8
  η9 = 0.9
  L  = 3
  M  = 30 + L
  α  = 1.0
  v_ext = vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L))

  # Generiere System zum testen
  sys1 = make_RodLat(L, M, η1, v_ext)
  sys2 = make_RodLat(L, M, η2, v_ext)
  sys3 = make_RodLat(L, M, η3, v_ext)
  sys4 = make_RodLat(L, M, η4, v_ext)
  sys5 = make_RodLat(L, M, η5, v_ext)
  sys6 = make_RodLat(L, M, η6, v_ext)
  #sys7 = make_RodLat(L, M, η7, v_ext)
  #sys8 = make_RodLat(L, M, η8, v_ext)
  #sys9 = make_RodLat(L, M, η9, v_ext)

  # Generiere Entwicklung bzgl. verschiedener η
  A1 = ρ_step(sys1, α)
  A2 = ρ_step(sys2, α)
  A3 = ρ_step(sys3, α)
  A4 = ρ_step(sys4, α)
  A5 = ρ_step(sys5, α)
  A6 = ρ_step(sys6, α)
  #A7 = ρ_step(sys7, α)
  #A8 = ρ_step(sys8, α)
  #A9 = ρ_step(sys9, α)

  plot((eachindex(A1.ρ) |> collect)[3:30], [A1.ρ[3:30], A2.ρ[3:30], A3.ρ[3:30], A4.ρ[3:30], A5.ρ[3:30], A6.ρ[3:30]], label=[L"\eta_0 = 0.1" L"\eta_0 = 0.2" L"\eta_0 = 0.3" L"\eta_0 = 0.4" L"\eta_0 = 0.5" L"\eta_0 = 0.6"], xlabel = "s", ylabel = L"\rho", title = "L = 3", legend = :topright, dpi = 300)
  #plot((eachindex(A1.ρ) |> collect)[3:30], [A1.ρ[3:30], A2.ρ[3:30], A3.ρ[3:30], A4.ρ[3:30], A5.ρ[3:30], A6.ρ[3:30], A7.ρ[3:30], A8.ρ[3:30], A9.ρ[3:30]], label=[L"\eta_0 = 0.1" L"\eta_0 = 0.2" L"\eta_0 = 0.3" L"\eta_0 = 0.4" L"\eta_0 = 0.5" L"\eta_0 = 0.6" L"\eta_0 = 0.7" L"\eta_0 = 0.8" L"\eta_0 = 0.9"], xlabel = "s", ylabel = L"\rho", title = "L = 3", legend = :topright, dpi = 300)
end

test_z()

#test()