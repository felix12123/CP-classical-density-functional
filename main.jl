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
  #η = 0.1:0.1:0.9
  L  = 10
  M  = 500 + L
  α  = 0.001
  ϵ  = 1e-10
  ϵ8 = 1e-9     # FIXME: Wieso ist das nicht ~10e-10?
  ϵ9 = 1e-5     # FIXME: Wieso ist das nicht ~10e-10?
  v_ext = vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L))

  # Generiere Systeme
  sys1 = make_RodLat(L, M, η1, v_ext)
  sys2 = make_RodLat(L, M, η2, v_ext)
  sys3 = make_RodLat(L, M, η3, v_ext)
  sys4 = make_RodLat(L, M, η4, v_ext)
  sys5 = make_RodLat(L, M, η5, v_ext)
  sys6 = make_RodLat(L, M, η6, v_ext)
  sys7 = make_RodLat(L, M, η7, v_ext)
  sys8 = make_RodLat(L, M, η8, v_ext)
  sys9 = make_RodLat(L, M, η9, v_ext)

  println("Stats for Sys1")
  sys1 = ρ_steps(sys1, α, ϵ)
  println("\n Stats for Sys2")
  sys2 = ρ_steps(sys2, α, ϵ)
  println("\n Stats for Sys3")
  sys3 = ρ_steps(sys3, α, ϵ)
  println("\n Stats for Sys4")
  sys4 = ρ_steps(sys4, α, ϵ)
  println("\n Stats for Sys5")
  sys5 = ρ_steps(sys5, α, ϵ)
  println("\n Stats for Sys6")
  sys6 = ρ_steps(sys6, α, ϵ)
  println("\n Stats for Sys7")
  sys7 = ρ_steps(sys7, α, ϵ)
  println("\n Stats for Sys8")
  sys8 = ρ_steps(sys8, α, ϵ8)
  println("\n Stats for Sys9")
  sys9 = ρ_steps(sys9, α, ϵ9)
  println("\n")

  #plot((eachindex(sys1.ρ) |> collect)[3:30], [sys1.ρ[3:30], sys2.ρ[3:30], sys3.ρ[3:30], sys4.ρ[3:30], sys5.ρ[3:30], sys6.ρ[3:30], sys7.ρ[3:30], sys8.ρ[3:30], sys9.ρ[3:30]], label=[L"\eta_0 = 0.1" L"\eta_0 = 0.2" L"\eta_0 = 0.3" L"\eta_0 = 0.4" L"\eta_0 = 0.5" L"\eta_0 = 0.6" L"\eta_0 = 0.7" L"\eta_0 = 0.8" L"\eta_0 = 0.9"], xlabel = "s", ylabel = L"\rho", title = "L = 3", legend = :outerright, dpi = 300)
  #savefig("L=3_Dichteprofil.png")
  plot((eachindex(sys1.ρ) |> collect)[10:100], [sys1.ρ[10:100], sys2.ρ[10:100], sys3.ρ[10:100], sys4.ρ[10:100], sys5.ρ[10:100], sys6.ρ[10:100], sys7.ρ[10:100], sys8.ρ[10:100], sys9.ρ[10:100]], label=[L"\eta_0 = 0.1" L"\eta_0 = 0.2" L"\eta_0 = 0.3" L"\eta_0 = 0.4" L"\eta_0 = 0.5" L"\eta_0 = 0.6" L"\eta_0 = 0.7" L"\eta_0 = 0.8" L"\eta_0 = 0.9"], xlabel = "s", ylabel = L"\rho", title = "L = 10", legend = :outerright, dpi = 300)
  savefig("L=10_Dichteprofil.png")
end

test_z()

#test()