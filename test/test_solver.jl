using Test


function test_solver()
  @testset "solver" begin
    v_ext = vcat(zeros(Float64, L), ones(Float64, M-L-1), zeros(Float64, L))
    test_sys = make_RodLat(3, 33, 0.2, v_ext)
    # @test 
  end
end