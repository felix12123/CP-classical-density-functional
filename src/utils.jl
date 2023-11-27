
import Base.println
function println(sys::RodLat)
  println("L  = $(sys.L); \t M  = $(sys.M)")
  println("η0 = $(sys.η0);\t ρ0 = $(sys.ρ0)")
end

function visualize_RodLat(sys::RodLat)
  x  = eachindex(sys.ρ)
  y1 = sys.ρ
  plot(x, y1)
end