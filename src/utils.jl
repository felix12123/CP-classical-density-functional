
import Base.println
function println(sys::RodLat)
  println("L  = $(sys.L); \t M  = $(sys.M)")
  println("η0 = $(sys.η0);\t ρ0 = $(sys.ρ0)")
end