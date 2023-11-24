# Berechnung der n_s(0) bzw. n_(1)
function n(start::Integer, stop::Integer, ρs::Vector{Float64})
    
    len_ρs =size(ρs, 1)

    sum = 0
    for i in max(1, start):min(len_ρs, stop)
        sum += ρs[i]
    end
    return sum

end


# Berechne µ_ex
function μ_ex(s::Integer, L::Integer, ρs::Vector{Float64})
   len_ρs = size(ρs, 1)
   µs1 = 0
   µs2 = 0

   # Berechne µs1
   for i in (s-L+1):s
    µs1 += -log(1-n((i-L+1), i, ρs))
   end

   # Berechne µs2
   for i in (s-L+11):(s-1)
    println("prob hier")
    -log(1-n((i-L+1), (i-1), ρs))
   end

   # Gebe µ_ex zurück
   return µs1 - µs2 
end


# Berechne den Vorfaktor für ρs
function prefactor(L::Integer, ρ0::Float64)
    return ρ0 * exp(L * log(1 + (ρ0/(1 - L * ρ0))) - log(1 + ρ0 + L * ρ0))
end


function ρ_step(sys::RodLat, α::Float64)
    
    # Ziehe Systemdaten
    L     = sys.L
    M     = sys.M
    ρ0    = sys.ρ0
    v_ext = sys.v_ext
    ρs    = sys.ρ
    test=deepcopy(ρs)

    # Ziehe Parameter für Berechnung der neuen Dichte-Werte nach Picard & mixing
    len_ρs = size(ρs, 1)
    ρs_new = deepcopy(ρs)

    # Berechne für jeden index s in ρs den neuen Wert 
    for s in 1:len_ρs
        value    = prefactor(L, ρ0) * exp(-1 * μ_ex(s, L, ρs)) * v_ext[s]
        ρs_new[s] = (1 - α) * ρs[s] + α * value
    end

    # Update System
    sys.ρ = ρs_new

    println(test .- ρs_new)

    # Gebe System zurück
    return sys
end