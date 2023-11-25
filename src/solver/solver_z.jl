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

    res = 0

   # Berechne µs1
   for i in s:(s+L-1)
    res += -log(1-n((i-L+1), i, ρs))
   end

   # Berechne µs2
   for i in (s+1):(s+L-1)
    res -= -log(1-n((i-L+1), (i-1), ρs))
   end
   
   return res
   
end


# Alternativer Weg für die Berechnung von µ_ex

function μ_ex_alternative(s::Integer, L::Integer, ρs::Vector{Float64})
    
    res = -log(1 - n((s-L+1), s, ρs))

    # Berechne Summe
    for i in (s+1):(s+L-1)
        res += log((1 - n((i-L+1), (i-1), ρs)) / (1 - n((i-L+1), i, ρs)))
    end

    return res

end


# Berechne den Vorfaktor für ρs
function prefactor(L::Integer, ρ0::Float64)

    return ρ0 * exp(L * log(1 + (ρ0/ (1 - L * ρ0))) - log(1 + ρ0 - L * ρ0))
    #return (η0 / L) * exp(L * log(1 + (η0/ (L * (1 - η0)))) - log(1 + (η0 / L) - η0))

end


function ρ_step(sys::RodLat, α::Float64)
    
    # Ziehe Systemdaten
    L     = sys.L
    ρ0    = sys.ρ0
    v_ext = sys.v_ext
    ρs    = sys.ρ

    # Ziehe Parameter für Berechnung der neuen Dichte-Werte nach Picard & mixing
    len_ρs = size(ρs, 1)
    ρs_new = deepcopy(ρs)

    # Berechne für jeden index s in ρs den neuen Wert 
    for s in 1:len_ρs
        value     = prefactor(L, ρ0) * exp(-1 * μ_ex(s, L, ρs_new)) * v_ext[s]
        ρs_new[s] = (1 - α) * ρs_new[s] + α * value
    end

    # Update System
    sys.ρ = ρs_new
    sys.ρ0 = sys.ρ[div(end, 2)]
    sys.η0 = sys.ρ0 * sys.L

    # Gebe System zurück
    return sys

end