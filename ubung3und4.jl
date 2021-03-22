# Simulationsparameter
m = 0.4     # Masse des freien Teilchens
d_t = 1     # Räumliche Dimensionen
d_s = 1     # Temporale Dimensionen
N_t = 8     # Anzahl räumlicher Gitterpunkte
N_s = 8     # Anzahl temporaler Gitterpunkte
m = 40      # Anzahl Wiederholungen in Metropolis
r = 0.2     # Schrittweite in Metropolis

# Startkonfiguration generieren
# Der erste Index steht für temporale, der zweite für räumliche Koordinaten
x_0 = []
for i = 1:N_t
    q = 2*(rand(N_s) - 0.5*ones(N_s))
    push!(x_0, q)
end
x_0 = 10*x_0

# Funktion um Wirkung einer Konfiguration auszurechnen (vgl. Gl(56) im Skript)
# Beachte die periodischen Randbedingungen
function S(x)
    a = 0
    for i = 1:N_t-1
        for j in 1:N_s-1
            a = a + (0.5*m^2 + d_s+d_t) * x[i][j]^2 - x[i][j] * (x[i+1][j] + x[i][j+1])
        end
    end
    a = a + (0.5*m^2 + d_s+d_t) * x[N_t][N_s]^2 - x[N_t][N_s] * (x[1][N_s] + x[N_t][1])
    a
end

# Der Metropolis Algorithmus
configs = [x_0]
accepted = 0
for i = 1:m
    q = []
    for i = 1:N_t
        b = 2*(rand(N_s) - 0.5*ones(N_s))
        push!(q, b)
    end

    proposal = configs[i] + r*q
    Delta = S(proposal) - S(configs[i])
    a = rand(1)
    if a[1] < exp(-Delta)
        push!(configs, proposal)
        accepted = accepted + 1
    else
        push!(configs, configs[i])
    end
end

actions = []
for i = 1:m+1
    push!(actions, S(configs[i]))
end

println("Es wurden $accepted von $m Updates akzepiert, also ca. ", 100*accepted/m, "%" )


# Funktion für diskrete Fourier Transformation einer Konfiguration x[l],
# wobei n ∈ {0,1,...,N_s/2} und τ (tau) die temporale Koordinate
function discrete_fourier(l, tau, n)
    a = 0
    p = n*(2*pi/N_s)
    for i = 1:N_s
        a = a + exp(-im*p*configs[l][tau][i]) * configs[l][tau][i]
    end
    a
end
