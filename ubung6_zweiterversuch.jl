using Plots
using Statistics
#=
Die Phasen von Eichfeld-Konfigurationen Uμ( ⃗n) werden in einem Array
config[μ, nₜ, nₛ] gespeichert, wobei ⃗n = (nₜ, nₛ). Im weiteren soll
μ = 1 für die zeitliche und μ = 2 für die räumliche Dimension stehen.
(es wird auch μ ≡ 1 für temp. und ν ≡ 2 für räuml. synonym verwendet)
Phasen der Eichtransformationen Ω( ⃗n) werden als trans[nₜ, nₛ] gespeichert.
=#
N_t = 12        # Anzahl temporaler Gitterpunkte
N_s = 12    	# Anzahl räumlicher Gitterpunkte
beta = 1.8      # 1/g², wobei g die Eich-Kopplung ist
m = 100000          # Updates im Metropolis
r = 2*pi*0.1      # Schrittweite im Metropolis (nur in pos. Richtung, s.u.)
cut = Int64(0.6*m)


# Die erste Feldkonfig. wird per Zufall generiert, aber mit Werten ∈ [0,2π]
config_1 = []
for mu = 1:2
    push!(config_1, [])
    for i = 1:N_t
        a = (2*pi).*rand(N_s)
        push!(config_1[mu], a)
    end
end

#=
# Diese Staple Funktion diente nur zum Vergleich mit der unteren, in der
# keine if-condition verwendet wurde. Eine Verbesserung konnte nach
# halb-professionellem Ausprobieren nicht festgestellt werden
function staple_dagger(mu, nt, ns, config)
    a = 0
    b = 0

    NT = nt%N_t +1
    NS = ns%N_s +1
    TN = (nt + N_t -2)%N_t +1
    SN = (ns + N_s -2)%N_s +1

    if mu == 1
        a = -config[2][nt][ns] - config[1][nt][NS] + config[2][NT][ns]
        b = config[2][nt][SN] - config[1][nt][SN] - config[2][NT][SN]
    elseif mu == 2
        a = -config[1][nt][ns] - config[2][NT][ns] + config[1][nt][NS]
        b = config[1][TN][ns] - config[2][TN][ns] - config[1][TN][NS]
    end

    exp(im*a) + exp(im*b)
end
=#


# Die Lesbarkeit der Indizes c-j ist ein wenig der Effizienz zum Opfer gefallen.
# Daher steht die jeweilige Bedeutung neben jedem Index im Kommentar.
function staple_dagger(mu, nt, ns, config)
    MU = mu%2 +1    # MU ∈ {1,2}, MU ≠ mu
    V = (mu+1)%2

    c = (nt + V -1)%N_t +1   # nt für mu = 1, sonst nt+1 (mit per. Randb.)
    d = (ns + MU -2)%N_s +1  # ns+1 für mu = 1, sonst ns
    e = (nt + MU -2)%N_t +1  # nt+1 für mu = 1, sonst nt
    f = (ns + V -1)%N_s +1   # ns für mu = 1, sonst ns+1

    g = (nt + N_t - V -1)%N_t +1    # nt für mu = 1, sonst nt-1
    h = (ns + N_s - MU)%N_s +1      # ns-1 für mu = 1, sonst ns
    i = Int64((nt + 2*(MU -1.5) +N_t -1)%N_t +1)    # nt+1 für mu = 1, sonst nt-1
    j = Int64((ns + 2*(V -0.5) +N_s -1)%N_s +1)     # ns-1 für mu = 1, sonst ns+1

    a = -config[MU][nt][ns] - config[mu][c][d] + config[MU][e][f]
    b = config[MU][g][h] - config[mu][g][h] - config[MU][i][j]

    exp(im*a) + exp(im*b)
end


link_accepted = 0
raw_configs = [config_1]
for k = 1:m
    test = deepcopy(raw_configs[k])

    for i = 1:N_t
        for j = 1:Int64(0.5*N_s)
            J = 2*j
            a = test[1][i][J]
            b = staple_dagger(1,i,J,test)
            proposal = deepcopy(a)
            q = 2*(rand(1)[1] -0.5)
            proposal += r*q

            Delta = -beta * real((exp(im*proposal)- exp(im*a)) * b)
            P = rand(1)[1]
            if P < exp(-Delta)
                test[1][i][J] = proposal
                link_accepted += 1
            end
        end
    end

    for i = 1:N_t
        for j = 1:Int64(0.5*N_s)
            J = 2*j-1
            a = test[1][i][J]
            b = staple_dagger(1,i,J,test)
            proposal = deepcopy(a)
            q = 2*(rand(1)[1] -0.5)
            proposal += r*q

            Delta = -beta * real((exp(im*proposal)- exp(im*a)) * b)
            P = rand(1)[1]
            if P < exp(-Delta)
                test[1][i][J] = proposal
                link_accepted += 1
            end
        end
    end

    for j = 1:N_s
        for i = 1:Int64(0.5*N_t)
            I = 2*i
            a = test[2][I][j]
            b = staple_dagger(2,I,j,test)
            proposal = deepcopy(a)
            q = 2*(rand(1)[1] -0.5)
            proposal += r*q

            Delta = -beta * real((exp(im*proposal)- exp(im*a)) * b)
            P = rand(1)[1]
            if P < exp(-Delta)
                test[2][I][j] = proposal
                link_accepted += 1
            end
        end
    end

    for j = 1:N_s
        for i = 1:Int64(0.5*N_t)
            I = 2*i-1
            a = test[2][I][j]
            b = staple_dagger(2,I,j,test)
            proposal = deepcopy(a)
            q = 2*(rand(1)[1] -0.5)
            proposal += r*q

            Delta = -beta * real((exp(im*proposal)- exp(im*a)) * b)
            P = rand(1)[1]
            if P < exp(-Delta)
                test[2][I][j] = proposal
                link_accepted += 1
            end
        end
    end

    push!(raw_configs, test)
end



accepted = 0
for k = 2:m+1
    if raw_configs[k] != raw_configs[k-1]
        accepted += 1
    end
end


println("", "Es wurden $accepted von $m Updates akzeptiert, also ",
100*(accepted/m), "%")
println("", "Es wurden $link_accepted von ", m*2*N_t*N_s, " Link-Updates akzeptiert, also ",
100*(link_accepted/(m*2*N_t*N_s)), "%")

# Funktion, um die Phase einer Plaquette Pμν( ⃗n) zu berechnen
function P12(nt, ns, config)
    NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
    NS = ns%N_s +1  # nₛ+1 mit periodische Randbed.
    config[1][nt][ns] + config[2][NT][ns] - config[1][nt][NS] - config[2][nt][ns]
end

Plaquettes_avg = []
for k = cut:m+1
    a = 0
    for i = 1:N_t
        for j = 1:N_s
            a += exp(im*P12(i, j, raw_configs[k]))
        end
    end
    a = (1/(N_t*N_s)) * a
    push!(Plaquettes_avg, a)
end

println("", "And the average plaquette is ", real(mean(Plaquettes_avg)), " ± ",
std(real.(Plaquettes_avg)))
println("  ")

# Für m = 10⁵:
# Simulationsdauer if-staple:       53 sek.
# Ergebnis if-staple:               0.6615 ± 0.0369

# Simulationsdauer anderer Staple:  73 sek.
# Ergebnis anderer staple:          0.6613 ± 0.0361

# Simulationsdauer aktu. Staple:    80 sek.
# Ergebnis anderer aktu. Staple:    0.6612 ± 0.0369
