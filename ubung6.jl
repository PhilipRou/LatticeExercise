using Plots
using StatsBase
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
m = 5000          # Updates im Metropolis
r = 2*pi*0.05      # Schrittweite im Metropolis (nur in pos. Richtung, s.u.)
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


# Funktion um den adj. Staple A†μ( ⃗n) zu berechnen
function staple_dagger(mu, nt, ns, config)
    MU = mu%2 +1    # MU ∈ {1,2}, MU ≠ mu
    NT = nt%N_t +1  # ⃗n + ̂μ  mit periodischer Randbed.
    NS = ns%N_s +1  # ⃗n + ̂ν
    TN = (nt+N_t-2)%N_t +1  # ⃗n - ̂μ  mit per. Randb. (TN ≡ NT rückwärts)
    SN = (nt+N_s-2)%N_s +1  # ⃗n - ̂ν

    a = -config[MU][nt][ns] - config[mu][nt][NS] + config[MU][NT][ns]
    b = config[MU][nt][SN] - config[mu][nt][SN] - config[MU][NT][SN]

    exp(im*a) + exp(im*b)
end

#=
# Metropolis-Algorithmus:
# Hier werden erst alle Links mit geradem Index, dann alle mit ungeradem Index
# geupdatet (und das zuerst für μ, dann für ν).
raw_configs = [config_1]
for k = 1:m
    test = deepcopy(raw_configs[k])
    # Zunächst werden alle temporalen Links geupdatet (also μ=1)
    # Zuerst die mit geradem Index (für jeden Zeitpunkt)
    for i = 1:N_t
        for j = 1:Int64(0.5*N_s)
            J = 2*j
            proposal = deepcopy(test[1][i][J])
            q = 2*(rand(1)[1] - 0.5)
            proposal += r*q

            a = proposal
            b = test[1][i][J]
            c = staple_dagger(1, i, J, test)
            Delta = -beta * real((exp(im*a)-exp(im*b)) * c)

            P = rand(1)[1]
            if P < exp(-Delta)
                test[1][i][J] = proposal
            end
        end
    end
    # Anschließend die mit ungeradem Index (für jeden Zeitpunkt)
    for i = 1:N_t
        for j = 1:Int64(0.5*N_s)
            J = 2*j-1
            proposal = deepcopy(test[1][i][J])
            q = 2*(rand(1)[1] - 0.5)
            proposal += r*q

            a = proposal
            b = test[1][i][J]
            c = staple_dagger(1, i, J, test)
            Delta = -beta * real((exp(im*a)-exp(im*b)) * c)

            P = rand(1)[1]
            if P < exp(-Delta)
                test[1][i][J] = proposal
            end
        end
    end
    # Nun werden die räumlichen Links geupdatet (μ=2)
    # Zunächst die mit geradem Index (für jeden Raumpunkt)
    for j = 1:N_s
        for i = 1:Int64(0.5*N_t)
            I = 2*i
            proposal = deepcopy(test[2][I][j])
            q = 2*(rand(1)[1] - 0.5)
            proposal += r*q

            a = proposal
            b = test[2][I][j]
            c = staple_dagger(2, I, j, test)
            Delta = -beta * real((exp(im*a)-exp(im*b)) * c)

            P = rand(1)[1]
            if P < exp(-Delta)
                test[2][I][j] = (proposal)
            end
        end
    end
    # Anschließend die mit ungeradem Index (für jeden Raumpunkt)
    for j = 1:N_s
        for i = 1:Int64(0.5*N_t)
            I = 2*i-1
            proposal = deepcopy(test[2][I][j])
            q = 2*(rand(1)[1] - 0.5)
            proposal += r*q

            a = proposal
            b = test[2][I][j]
            c = staple_dagger(2, I, j, test)
            Delta = -beta * real((exp(im*a)-exp(im*b)) * c)

            P = rand(1)[1]
            if P < exp(-Delta)
                test[2][I][j] = (proposal)
            end
        end
    end

    push!(raw_configs, test)
end
=#



# Funktion um Wirkungsdifferenz ΔS mit adj. Staples zu berechnen
function Delta_staple(new_config, config, mu)
    a = 0
    for i = 1:N_t
        for j = 1:N_s
            b = new_config[mu][i][j]
            c = config[mu][i][j]
            if b != c
                d = staple_dagger(mu, i, j, config)
                a += -beta * real((exp(im*b)-exp(im*c)) *d)
            end
        end
    end
    a
end

#=
# Noch ein Versuch eines Metropolis:
raw_configs = [config_1]
for k = 1:m
    test = deepcopy(raw_configs[k])

    for mu = 1:2
        proposal = deepcopy(test)
        for i = 1:N_t
            for j = 1:Int64(0.5*N_s)
                J = 2*j
                q = 2*(rand(1)[1] -0.5)
                proposal[mu][i][J] += r*q
            end
        end
        P = rand(1)[1]
        Delta = Delta_staple(proposal, test, mu)
        if P < exp(-Delta)
            test = proposal
        end
    end

    for mu = 1:2
        proposal = deepcopy(test)
        for i = 1:N_t
            for j = 1:Int64(0.5*N_s)
                J = 2*j-1
                q = 2*(rand(1)[1] -0.5)
                proposal[mu][i][J] += r*q
            end
        end
        P = rand(1)[1]
        Delta = Delta_staple(proposal, test, mu)
        if P < exp(-Delta)
            test = proposal
        end
    end

    push!(raw_configs, test)
end
=#


#=
# Sogar noch ein Versuch eines Metropolis, bringt aber auch acceptence = 100
raw_configs = [config_1]
for k = 1:m
    test = deepcopy(raw_configs[k])

        proposal = deepcopy(test)
        for i = 1:N_t
            for j = 1:Int64(0.5*N_s)
                J = 2*j
                q = 2*(rand(1)[1] -0.5)
                proposal[1][i][J] += r*q
            end
        end
        P = rand(1)[1]
        Delta = Delta_staple(proposal, test, 1)
        if P < exp(-Delta)
            test = proposal
        end

        proposal = deepcopy(test)
        for i = 1:N_t
            for j = 1:Int64(0.5*N_s)
                J = 2*j-1
                q = 2*(rand(1)[1] -0.5)
                proposal[1][i][J] += r*q
            end
        end
        P = rand(1)[1]
        Delta = Delta_staple(proposal, test, 1)
        if P < exp(-Delta)
            test = proposal
        end

        proposal = deepcopy(test)
        for j = 1:N_s
            for i = 1:Int64(0.5*N_t)
                I = 2*i
                q = 2*(rand(1)[1] -0.5)
                proposal[2][I][j] += r*q
            end
        end
        P = rand(1)[1]
        Delta = Delta_staple(proposal, test, 1)
        if P < exp(-Delta)
            test = proposal
        end

        proposal = deepcopy(test)
        for j = 1:N_s
            for i = 1:Int64(0.5*N_t)
                I = 2*i-1
                q = 2*(rand(1)[1] -0.5)
                proposal[2][I][j] += r*q
            end
        end
        P = rand(1)[1]
        Delta = Delta_staple(proposal, test, 1)
        if P < exp(-Delta)
            test = proposal
        end

    push!(raw_configs, test)
end
=#


# Funktion, um die Phase einer Plaquette Pμν( ⃗n) zu berechnen
function P12(nt, ns, config)
    NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
    NS = ns%N_s +1  # nₛ+1 mit periodische Randbed.
    config[1][nt][ns] + config[2][NT][ns] - config[1][nt][NS] - config[2][nt][ns]
end

# Funktion, um die Wirkung einer ganzen Konfig zu bestimmen
function S(config)
    a = 0
    for i = 1:N_t
        for j = 1:N_s
            b = P12(i,j,config)
            a += real(1 - exp(im*b))
        end
    end
    a = beta*a
end




# Metropolis Algorithmus
# Hier wird nicht der Staple verwendet, sondern auf die altmodische Weise
# eine Konfig geupdatet, deren Wirkung gemessen usw. Klappt aber super!
# (aber schon nach einer Art Schachbrettmuster)
raw_configs = [config_1]
for k = 1:m
    test = deepcopy(raw_configs[k])

    for mu = 1:2
        proposal = deepcopy(test)
        for i = 1:N_t
            for j = 1:Int64(0.5*N_s)
                J = 2*j
                proposal[mu][i][J] += r*rand(1)[1]
            end
        end
        P = rand(1)[1]
        Delta = S(proposal) - S(test)
        if P < exp(-Delta)
            test = proposal
        end
    end

    for mu = 1:2
        proposal = deepcopy(test)
        for i = 1:N_t
            for j = 1:Int64(0.5*N_s)
                J = 2*j-1
                for mu = 1:2
                    proposal[mu][i][J] += r*rand(1)[1]
                end
            end
        end
        P = rand(1)[1]
        Delta = S(proposal) - S(test)
        if P < exp(-Delta)
            test = proposal
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

println("", "And the average plaquette is ", mean(Plaquettes_avg), "\n")
