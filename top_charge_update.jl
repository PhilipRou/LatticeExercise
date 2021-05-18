using Plots
using Statistics
#=
Die Phasen von Eichfeld-Konfigurationen Uμ( ⃗n) werden in einem Array
config[μ, nₜ, nₛ] gespeichert, wobei ⃗n = (nₜ, nₛ). Im weiteren soll
μ = 1 für die zeitliche und μ = 2 für die räumliche Dimension stehen.
(es wird auch μ ≡ 1 für temp. und ν ≡ 2 für räuml. synonym verwendet)
Phasen der Eichtransformationen Ω( ⃗n) werden als trans[nₜ, nₛ] gespeichert.
=#
N_t = 32        # Anzahl temporaler Gitterpunkte
N_s = 32    	# Anzahl räumlicher Gitterpunkte
beta = 12.8      # 1/g², wobei g die Eich-Kopplung ist
r = 2*pi*0.05      # Schrittweite im Metropolis (nur in pos. Richtung, s.u.)
nb_sim = 1000   # Anzahl der Update-Durchläufe. Pro Durchlauf werden...
nb_metro = 1    # nb_metro Metropolis-Updates und
nb_top = 1      # nb_top topologische-Updates durchgeführt
Nb_updates = nb_sim*(nb_metro+nb_top)   # Anzahl insg. durchgeführter Updates
cut = Int64(0.6*nb_sim) # Thermalization
N_skip = 20         # Nur jede N_skip-te Konfig wird gemessen
link_accepted = 0   # Misst die Akzeptanz-Rate der Link-Updates im Metropolis
R = r/(2*pi)


# Die erste Feldkonfig. wird per Zufall generiert, aber mit Werten ∈ [0,2π]
config_1 = []
for mu = 1:2
    push!(config_1, [])
    for i = 1:N_t
        a = (2*pi).*rand(N_s)
        push!(config_1[mu], a)
    end
end

# Funktion um einen adjungierten Staple in mu-Richtung an der Stelle
# ⃗n = (nₜ,nₛ) zu bauen
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


# Metropolis einer Konfiguration "config"
# Das Ergebnis ist ein Array, dessen erster Eintrag die geupdatete Konfig
# und dessen zweiter Eintrag die Anzahl der akzeptierten Link-Updates ist
function metro_update(config)
    nb_updates = 0
    test = deepcopy(config)
    for eo = 0:1
        for i = 1:N_t
            for j = 1:Int64(0.5*N_s)
                J = 2*j -eo
                a = test[1][i][J]
                b = staple_dagger(1,i,J,test)
                proposal = deepcopy(a)
                q = 2*(rand(1)[1] -0.5)
                proposal += r*q

                Delta = -beta * real((exp(im*proposal)- exp(im*a)) * b)
                P = rand(1)[1]
                if P < exp(-Delta)
                    test[1][i][J] = proposal
                    nb_updates += 1
                end
            end
        end
        for j = 1:N_s
            for i = 1:Int64(0.5*N_t)
                I = 2*i -eo
                a = test[2][I][j]
                b = staple_dagger(2,I,j,test)
                proposal = deepcopy(a)
                q = 2*(rand(1)[1] -0.5)
                proposal += r*q

                Delta = -beta * real((exp(im*proposal)- exp(im*a)) * b)
                P = rand(1)[1]
                if P < exp(-Delta)
                    test[2][I][j] = proposal
                    nb_updates += 1
                end
            end
        end
    end
    result = [test, nb_updates]
    result
end

# Funktion um die Plaquette einer Konfig bei ⃗n = (nₜ,nₛ) zu bestimmen
function P12(nt, ns, config)
    NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
    NS = ns%N_s +1  # nₛ+1 mit periodische Randbed.
    config[1][nt][ns] + config[2][NT][ns] - config[1][nt][NS] - config[2][nt][ns]
end

# Funktion um Wirkung einer gegebenen Konfig zu messen
function action(config)
    S = 0
    #for mu = 1:2
        for i = 1:N_t
            for j = 1:N_s
                plaq = exp(im*P12(i,j,config))
                S += beta*real(1-plaq)
            end
        end
    #end
    S
end

# Funktion um die topologische Ladung einer Konfig zu messen
function charge(config)
    Q = 0
    for i = 1:N_t
        for j = 1:N_s
            plaq = exp(im*P12(i,j,config))
            Q += (imag(log(plaq))) / (2*pi)
        end
    end
    Q
end

# Funktion um eine n-Instanton-Konfig zu bauen
function instanton(n)
    new_config = [[],[]]
    for i = 1:N_t
        push!(new_config[1], zeros(N_s))
        push!(new_config[2], n*i*2*pi/(N_t*N_s).*ones(N_s))
    end
    for j = 1:N_s
        new_config[1][N_t][j] = -n*j*2*pi/N_s
    end
    new_config
end


# Funktion für ein Update der topologischen Ladung
function top_update(config)
    #new_config = [[],[]]
    if rand() < 0.5
        insta_1 = instanton(1)
    else
        insta_1 = -instanton(1)
    end
    new_config = config.+insta_1
    #=
    for mu = 1:2
        for i = 1:N_t
            push!(new_config[mu], [])
            for j = 1:N_s
                a = insta_1[mu][i][j] + config[mu][i][j]
                push!(new_config[mu][i],a)
            end
        end
    end
    =#

    Delta = action(config) - action(new_config)
    if rand() < exp(Delta)
        return new_config
    else
        return config
    end
end


#                  ##### Die eigentliche Simulation #####
# Zunächst werden die "rohen" Konfigs durch Metropolis- und top. Updates erstellt
raw_configs = [config_1]
for i = 1:nb_sim
    for j = 1:nb_metro
        J = (nb_metro+nb_top)*(i-1) + j
        a = metro_update(raw_configs[J])
        link_accepted += a[2]
        push!(raw_configs, a[1])
    end
    for j = 1:nb_top
        J = (nb_metro+nb_top)*(i-1) + nb_metro + j
        push!(raw_configs, top_update(raw_configs[J]))
    end
end

acceptance = round(link_accepted/(nb_sim*nb_metro*N_t*N_s*2), digits = 3)

# Gemessen wird nur jede N_skip-te Konfig, nach der Thermalization
configs = []
for i = cut:size(raw_configs, 1)
    if i%N_skip == 1
        push!(configs, raw_configs[i])
    end
end
nb_configs = size(configs, 1)

# Plots
actions = []
charges = []
x = []
for i = 1:nb_configs
    push!(actions, action(configs[i]))
    push!(charges, charge(configs[i]))
    push!(x, i)
end


histo1 = histogram(actions, bins = 35,
label = "Wilson Gauge Actions",
title = "Nₛ×Nₜ = $N_s × $N_t,   β = $beta,   nb_sim⋅(nb_metro+nb_top) = $nb_sim⋅($nb_metro+$nb_top),
cut = $cut,   N_skip = $N_skip,   r/2π = $R,   acceptance = $acceptance " )
histo2 = histogram(charges, bins = 130,
label = "Topological Charges")
image_histo = plot(histo1, histo2, layout=(2,1), size=(800,600))
display(image_histo)

plot1 = plot(x, charges)

#savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\
#Lattice Gauge Theory Varnhorst\\topup_32x32")
