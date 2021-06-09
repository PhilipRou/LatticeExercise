using Plots
using Statistics
#=
Die Phasen von Eichfeld-Konfigurationen Uμ( ⃗n) werden in einem Array
config[μ, nₜ, nₛ] gespeichert, wobei ⃗n = (nₜ, nₛ). Im weiteren soll
μ = 1 für die zeitliche und μ = 2 für die räumliche Dimension stehen.
(es wird auch μ ≡ 1 für temp. und ν ≡ 2 für räuml. synonym verwendet)
Phasen der Eichtransformationen Ω( ⃗n) werden als trans[nₜ, nₛ] gespeichert.
=#
N_t = 20        # Anzahl temporaler Gitterpunkte
N_x = 20    	# Anzahl räumlicher Gitterpunkte
β = (N_x*N_t)/80      # 1/a²g², wobei g die Eich-Kopplung ist
m = 25000          # Updates im Metropolis
r = 2*pi*0.08      # Schrittweite im Metropolis (nur in pos. Richtung, s.u.)
cut = Int64(0.6*m)  # Thermalization
N_skip = 20         # Nur jede N_skip-te Konfig wird gemessen


# Die erste Feldkonfig. wird per Zufall generiert, aber mit Werten ∈ [0,2π]
config_1 = []
for mu = 1:2
    push!(config_1, [])
    for i = 1:N_t
        a = (2*pi).*rand(N_x)
        push!(config_1[mu], a)
    end
end

# Funktion, um den adjungierten Staple in Richtung mu bei ⃗n = (nₜ,nₓ) einer
# gegebenen Konfig zu berechnen
function staple_dagger(mu, nt, nx, config)
    a = 0
    b = 0

    NT = nt%N_t +1
    NX = nx%N_x +1
    TN = (nt + N_t -2)%N_t +1
    XN = (nx + N_x -2)%N_x +1

    if mu == 1
        a = -config[2][nt][nx] - config[1][nt][NX] + config[2][NT][nx]
        b = config[2][nt][XN] - config[1][nt][XN] - config[2][NT][XN]
    elseif mu == 2
        a = -config[1][nt][nx] - config[2][NT][nx] + config[1][nt][NX]
        b = config[1][TN][nx] - config[2][TN][nx] - config[1][TN][NX]
    end

    exp(im*a) + exp(im*b)
end


# Der Metropolis-Algorithmus
# Zunächst (eo=0) wird zu jeder Zeit jeder örtlich geradzahlige Link in
# zeitliche Richtung geupdatet, danach an jedem Ort jeder zeitlich geradzahlige
# Link in Ortsrichtung. Anschließend (eo=1) das gleiche nochmal mit dem jeweils
# ungeradzahligen Link.
link_accepted = 0
raw_configs = [config_1]
for k = 1:m
    test = deepcopy(raw_configs[k])

    for eo = 0:1  # even-odd
        for i = 1:N_t
            for j = 1:Int64(0.5*N_x)
                J = 2*j -eo
                a = test[1][i][J]
                proposal = deepcopy(a)
                q = 2*(rand() -0.5)
                proposal += r*q

                b = staple_dagger(1,i,J,test)
                Δ = -β * real((exp(im*proposal)- exp(im*a)) * b)
                if rand() < exp(-Δ)
                    test[1][i][J] = proposal
                    link_accepted += 1
                end
            end
        end
        for j = 1:N_x
            for i = 1:Int64(0.5*N_t)
                I = 2*i -eo
                a = test[2][I][j]
                proposal = deepcopy(a)
                q = 2*(rand() -0.5)
                proposal += r*q

                b = staple_dagger(2,I,j,test)
                Δ = -β * real((exp(im*proposal)- exp(im*a)) * b)
                if rand() < exp(-Δ)
                    test[2][I][j] = proposal
                    link_accepted += 1
                end
            end
        end
    end

    push!(raw_configs, test)
end

println("", "Es wurden $link_accepted von ", m*2*N_t*N_x, " Link-Updates akzeptiert, also ",
100*(link_accepted/(m*2*N_t*N_x)), "%")


# Funktion, um die PHASE einer Plaquette Pμν( ⃗n) zu berechnen
function P12(nt, nx, config)
    NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
    NX = nx%N_x +1  # nₛ+1 mit periodische Randbed.
    config[1][nt][nx] + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
end


# Es wird nur jede N_skip-te Konfig gemessen (Autokorrelation)
configs = []
for k = cut:m+1
    if k%N_skip == 0
        push!(configs, raw_configs[k])
    end
end
nb_configs = size(configs, 1)

# Funktion um eine n-Instanton-Konfig zu generieren
function instanton(n)
    new_config = [[],[]]
    for i = 1:N_t
        push!(new_config[1], zeros(N_x))
        push!(new_config[2], n*i*2*pi/(N_t*N_x).*ones(N_x))
    end
    for j = 1:N_x
        new_config[1][N_t][j] = -n*j*2*pi/N_x
    end
    new_config
end

# Funktion um die ganzzahlige topologische Ladung einer Konfig zu messen
function int_charge(config)
    Q = 0
    for i = 1:N_t
        for j = 1:N_x
            plaq = exp(im*P12(i,j,config))
            Q += (imag(log(plaq))) / (2*pi)
        end
    end
    Q
end

# Funktion um die "kontinuierliche" topologische Ladung einer Konfig zu messen
function cont_charge(config)
    Q = 0
    for i = 1:N_t
        for j = 1:N_x
            plaq = exp(im*P12(i,j,config))
            Q += (imag(plaq)) / (2*pi)
        end
    end
    Q
end

# Für die Plots: Listen erstellen mit den ganzzahligen und kontinuierlichen
# Ladungen
int_charges = []
cont_charges = []
for k = 1:nb_configs
    push!(int_charges, int_charge(configs[k]))
    push!(cont_charges, cont_charge(configs[k]))
end
# Die ganzzahligen werden auf ganze Zahlen gerundet (begrenzte
# Maschinenpräzision führt zu Abweichungen der Ordnung e-10 und kleiner)
int_charges = round.(Int, int_charges)
# Ganzzahlige Ladungen werden für das Histogram so verschoben, dass die
# Säulen auf den Zahlen liegen (s.u.)
int_charges = int_charges.-0.25


# Das 2D-Histogram zum Vergleich der diskreten und kont. Ladungen
histo2d = histogram2d(
int_charges, cont_charges, bins = 100, c = :reds, background_color = :grey)

histo2d = plot!(
size = (750,600),
title = "Different top. charges
# $N_t×$N_x, β = $β, $m-$cut sweeps, N_skip = $N_skip",
xticks = yticks = -6:1:6,
xlabel = "Discrete Q",
ylabel = "Continuous Q")

display(histo2d)


# Histogram der diskreten Ladungsverteilung, muss bei jedem Run angepasst werden
histo_int = histogram(
int_charges, bins = -4.25:0.5:5, c = :salmon, background_color = :grey,
bar_edges = :false)

histo_int = plot!(
size = (750,600), title = "Discrete top. charges
$N_t×$N_x, β = $β, $m-$cut sweeps, N_skip = $N_skip",
xticks = -7:1:7,
xlabel = "Discrete Q",
legend = :false)

display(histo_int)


# Histogram der kontinuierlichen Ladungsverteilung, muss bei jedem Run angepasst
histo_cont = histogram(
cont_charges, bins = -3:0.06:3, c = :salmon, background_color = :grey,
bar_edges = :false)

histo_cont = plot!(
size = (750,600), title = "Continuous top. charges
$N_t×$N_x, β = $β, $m-$cut sweeps, N_skip = $N_skip",
xticks = -6:1:6,
xlabel = "Continuous Q",
legend = :false)

display(histo_cont)


# Zeitserie der Ladungen (kont. oder diskret muss jedes Mal angepasst werden)
# Der Cut ist mit inbegriffen, die Thermalization also nicht durchgeführt

# nb_raw = size(raw_configs, 1)
# x = []
# raw_cont_charges = []
# raw_int_charges = []
# for k = 1:nb_raw
#     push!(x, k)
#     push!(raw_cont_charges, cont_charge(raw_configs[k]))
#     push!(raw_int_charges, int_charge(raw_configs[k]))
# end
#
# image = plot(
# x, raw_cont_charges, seriescolor = :salmon,
# seriestype = :scatter
# )
#
# image = plot!(
# size = (750, 600), title = "Time series of continuous charges (no thermalization)
# $N_t×$N_x, β = $β, $m-$cut sweeps, N_skip = $N_skip",
# ylabel = "cont. charges",
# # yticks = -8:1:2,
# legend = :false,
# background_color = :grey
# )


#savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\compare_20x20_larger")
