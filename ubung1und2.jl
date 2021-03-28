using Plots
using StatsBase
using LsqFit

# Startparameter
N = 20                  # Dimension der Konfigurationen (Anzahl Gitterpunkte)
ω = 0.5                 # Oszillatorfrequenz, eigentlich omega-hut ̂ω
m = 200000               # Anzahl der Konfigurationen (ohne Startkonfiguration)
cut = 100000              # Anzahl der ersten Konfigurationen, die nicht beachtet werden
N_skip = 200            # Anzahl der Konfig., die nicht gemessen werden
r = 0.1                 # Schrittweite
x_0 = 5*(rand(N)-0.5*ones(N))   #Startkonfiguration


# Funktion um Wirkung S einer Konfiguration zu bestimmen
function S(array)
    a = 0
    for i = 1:N-1
        b = i%N +1  # Index für periodische Randbed.
        a = a + (1/ω)*(array[b]-array[i])^2 + ω*array[i]^2
    end
    a
end

# Dieses Array wird später die Konfigurationen beherbergen
# Allerdings wird nicht jede, sondern nur jede N_skip-te gespeichert, s.u.
raw_configs = [x_0]

accepted = 0

# Der eigentliche Metropolis-Algorithmus, der raw_configs mit Konfigurationen füllt
for i = 1:m
    q = 2*(rand(N)-0.5*ones(N))
    x = raw_configs[i] + r*q
    Delta = S(x) - S(raw_configs[i])
    a = rand(1)[1]

    if a < exp(-Delta)
        push!(raw_configs, x)
        accepted = accepted + 1
    else
        push!(raw_configs, raw_configs[i])
    end
end

# Da nur jede N_skip-te Konfiguration gemessen werden soll, wird hier
# jede N_skipte-te Konfiguration in einem Array gespeichert, allerdings erst
# ab der cut-ten Konfig.
configs = []
for i = cut:m+1
    if i%N_skip == 1
        push!(configs, raw_configs[i])
    else
        nothing
    end
end
nb_configs = size(configs, 1)       # number of configs

println("")
println("Für N = $N,  ω = $ω,  m = $m,  r = $r")
println("wurden $accepted von $m Updates akzeptiert, also ca. ", 100*accepted/m, "%")


# Nun zum Plotten:
positions = Real[]

for i = 1:nb_configs
    for j = 1:N
        push!(positions, configs[i][j])
    end
end

H = StatsBase.fit(Histogram, positions; nbins=60)

model(x, p) = p[1]*exp.(-p[2]^(2).*(x).^2)
p0 = [3000., 1.]
xvalues = Array(H.edges[1])
yvalues = H.weights + zeros(length(H.weights))
pop!(xvalues)
fit = LsqFit.curve_fit(model, xvalues, yvalues, p0)
p = fit.param
image = plot(H, label = "Häufigkeiten", fmt = :png)
image = plot!(-2:0.05:2, model(-2:0.05:2, p), label = "Gauß-Fit")
image = plot!(title = "N = $N,  ω = $ω,  m = $m,  cut = $cut,  N_skip = $N_skip,  r = $r")
image = xlabel!("Position x")


###################### Übung 2 ########################

# Funktion um Korrelation ⟨x(0)⋅x(̂τ)⟩ zu berechnen, wobei ̂τ = n⋅Δτ
function correlation(n)
    a = 0
    b = n%N +1
    for i = 1:nb_configs
        a = a + configs[i][1]*configs[i][b]
    end
    a/(m+2-cut)
end

# Funktion für durchschnittliche Korrelation, also Mitteln über alle i mit
# C(τ) = ⟨x(i)⋅x(i+τ)⟩
function bigger_correlation(n)
    a = 0
    for i = 1:nb_configs
        for j = 1:N
            b = (j+n-1)%N +1
            a = a + configs[i][j]*configs[i][b]
        end
    end
    a/(nb_configs*N)
end

bigger_corrs = []
for i = 0:N
    push!(bigger_corrs, bigger_correlation(i))
end


# Funktionen um E_eff^exp(̂τ + 1/2 Δτ) zu berechnen
function E_exp(n)
    log(corrs[n]/corrs[n+1])
end

function E_exp_big(n)
    log(bigger_corrs[n]/bigger_corrs[n+1])
end

# Funktion um E_eff^3pt(̂τ) zu berechnen
# Die Indizes der Arrays sind um +1 verschoben, da Δτ=0 möglich ist,
# was aber Eintrag Nr. [1] entspricht
function E_pt(n)
    b = (n+1)%N +1
    acosh((bigger_corrs[b] + bigger_corrs[n]) / (2*bigger_corrs[n+1]))
end


# Funktion um Autokorrelation von x(̂τ) zu bestimmen, ̂τ = n⋅Δτ
# Wird gebraucht, um Autokorrelationszeit zu bestimmen
function autocorrelation(n)
    a = 0                   # ⟨x(0)⟩
    b = 0                   # ⟨x(̂τ)⟩
    c = n%N +1              # Index für periodische Randbed.
    for i = 1:nb_configs
        a = a + configs[i][1]
        b = b + configs[i][c]
    end
    a = a/(m+2-cut)
    b = b/(m+2-cut)
    bigger_correlation(n) - a*b    # ⟨x(0)⋅x(̂τ)⟩ - ⟨x(0)⟩⋅⟨x(̂τ)⟩
end

function norm_autocorr(n)
    autocorrelation(n)/autocorrelation(0)
end

# Funktion um die Autokorrelationszeit τₐ auszurechnen
function autocorrelation_time()
    a = 0
    for i = 1:N
        b = norm_autocorr(i)
        if b > 0
            a = a + b
        else
            break
        end
    end
    0.5 + a
end

# Funktion zur Fehlerabschätzung σ von correlation(n) durch Jackknife.
# J steht für die Länge der Datenblöcke, die jeweils rausgeschnitten werden.
# Laut Eichhorn, Frech ist J ≈ 2⋅τₐ eine gute Wahl
function Jackknife(n, J)
    R = (nb_configs-J)
    A = floor(Int64, R/J)       #Anzahl der Datenblöcke
    Jack_corrs = Real[]
    mean_corr = bigger_correlation(n)

    # Für jeden Block wird nun die correlation ausgerechnet, wobei zwei
    # Schleifen nötig sind, um einen Datenblock jeweils auszuschneiden
    for i = 1:A     # für jeden Datenblock Nr.i
        a = 0
        for j = 1:i*J           # bis zur (i⋅J)-ten Konfiguration
            for k = 1:N             # Ausrechnen der unnorm. symm. Korr.
                b = (k+n-1)%N +1
                a = a + configs[j][k]*configs[j][b]
            end
        end
        for j = (i+1)*J+1:A*J   # ab der ersten Konfig nach ausgeschn. Block
            for k = 1:N             # Ausrechnen der unnorm. symm. Korr.
                b = (k+n-1)%N +1
                a = a + configs[j][k]*configs[j][b]
            end
        end
        a = a/(R*N)                 # Normieren der noch unnorm. symm. Korr.
        push!(Jack_corrs, a)
    end

    # Hier wird der eigentliche Jackknife-Fehler ausgerechnet, also σ = sqrt(Var)
    b = 0
    for i = 1:A
        b = b + (mean_corr - Jack_corrs[i])^2
    end
    sqrt(b/((A-1)*A))
end

# Funktion um die Fehler der Energien E_exp auszurechnen, basierend auf
# Gaußscher Fehlerfortpflanzung, nutzt Δcorrelation(n) = Jackknife(n, J).
function error_E_exp(n, J)
    j = n%N +1
    A = bigger_corrs[j]
    B = bigger_corrs[j+1]
    a = Jackknife(n, J)
    b = Jackknife(n+1, J)
    sqrt((a/A)^2 + (b/B)^2)
end

#function E_pt(n)
#    b = (n+1)%N +1
#    acosh((bigger_corrs[b] + bigger_corrs[n]) / (2*bigger_corrs[n+1]))
#end

#=
function error_E_pt(n,J)
    j = (n+1)%N +1
    A = bigger_corrs[j]
    B = bigger_corrs[n]
    C = bigger_corrs[n+1]
    a = Jackknife(j, J)
    b = Jackknife(n, J)
    c = Jackknife(n+1, J)

    sqrt(a^2/(A^2-4*C^2) + b^2/(B^2 - 4*C^2) + (c*(A+B/2*C))^2 / ((A+B)^2-4*C^2) )
end
println(error_E_pt(1,8))
=#

# Da laut Eichhorn, Frech J ≈ 2⋅τₐ eine gute Wahl für Jackknife ist, wird
# dies im folgenden berechnet und für die Fehler der Energien verwendet
autotime = autocorrelation_time()
T = floor(Int64, 2*autotime)


errors_exp = Real[]
for i = 1:Int64(N/2 +1)
    push!(errors_exp, error_E_exp(i, T))
end


# Nun zum Plotten, schon wieder :D
# Erst werden Arrays für die Propagatoren, dann für die eff. Energien gemacht

y_corrs= bigger_corrs
x_corrs = Real[]
for i = 0:N
    push!(x_corrs, i)
end

y_energy = Real[]   # Hier werden die eff. Energien gespeichert
x_energy = Real[]
for i = 1:Int64(N/2 +1)
    a = bigger_corrs[i]/bigger_corrs[i+1]

    # folgendes verhindert Auftreten neg. Argumente im log()
    if a > 0
        push!(y_energy, E_exp_big(i))
        push!(x_energy, i-1)
    else
        nothing
    end
end

y_energy_pt = Real[]
for i = 1:Int64(N/2 +1)
    b = (i+1)%N +1
    a = (bigger_corrs[b] + bigger_corrs[i]) / (2*bigger_corrs[i+1])

    # folgendes verhindert Auftreten zu kleinen Argumenten im acosh()
    if a > 1
        push!(y_energy_pt, E_pt(i))
    else
        nothing
    end
end


image2 = plot(x_corrs, y_corrs, label = "symm Propagatoren C(τ)",
              seriestype = :scatter)
image2 = plot!(x_energy, y_energy, yerror = errors_exp,
            label = "Effektive Energie E_exp_big", seriestype = :scatter)
image2 = plot!(x_energy, y_energy_pt,
            label = "Effektive Energie E_pt", seriestype = :scatter)
image2 = xlabel!("Zeit-Differenz τ")

big_image = plot(image, image2, layout = (2,1), size=(750,600), legend = :topright)
display(big_image)


#savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\
#Lattice Gauge Theory Varnhorst\\ubung2_test_plot3")
