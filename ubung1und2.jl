using Plots
using StatsBase
using LsqFit

# Startparameter
N = 20                  # Dimension der Konfigurationen
ω = 0.5                 # Oszillatorfrequenz, omega-Hut
m = 10000               # Anzahl der Konfigurationen (ohne Startkonfiguration)
cut = 1000              # Anzahl der ersten Konfigurationen, die nicht beachtet werden
r = 0.1                 # Schrittweite
x_0 = 5*(rand(N)-0.5*ones(N))   #Startkonfiguration


# Funktion um Wirkung S einer Konfiguration zu bestimmen
function S(array)
    a = 0
    for i = 1:N-1
        a = a + (1/ω)*(array[i+1]-array[i])^2 + ω*array[i]^2
    end
    a + (1/ω)*(array[1]-array[N])^2 + ω*array[N]^2
end

# Dieses Array wird später die Konfigurationen beherbergen
configs = []
push!(configs, x_0)

accepted = 0

# Der eigentliche Metropolis-Algorithmus, der configs mit Konfigurationen füllt
for i = 1:m
    q = 2*(rand(N)-0.5*ones(N))
    x = configs[i] + r*q
    Delta = S(x) - S(configs[i])
    #println(Delta)
    a = rand(1)

    if a[1] < exp(-Delta)
        push!(configs, x)
        accepted = accepted + 1
    else
        push!(configs, configs[i])
    end

end


println("")
println("Es wurden $accepted von $m Updates akzeptiert, also ca. ", 100*accepted/m, "%")


# Nun zum Plotten:
positions = Real[]

for i = cut:m+1
    for j = 1:N
        push!(positions, configs[i][j])
    end
end


H = StatsBase.fit(Histogram, positions; closed=:left, nbins=60)

model(x, p) = p[1]*exp.(-p[2]^(2).*(x).^2)
p0 = [3000., 1.]
xvalues = Array(H.edges[1])
yvalues = H.weights + zeros(length(H.weights))
pop!(xvalues)
fit = LsqFit.curve_fit(model, xvalues, yvalues, p0)
p = fit.param
image = plot(H, label = "Häufigkeiten", fmt = :png)
image = plot!(-2:0.05:2, model(-2:0.05:2, p), label = "Gauß-Fit")
image = plot!(title = "N = $N,  ω = $ω,  m = $m,  r = $r")
image = xlabel!("Position x")
display(image)
#savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\Lattice Gauge Theory Varnhorst\\ubung1_plot")


###################### Übung 2 ########################

# Funktion um Korrelation ⟨x(0)x(̂τ)⟩ zu berechnen, wobei ̂τ = n⋅Δτ
function correlation(n)
    a = 0
    for i = cut:m+1
        a = a + configs[i][1]*configs[i][n]
    end
    a/(m+2-cut)
end

# Funktion um Autokorrelation von x(̂τ) zu bestimmen, ̂τ = n⋅Δτ
function autocorrelation(n)
    a = 0
    b = 0
    for i = cut:m+1
        a = a + configs[i][1]
        b = b + configs[i][n]
    end
    a = a/(m+2-cut)
    b = b/(m+2-cut)
    correlation(n) - a*b
end

# Funktion um E_eff^exp(̂τ + 1/2 Δτ) zu berechnen,wobei die Fallunterscheidung
# wegen der periodischen Randbedingungen gemacht wird
function E_exp(n)
    a = 0
    # Wegen periodischer Randbedingungen
    if n == N
        a = correlation(n)/correlation(1)
    else
        a = correlation(n)/correlation(n+1)
    end
    log(a)
end

function E_pt(n)
    a = 0
    # Wegen periodischer Randbedingungen
    if n == N
        a = (correlation(1) + correlation(n-1)) / (2*correlation(n))
    elseif n == 1
        a = (correlation(n+1) + correlation(N)) / (2*correlation(n))
    else
        a = (correlation(n+1) + correlation(n-1)) / (2*correlation(n))
    end
    acosh(a)
end

#=
Bei der Berechnung der Energien treten Fehler auf, da bei Startkonfigurationen
mit negativen Einträgen auch negative Korrelationen heraus kommen können,
sodass das Argument im log() negativ wird (bzw. das Argument
von achosh() kleiner als 1
=#

# Funktion zur Fehlerabschätzung von correlation(n) durch Jackknife.
# J steht für die Breite der Datenblöcke, die jeweils rausgeschnitten werden,
function Jackknife(n, J)
    R = m+2-cut-J
    A = floor(Int64, R/J)       #Anzahl der Datenblöcke
    Jack_corrs = Real[]
    mean_corr = correlation(n)

    # Für jeden Block wird nun die correlation ausgerechnet, wobei zwei
    # Schleifen nötig sind, um einen Datenblock jeweils auszuschneiden
    for i = 1:A
        a = 0
        for j = cut:cut+(i-1)*J-1
            a = a + configs[j][1]*configs[j][n]
        end
        for k = cut+i*J:cut+A*J
            a = a + configs[k][1]*configs[k][n]
        end
        a = a/R
        push!(Jack_corrs, a)
    end

    # Hier wird nun der eigentliche Jackknife-Fehler ausgerechnet
    b = 0
    for i = 1:A
        b = b + (mean_corr - Jack_corrs[i])^2
    end

    b/(A*(A-1))
end

# Funktion um die Fehler der Energien E_exp auszurechnen, basierend auf
# Gaußscher Fehlerfortpflanzung, nutzt Δcorrelation(n) = Jackknife(n, J)

function error_E_exp(n, J)
    A = correlation(n)
    B = correlation(n+1)
    a = Jackknife(n, J)
    b = Jackknife(n+1, J)
    sqrt((a/A)^2 + (b/B)^2)
end


# Listen, die man sich evtl ausgeben lassen möchte
actions = Real[]
correlations = Real[]
autocorrelations = Real[]
exp_energies = Real[]
pt_energies = Real[]

for i = 1:N
    push!(actions, S(configs[i]))
    push!(correlations, correlation(i))
    push!(autocorrelations, autocorrelation(i))
    #push!(exp_energies, E_exp(i))              #klappt nicht, s.o.
    #push!(pt_energies, E_pt(i))                #klappt nicht, s.o.
end
