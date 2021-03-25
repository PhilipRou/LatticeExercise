using Plots
using StatsBase
using LsqFit

# Startparameter
N = 20                  # Dimension der Konfigurationen
ω = 0.1                 # Oszillatorfrequenz, eigentlich omega-hut ̂ω
m = 80000               # Anzahl der Konfigurationen (ohne Startkonfiguration)
cut = 40000              # Anzahl der ersten Konfigurationen, die nicht beachtet werden
r = 0.08                 # Schrittweite
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
configs = [x_0]

accepted = 0

# Der eigentliche Metropolis-Algorithmus, der configs mit Konfigurationen füllt
for i = 1:m
    q = 2*(rand(N)-0.5*ones(N))
    x = configs[i] + r*q
    Delta = S(x) - S(configs[i])
    a = rand(1)[1]

    if a < exp(-Delta)
        push!(configs, x)
        accepted = accepted + 1
    else
        push!(configs, configs[i])
    end

end


println("")
println("Für N = $N,  ω = $ω,  m = $m,  r = $r")
println("wurden $accepted von $m Updates akzeptiert, also ca. ", 100*accepted/m, "%")


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


###################### Übung 2 ########################

# Funktion um Korrelation ⟨x(0)⋅x(̂τ)⟩ zu berechnen, wobei ̂τ = n⋅Δτ
function correlation(n)
    a = 0
    b = n%N +1
    for i = cut:m+1
        a = a + configs[i][1]*configs[i][b]
    end
    a/(m+2-cut)
end

# Funktion für durchschnittliche Korrelation, also Mitteln über alle i mit
# C(τ) = ⟨x(i)⋅x(i+τ)⟩
function bigger_correlation(n)
    a = 0
    for i = cut:m+1
        for j = 1:N
            b = (j+n-1)%N +1
            a = a + configs[i][j]*configs[i][b]
        end
    end
    a/((m+2-cut)*N)
end

corrs = []
bigger_corrs = []
for i = 0:N
    push!(corrs, correlation(i+1))
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
function E_pt(n)
    acosh((bigger_corrs[n+1] + bigger_corrs[n-1]) / bigger_corrs[n])
end


# Funktion um Autokorrelation von x(̂τ) zu bestimmen, ̂τ = n⋅Δτ
# Wird gebraucht, um Autokorrelationszeit zu bestimmen
function autocorrelation(n)
    a = 0                   # ⟨x(0)⟩
    b = 0                   # ⟨x(̂τ)⟩
    c = n%N +1              # Index für periodische Randbed.
    for i = cut:m+1
        a = a + configs[i][1]
        b = b + configs[i][c]
    end
    a = a/(m+2-cut)
    b = b/(m+2-cut)
    correlation(n) - a*b    # ⟨x(0)⋅x(̂τ)⟩ - ⟨x(0)⟩⋅⟨x(̂τ)⟩
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
# J steht für die Breite der Datenblöcke, die jeweils rausgeschnitten werden.
# Laut Eichhorn, Frech ist J ≈ 2⋅τₐ eine gute Wahl
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
            b = n%N +1
            a = a + configs[j][1]*configs[j][b]
        end
        for k = cut+i*J:cut+A*J
            b = n%N +1
            a = a + configs[k][1]*configs[k][b]
        end
        a = a/R
        push!(Jack_corrs, a)
    end

    # Hier wird der eigentliche Jackknife-Fehler ausgerechnet, also σ = sqrt(Var)
    b = 0
    for i = 1:A
        b = b + (mean_corr - Jack_corrs[i])^2
    end
    sqrt((b*(A-1))/A)
end

# Funktion um die Fehler der Energien E_exp auszurechnen, basierend auf
# Gaußscher Fehlerfortpflanzung, nutzt Δcorrelation(n) = Jackknife(n, J).
# Kann in Zukunft noch für bigger_correlation gemacht werden
function error_E_exp(n, J)
    A = correlation(n)
    B = correlation(n+1)
    a = Jackknife(n, J)
    b = Jackknife(n+1, J)
    sqrt((a/A)^2 + (b/B)^2)
end


yplot_big = bigger_corrs
yplot = corrs
xplot = Real[]
for i = 0:N
    push!(xplot, i)
end

y_energy_big = Real[]
y_energy = Real[]
#y_energy_pt = Real[]                           # klappt nicht ❌
x_energy = Real[]

for i = 1:N/2
    push!(y_energy, E_exp(Int64(i)))
    push!(y_energy_big, E_exp_big(Int64(i)))
    #push!(y_energy_pt, E_pt(Int64(i)))         # klappt nicht ❌
    push!(x_energy, i-1)
end

image2 = plot(xplot, yplot_big, label = "symm Propagatoren C(τ)")
image2 = plot!(x_energy, y_energy_big, label = "Effektive Energie E_exp_big")
image2 = plot!(xplot, yplot, label = "Propagatoren C(τ)")
image2 = plot!(x_energy, y_energy, label = "Effektive Energie E_exp")
#image2 = plot!(x_energy, y_energy_pt, label = "Effektive Energie E_3pt")
#image2 = plot!(title = "N = $N,  ω = $ω,  m = $m,  r = $r", legend = :right)
image2 = xlabel!("Zeit-Differenz τ")

big_image = plot(image, image2, layout = (2,1), size=(750,600), legend = :bottomright)
display(big_image)
#savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\
#Lattice Gauge Theory Varnhorst\\ubung2_test_plot2")
