using Plots
using StatsBase
using LsqFit

# Simulationsparameter
M = 0.4     # Masse des freien Teilchens
d_t = 1     # Räumliche Dimensionen
d_s = 1     # Temporale Dimensionen
N_t = 8     # Anzahl räumlicher Gitterpunkte
N_s = 8     # Anzahl temporaler Gitterpunkte
m = 100000     # Anzahl Wiederholungen in Metropolis, also ∃ m+1 Konfigurationen
r = 0.1     # Schrittweite in Metropolis
cut = 60000   # Anzahl der ersten Konfig., die verworfen werden (Thermalisation)
N_skip = 100  # Es wird nur jede N_skip-te Konfig. gemessen


# Startkonfiguration generieren
# Der erste Index steht für temporale, der zweite für räumliche Koordinaten
x_0 = []
for i = 1:N_t
    q = 5*(rand(N_s) - 0.5*ones(N_s))
    push!(x_0, q)
end
x_0 = 4*x_0

# Funktion um Wirkung einer Konfiguration auszurechnen (vgl. Gl(56) im Skript)
function S(x)
    a = 0
    for i = 1:N_t
        I = (i)%N_t +1      # period. Randbed., s.u.

        for j in 1:N_s
            J = (j)%N_s +1  # period. Randbed., nämlich:           ↓            ↓
            a = a + (0.5*M^2 + d_s+d_t) * x[i][j]^2 - x[i][j] * (x[I][j] + x[i][J])
        end
    end
    a
end

#=
# Der Metropolis Algorithmus
raw_configs = [x_0]
accepted = 0
accepted = 0
for i = 1:m
    q = []
    for i = 1:N_t
        b = 2*(rand(N_s) - 0.5*ones(N_s))
        push!(q, b)
    end

    proposal = raw_configs[i] + r*q
    Delta = S(proposal) - S(raw_configs[i])
    P = rand(1)[1]
    if P < exp(-Delta)
        push!(raw_configs, proposal)
        accepted = accepted + 1
    else
        push!(raw_configs, raw_configs[i])
    end
end
=#

#=
# Zeilen-Metropolis:
# Variant des Metropolis: anstatt alle 8×8 Einträge zu ändern, wird
# jede ZEILE einzeln verändert, wobei diese Veränderung mit der gleichen
# W-keit wie oben akzeptiert wird. Wenn dies für alle N_t Zeilen einer Konfig
# gemacht wurde, wird die somit neue entstandene Konfig in raw_configs gepusht

raw_configs = [x_0]
accepted = 0
for i = 1:m
    test = copy(raw_configs[i])

    for j = 1:N_t
        q = 2*(rand(N_s)-0.5*ones(N_s))
        proposal = copy(test)
        proposal[j] = proposal[j] + r*q
        Delta = S(proposal) - S(test)
        P = rand(1)[1]

        if P < exp(-Delta)
            test = copy(proposal)
        end
    end

    push!(raw_configs, test)

    if test != raw_configs[i]
        accepted = accepted + 1
    end
end
=#

#=
# Einzelner-Eintrag-Metropolis:
# Variante des Metropolis: anstatt alle 8×8 Einträge zu ändern, wird
# jeder EINTRAG einzeln verändert, wobei diese Veränderung mit der gleichen
# W-keit wie oben akzeptiert wird. Wenn dies für alle Einträge einer Konfig
# gemacht wurde, wird die somit neue entstandene Konfig in raw_configs gepusht
raw_configs = [x_0]
accepted = 0
for i = 1:m
    test = copy(raw_configs[i])

    for j = 1:N_t
        for k = 1:N_s
            q = 2*(rand(1)[1]-0.5)
            proposal = deepcopy(test)
            b = test[j][k]
            proposal[j][k] = proposal[j][k] + r*q
            a = proposal[j][k]

            J = j%N_t +1            # ̂τ+Δτ, jeweils mit period. Randbed.
            K = k%N_s +1            # ⃗x+̂μ
            Y = (N_t+j-2)%N_t +1    # ̂τ-Δτ
            G = (N_s+j-2)%N_s +1    # ⃗x-̂μ

            Delta = (0.5*M^2+d_s+d_t)*(a^2-b^2) -
            (a-b)*(test[j][K] + test[J][k] + test[j][G] + test[Y][k])

            P = rand(1)[1]
            if P < exp(-Delta)
                test = copy(proposal)
            end
        end
    end

    push!(raw_configs, test)
    if test != raw_configs[i]
        accepted = accepted + 1
    end
end
=#


# Schachbrett-Metropolis:
# Man kann sich eine Konfiguration wie ein Schachbrett mit abwechselnd schwarzen
# und weißen Kacheln vorstellen. In dieser Variante werden erst die "weißen"
# Einträge mit einer Zufallszahl addiert (weiß zieht zuerst) und die Wirkung
# dieser neuen Konfiguration gemessen und dem Metropolis-Schritt unterzogen.
# Anschließend wird das Gleiche für die "schwarzen" Einträge wiederholt.
# Erst dann ist "ein ganzes Update" durchgeführt worden, dies wird unten in den
# raw_configs berücksichtigt
raw_raw_configs = [x_0]
for i = 1:2*m

    if i%2 == 1 # Hier wird "weiß" geupdatet, man braucht nur ½NₜNₛ Zuf-Zahlen
        test = deepcopy(raw_raw_configs[i])
        proposal = deepcopy(test)
        D = Int64(0.5*N_t*N_s)
        Q = 2*(rand(D)-0.5*ones(D))

        #=
        Die Idee ist: ungerade Zeilen bekommen Zuf-Zahlen mit ungeraden Indizes,
        die erste Zeile bekommt Q[1] bis Q[Nₛ-1],
        die dritte Zeile kriegt Q[Nₛ+1] bis Q[Nₛ+Nₛ-1], usw...
        Und umgekehrt: gerade Zeilen bekommen Zuf-Zahlen mit geraden Indizes,
        die zweite Zeile bekommt Q[2] bis Q[Nₛ],
        die vierte Zeile bekommt Q[Nₛ+2] bis Q[Nₛ+Nₛ], usw...
        =#
        for j = 1:N_t
            if j%2 == 1
                E = Int64(0.5*N_s -1) # Annahme: N_s ist gerade!
                for k = 0:E
                    F = Int64(1 + 2*k)     # F = 1, 3, 5, ..., N_s-1
                    G = Int64(0.5*N_s*(j-1) + F)
                    proposal[j][F] = proposal[j][F] + r*Q[G]
                end
            elseif j%2 == 0
                E = Int64(0.5*N_s)
                for k = 1:E
                    F = Int64(2*k)         # F = 2, 4, 6, ..., N_s
                    G = Int64(0.5*N_s*(j-2) + F)
                    proposal[j][F] = proposal[j][F] + r*Q[G]
                end
            end
        end

        Delta = S(proposal) - S(test)
        P = rand(1)[1]
        if P < exp(-Delta)
            push!(raw_raw_configs, proposal)
        else
            push!(raw_raw_configs, test)
        end
    end


    if i%2 == 0 # Hier wird "schwarz" geupdatet
        test = deepcopy(raw_raw_configs[i])
        proposal = deepcopy(test)
        D = Int64(0.5*N_t*N_s)
        Q = 2*(rand(D)-0.5*ones(D))

        # Nun bekommen die ungeraden Zeilen die geraden Indizes und umgekehrt
        for j = 1:N_t
            if j%2 == 1
                E = Int64(0.5*N_s)
                for k = 1:E
                    F = Int64(2*k)         # F = 2, 4, 6, ..., N_s
                    G = Int64(0.5*N_s*(j-1) + F)
                    proposal[j][F] = proposal[j][F] + r*Q[G]
                end
            elseif j%2 == 0
                E = Int64(0.5*N_s -1) # Annahme: N_s ist gerade!
                for k = 0:E
                    F = Int64(1 + 2*k)     # F = 1, 3, 5, ..., N_s-1
                    G = Int64(0.5*N_s*(j-2) + F)
                    proposal[j][F] = proposal[j][F] + r*Q[G]
                end
            end
        end

        Delta = S(proposal) - S(test)
        P = rand(1)[1]
        if P < exp(-Delta)
            push!(raw_raw_configs, proposal)
        else
            push!(raw_raw_configs, test)
        end
    end
end

# Wir übernehmen nur raw_raw-Konfigurationen, bei denen schwarz und weiß
# gleich oft dem Metropolis-Schritt unterzogen worden sind:
raw_configs = []
for i = 1:2*m+1
    if i%2 == 1
        push!(raw_configs, raw_raw_configs[i])
    end
end
accepted = 0
for i = 1:m
    if raw_configs[i+1] != raw_configs[i]
        accepted = accepted + 1
    end
end

println("Es wurden $accepted von $m Updates akzepiert, also ca. ", 100*accepted/m, "%" )


# Es wird nur jede N_skip-te Konfig zum Messen verwendet:
configs = []
for i = cut:m+1
    if i%N_skip == 1
        push!(configs, raw_configs[i])
    end
end
nb_configs = size(configs, 1)

# Funktion für diskrete Fourier Transformation einer l-ten Konfiguration x[l],
# wobei n_p ∈ {0,1,...,N_s/2} und τ (tau) die temporale Koordinate
function discrete_fourier(l, τ, n_p)
    a = 0
    p = n_p*(2*pi/N_s)
    for i = 1:N_s
        a = a + exp(-im*p*i) * configs[l][τ][i]
    end
    a
end

# Nun zu den symmetrisierten Korrelatoren Cₚ = ⟨𝚽*(p,i)⋅𝚽(p,i+̂τ)⟩,
# wobei ̂τ = n⋅Δτ, n_p ∈ {0,1,...,N_s/2}
function correlation(n, n_p)
    a = 0
    for i = 1:nb_configs
        for j = 1:N_t
            J = (j+n-1)%N_t +1
            a = a + conj(discrete_fourier(i,j,n_p)) * discrete_fourier(i,J,n_p)
            a = a + conj(discrete_fourier(i,J,n_p)) * discrete_fourier(i,j,n_p)
        end
    end
    real(a/(2*nb_configs*N_t))
end



corrs = [[],[],[],[],[]]
for j = 0:4
    for i = 0:N_t
        push!(corrs[j+1], correlation(i, j))
    end
end
#corrs = Symm. Korrelatoren für [p = 0 | p = ¼π | p = ½π | p = ¾π | p = π]


function E_exp(n, n_p)
    I = Int64(n_p +1)           # [n_p +1] weil Arrays bei [1] beginnen
    exp_corrs = corrs[I]
    a = exp_corrs[n]/exp_corrs[n+1]
    if a > 0
        log(a)
    end
end


function exp_bootstrap(n, n_p, u)
    mean_energy = E_exp(n, n_p)
    bootstrap_corrs_1 = []
    bootstrap_corrs_2 = []
    for i = 1:u
        a = 0
        Q = rand(1:nb_configs, nb_configs)

        for i = 1:nb_configs
            for j = 1:N_t
                J = (j+n-1)%N_t +1
                a = a + conj(discrete_fourier(Q[i],j,n_p)) * discrete_fourier(Q[i],J,n_p)
                a = a + conj(discrete_fourier(Q[i],J,n_p)) * discrete_fourier(Q[i],j,n_p)
            end
        end
        a = real(a/(2*nb_configs*N_t))
        push!(bootstrap_corrs_1, a)

        b = 0
        for i = 1:nb_configs
            for j = 1:N_t
                J = (j+n)%N_t +1
                b = b + conj(discrete_fourier(Q[i],j,n_p)) * discrete_fourier(Q[i],J,n_p)
                b = b + conj(discrete_fourier(Q[i],J,n_p)) * discrete_fourier(Q[i],j,n_p)
            end
        end
        b = real(b/(2*nb_configs*N_t))
        push!(bootstrap_corrs_2, b)
    end

    c = 0
    for i = 1:u
        c = c + (mean_energy - log(bootstrap_corrs_1[i]/bootstrap_corrs_2[i]))^2
    end

    sqrt(c/u)
end


function E_pt(n, n_p)
    I = n_p +1              # [n_p +1] weil Arrays bei [1] beginnen
    pt_corrs = corrs[I]

    a = (n+1+N_t)%N_t +1    # τ+Δτ
    b = (n-1+N_t)%N_t +1    # τ-Δτ
    c = (n+N_t)%N_t   +1    # τ

    A = (pt_corrs[a] + pt_corrs[b]) / (2*pt_corrs[c])
    if A >= 1
        acosh(A)
    end
end


function autocorrelation(n, n_p)
    a = 0
    aa = 0
    b = 0
    bb = 0
    for i = 1:nb_configs
        for j = 1:N_t
            J = (j+n-1)%N_t +1
            a = a + conj(discrete_fourier(i,j,n_p))
            aa = aa + discrete_fourier(i,j,n_p)
            b = b + conj(discrete_fourier(i,J,n_p))
            bb = bb + discrete_fourier(i,J,n_p)
        end
    end
    a = real(a/(2*nb_configs*N_t))
    aa = real(aa/(2*nb_configs*N_t))
    b = real(b/(2*nb_configs*N_t))
    bb = real(bb/(2*nb_configs*N_t))

    correlation(n,n_p) - a*bb -aa*b
end

function autocorrelation_time(n_p)
    a = 0
    for i = 1:N_t
        b = autocorrelation(i,n_p)/autocorrelation(0,n_p)
        if b > 0
            a = a + b
        else
            break
        end
    end
    0.5 + a
end

x = []
for i = 0:N_t
    push!(x, i)
end

y_0 = abs.(corrs[1])
y_1 = abs.(corrs[2])
y_2 = abs.(corrs[3])
y_3 = abs.(corrs[4])
y_4 = abs.(corrs[5])

image = plot(x, y_0, size = (750,600))
image = plot!(x, y_1)
image = plot!(x, y_2)
image = plot!(x, y_3)
image = plot!(x, y_4)
image = plot!(yscale = :log)
display(image)


x_energies_0 = Real[]
x_energies_1 = Real[]
x_energies_2 = Real[]
x_energies_3 = Real[]
x_energies_4 = Real[]
exp_energies_0 = Real[]
exp_energies_1 = Real[]
exp_energies_2 = Real[]
exp_energies_3 = Real[]
exp_energies_4 = Real[]
for i = 1:N_t
    I = i%N_t +1         #  ↓  hier sind die period. Randbed.
    if corrs[1][i]/corrs[1][I] > 0
        push!(x_energies_0, i-0.5)
        push!(exp_energies_0, E_exp(i,0))
    end                  #  ↓
    if corrs[2][i]/corrs[2][I] > 0
        push!(x_energies_1, i-0.5)
        push!(exp_energies_1, E_exp(i,1))
    end                  #  ↓
    if corrs[3][i]/corrs[3][I] > 0
        push!(x_energies_2, i-0.5)
        push!(exp_energies_2, E_exp(i,2))
    end                  #  ↓
    if corrs[4][i]/corrs[4][I] > 0
        push!(x_energies_3, i-0.5)
        push!(exp_energies_3, E_exp(i,3))
    end                  #  ↓
    if corrs[5][i]/corrs[5][I] > 0
        push!(x_energies_4, i-0.5)
        push!(exp_energies_4, E_exp(i,4))
    end
end

y_exp_0 = (exp_energies_0)
y_exp_1 = (exp_energies_1)
y_exp_2 = (exp_energies_2)
y_exp_3 = (exp_energies_3)
y_exp_4 = (exp_energies_4)
#y_pt = (pt_energies_0)
image2 = plot(x_energies_0, y_exp_0, seriestype = :scatter, size = (750,600))
image2 = plot!(x_energies_1, y_exp_1, seriestype = :scatter)
image2 = plot!(x_energies_2, y_exp_2, seriestype = :scatter)
image2 = plot!(x_energies_3, y_exp_3, seriestype = :scatter)
image2 = plot!(x_energies_4, y_exp_4, seriestype = :scatter)
image2 = plot!(legend = :topright)
display(image2)


#=
# Histogram über alle Feldwerte
field_values = Real[]
for i = 1:nb_configs
    for j = 1:N_t
        for k = 1:N_s
            push!(field_values, configs[i][j][k])
        end
    end
end
H = StatsBase.fit(Histogram, field_values; nbins=60)
image_H = plot(H, label = "Häufigkeiten Feldwerte")

model(x, p) = p[1]*exp.(-p[2]^(2).*(x).^2)
p0 = [3000., 1.]
xvalues = Array(H.edges[1])
yvalues = H.weights + zeros(length(H.weights))
pop!(xvalues)
fit = LsqFit.curve_fit(model, xvalues, yvalues, p0)
p = fit.param
image_H = plot!(-0.00004:0.000001:0.00004, model(-0.00004:0.000001:0.00004, p), label = "Gauß-Fit")
display(image_H)
=#


#=
x_energies_pt_0 = Real[]
x_energies_pt_1 = Real[]
x_energies_pt_2 = Real[]
x_energies_pt_3 = Real[]
x_energies_pt_4 = Real[]
pt_energies_0 = Real[]
pt_energies_1 = Real[]
pt_energies_2 = Real[]
pt_energies_3 = Real[]
pt_energies_4 = Real[]

for i = 1:N_t
    a = (i+1+N_t)%N_t +1    # τ+Δτ
    b = (i-1+N_t)%N_t +1    # τ-Δτ
    c = i                   # τ

    if (corrs[1][a] + corrs[1][b]) / (2*corrs[1][c]) > 1
        push!(x_energies_pt_0, c)
        push!(pt_energies_0, E_pt(i,0))
    end
    #=
    if (corrs[2][a] + corrs[2][b]) / (2*corrs[2][c]) > 1
        push!(x_energies_pt_1, c)
        push!(pt_energies_1, E_pt(i,1))
    end
    if (corrs[3][a] + corrs[3][b]) / (2*corrs[3][c]) > 1
        push!(x_energies_pt_2, c)
        push!(pt_energies_2, E_pt(i,2))
    end
    if (corrs[4][a] + corrs[4][b]) / (2*corrs[4][c]) > 1
        push!(x_energies_pt_3, c)
        push!(pt_energies_3, E_pt(i,3))
    end
    =#
    if (corrs[5][a] + corrs[5][b]) / (2*corrs[5][c]) > 1
        push!(x_energies_pt_4, c)
        push!(pt_energies_4, E_pt(i,4))
    end
end

y_pt_0 = (pt_energies_0)
y_pt_1 = (pt_energies_1)
y_pt_2 = (pt_energies_2)
y_pt_3 = (pt_energies_3)
y_pt_4 = (pt_energies_4)

image3 = plot(x_energies_pt_0, y_pt_0, seriestype = :scatter, size = (750,600))
image3 = plot!(x_energies_pt_1, y_pt_1, seriestype = :scatter)
image3 = plot!(x_energies_pt_2, y_pt_2, seriestype = :scatter)
image3 = plot!(x_energies_pt_3, y_pt_3, seriestype = :scatter)
image3 = plot!(x_energies_pt_4, y_pt_4, seriestype = :scatter)
image3 = plot!(legend = false)
display(image3)


#=
function compare(x)
    4*sinh(x/2)^2
end


actions = []
for i = 1:nb_configs
    push!(actions, S(configs[i]))
end
