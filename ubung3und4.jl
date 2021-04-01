using Plots
using StatsBase

# Simulationsparameter
m = 0.4     # Masse des freien Teilchens
d_t = 1     # Räumliche Dimensionen
d_s = 1     # Temporale Dimensionen
N_t = 8     # Anzahl räumlicher Gitterpunkte
N_s = 8     # Anzahl temporaler Gitterpunkte
m = 50000     # Anzahl Wiederholungen in Metropolis, also ∃ m+1 Konfigurationen
r = 0.001     # Schrittweite in Metropolis
cut = 25000    # Anzahl der ersten Konfig., die verworfen werden (Thermalisation)
N_skip = 250  # Es wird nur jede N_skip-te Konfig. gemessen


# Startkonfiguration generieren
# Der erste Index steht für temporale, der zweite für räumliche Koordinaten
x_0 = []
for i = 1:N_t
    q = 5*(rand(N_s) - 0.5*ones(N_s))
    push!(x_0, q)
end
x_0 = 3*x_0

# Funktion um Wirkung einer Konfiguration auszurechnen (vgl. Gl(56) im Skript)
function S(x)
    a = 0
    for i = 1:N_t
        I = (i)%N_t +1      # period. Randbed., s.u.

        for j in 1:N_s
            J = (j)%N_s +1  # period. Randbed., nämlich:           ↓            ↓
            a = a + (0.5*m^2 + d_s+d_t) * x[i][j]^2 - x[i][j] * (x[I][j] + x[i][J])
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
# Variante des Metropolis: anstatt alle 8×8 Einträge zu ändern, wird
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
            b = (test[j][k])
            #println(test)
            proposal[j][k] = proposal[j][k] + r*q
            #println(test)
            a = (proposal[j][k])
            #println(test)

            J = j%N_t +1            # ̂τ+Δτ, jeweils mit period. Randbed.
            K = k%N_s +1            # ⃗x+̂μ
            Y = (N_t+j-2)%N_t +1    # ̂τ-Δτ
            G = (N_s+j-2)%N_s +1    # ⃗x-̂μ

            Delta = (0.5*m^2+d_s+d_t)*(a^2-b^2) -
            (a-b)*(test[j][K] + test[J][k] + test[j][G] + test[Y][k])

            P = rand(1)[1]
            #println(test, "\n")
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
    a/(2*nb_configs*N_t)
end


corrs_0 = []
corrs_1 = []
corrs_2 = []
corrs_3 = []
corrs_4 = []
for i = 0:N_t
    push!(corrs_0, correlation(i, 0))
    push!(corrs_1, correlation(i, 1))
    push!(corrs_2, correlation(i, 2))
    push!(corrs_3, correlation(i, 3))
    push!(corrs_4, correlation(i, 4))
end


function E_exp(n, n_p)
    if n_p == 0
        exp_corrs = abs.(corrs_0)
    elseif n_p == 1
        exp_corrs = abs.(corrs_1)
    elseif n_p == 2
        exp_corrs = abs.(corrs_2)
    elseif n_p == 3
        exp_corrs = abs.(corrs_3)
    elseif n_p == 4
        exp_corrs = abs.(corrs_4)
    else
        println("Falscher Wert für n_p, es muss n_p ∈ {0,1,...,N_s/2}")
    end

    a = exp_corrs[n]/exp_corrs[n+1]
    if a > 0
        log(a)
    end
end


function E_pt(n, n_p)
    if n_p == 0
        pt_corrs = abs.(corrs_0)
    elseif n_p == 1
        pt_corrs = abs.(corrs_1)
    elseif n_p == 2
        pt_corrs = abs.(corrs_2)
    elseif n_p == 3
        pt_corrs = abs.(corrs_3)
    elseif n_p == 4
        pt_corrs = abs.(corrs_4)
    end

    a = (n+1+N_t)%N_t +1    # τ+Δτ
    b = (n-1+N_t)%N_t +1    # τ-Δτ
    c = (n+N_t)%N_t   +1    # τ

    A = (pt_corrs[a] + pt_corrs[b]) / (2*pt_corrs[c])
    if A >= 1
        acosh(A)
    end
end


x = []
for i = 0:N_t
    push!(x, i)
end

y_0 = real(corrs_0)
y_1 = real(corrs_1)
y_2 = real(corrs_2)
y_3 = real(corrs_3)
y_4 = real(corrs_4)

image = plot(x, y_0)
image = plot!(x, y_1)
image = plot!(x, y_2)
image = plot!(x, y_3)
image = plot!(x, y_4)
#image = plot!(yscale = :ln)
display(image)


x_energies = Real[]
exp_energies_0 = Real[]
exp_energies_1 = Real[]
exp_energies_2 = Real[]
for i = 1:N_t
    push!(x_energies, i-1)
    push!(exp_energies_0, E_exp(i,0))
    push!(exp_energies_1, E_exp(i,1))
    push!(exp_energies_2, E_exp(i,2))
end

y_exp_0 = (exp_energies_0)
y_exp_1 = (exp_energies_1)
y_exp_2 = (exp_energies_2)
#y_pt = (pt_energies_0)
image2 = plot(x_energies, y_exp_0)
image2 = plot!(x_energies, y_exp_1)
image2 = plot!(x_energies, y_exp_2)
display(image2)


function compare(x)
    4*sinh(x/2)^2
end
