#=
Die Phasen von Eichfeld-Konfigurationen Uμ( ⃗n) werden in einem Array
config[μ, nₜ, nₛ] gespeichert, wobei ⃗n = (nₜ, nₛ). Im weiteren soll
μ = 1 für die zeitliche und μ = 2 für die räumliche Dimension stehen.
(es wird auch μ ≡ 1 für temp. und ν ≡ 2 für räuml. synonym verwendet)
Phasen der Eichtransformationen Ω( ⃗n) werden als trans[nₜ, nₛ] gespeichert.
=#
N_t = 8
N_s = 8
beta = 1

# Funktion um Plaquette Pμν( ⃗n) einer geg. Konfig. zu bestimmen
function P12(nt, ns, config)
    NT = nt%N_t +1  # periodische Randbed.
    NS = ns%N_s +1
    config[1][nt][ns] + config[2][NT][ns] - config[1][nt][NS] - config[2][nt][ns]
end


# Funktion um Plaquette Pνμ( ⃗n) einer geg. Konfig. zu bestimmen
function P21(nt, ns, config)
    NT = nt%N_t +1  # periodische Randbed.
    NS = ns%N_s +1
    config[2][nt][ns] + config[1][NT][ns] - config[2][nt][NS] - config[1][nt][ns]
end


# Funktion um eine gegebene Konfig. unter gegebener Trans. zu transformieren
function transform(config, trans)
    transformed = [[],[]]   # [[μ],[ν]]
    for i = 1:N_t
        push!(transformed[1], [])
        push!(transformed[2], [])
        for j = 1:N_s
            I = i%N_t +1
            J = j%N_s +1
            a = trans[i][j] + config[1][i][j] - trans[I][j]
            b = trans[i][j] + config[2][i][j] - trans[i][J]
            push!(transformed[1][i], a)
            push!(transformed[2][i], b)
        end
    end
    transformed
end

# Funktion um die Wilson Gauge Action einer gegebenen Konfig. zu bestimmen
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


# Aufgabenteil d)
test = []
for mu = 1:2
    push!(test, [])
    for i = 1:N_t
        push!(test[mu], ones(N_s))
    end
end

Omega = []
for i = 1:N_t
    a = (2*pi).*rand(N_s)
    push!(Omega, a)
end

a = transform(test, Omega)
println("", "Wirkung der 1er Konfig: ",  S(test), ",   Wirkung transf.:", S(a), "")

# Aufgabenteil e)
test2 = []
for mu = 1:2
    push!(test2, [])
    for i = 1:N_t
        a = (2*pi).*rand(N_s)
        push!(test2[mu], a)
    end
end

b = transform(test2, Omega)
println("", "Wirkung der 0-2π Konfig: ",  S(test2), ",   Wirkung transf.:", S(b), "")
