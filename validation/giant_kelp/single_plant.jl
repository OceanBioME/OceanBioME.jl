using CairoMakie, Oceananigans.Units, OceanBioME

model = GiantKelp()

const year = years = 365days

@inline PAR(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline temp(t) = 2.4 * cos(t * 2π / year + 50day) + 10

A, N, C = 10, 0.01, 0.1

t = 0.0

te = 1000days
Δt = 10minutes

nt = floor(Int, te / Δt)

A_ts, N_ts, C_ts = zeros(nt), zeros(nt), zeros(nt)

for n in 1:nt
    A_ts[n] = A
    N_ts[n] = N
    C_ts[n] = C

    global t += Δt

    global A += model(Val(:A), t, A, N, C, 0, 0, 0, temp(t), 5, 2, PAR(t) * exp(-0.2 * 10)) * Δt
    global N += model(Val(:N), t, A, N, C, 0, 0, 0, temp(t), 5, 2, PAR(t) * exp(-0.2 * 10)) * Δt
    global C += model(Val(:C), t, A, N, C, 0, 0, 0, temp(t), 5, 2, PAR(t) * exp(-0.2 * 10)) * Δt

    C < 0.01 && @error "C went below min"
end