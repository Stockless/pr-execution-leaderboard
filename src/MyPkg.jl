module MyPkg

t0 = 0.
tf = 3.
γ = 2 * π * 42.58 * 1e6  # rad s^-1 T^-1
T1 = 1.  # s
T2 = 0.5  # s
B = [0.; 0.; 1e-7]  # T
m0 = [1.; 0.; 0.]

struct Theoretical
end

"Returns the theoretical solution of the Bloch equations for a given initial magnetization m0, time step dt, maximum time tmax."
function solve(m0, dt, tmax, method::Theoretical)
    t = dt:dt:tmax
    M0 = m0[1]
    Bz = B[3]
    Mx = M0 * cos.(γ * Bz * t) .* exp.(-t / T2)
    My = -M0 * sin.(γ * Bz * t) .* exp.(-t / T2)
    Mz = M0 * (1 .- exp.(-t / T1))
    return [Mx'; My'; Mz']
end

"Simulates the Bloch equations for a given initial magnetization m0, time step dt, maximum time tmax and a simulation method."
function solve(m0, dt, tmax, method)
    Nsteps = Int(tmax / dt)
    m = m0
    mt = zeros(eltype(m), (length(m), Nsteps))
    for i in 1:Nsteps
        m = step(dt, m, method)
        mt[:, i] = m
    end
    return mt
end

struct ForwardEuler
end

"Returns the magnetization vector at the next time step using the Forward Euler method."
function step(dt, m, method::ForwardEuler)
    return m .+ dt * bloch(m)
end

"Returns the time derivative of the magnetization vector m."
function bloch(m)
    Mx, My, Mz = m
    Bx, By, Bz = B
    M0 = m0[1]
    crossmB = [My * Bz - Mz * By; Mz * Bx - Mx * Bz; Mx * By - My * Bx]
    return γ * crossmB - [Mx / T2; My / T2; (Mz - M0) / T1]
end

export solve, step, ForwardEuler, Theoretical

end # module MyPkg
