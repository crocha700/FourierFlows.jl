using PyPlot, FourierFlows, JLD2
import FourierFlows.VerticallyCosineBoussinesq

@load "testturb.jld2" q

nx, ny = size(q) 
L =  2π   # Domain
nu0, nnu0 = 1e-6,   1  # Viscosity
nu1, nnu1 = 1e-6,   1  # Viscosity
mu0, nmu0 =  0.0,   0
mu1, nmu1 =  0.0,   0
  f, N, m =  1.0, 1.0, 4.0
   uw, kw = 0.02, 16

# ε = U / σ L 
# σ = f sqrt[ 1 + (Nk/m)^2 ] = sqrt(2)-4
# tσ = 1/σ = 0.25 - 0.5
# tZ = Z^{-1}

 σ = f*sqrt(1 + (N*kw/m)^2)
tσ = 2π/σ
dt = tσ/100

@printf "σ/f: %.3f, ε: %.3f" σ/f maximum(abs.(q))/σ

function makesquare!(axs)
  for ax in axs
    ax[:set_aspect](1, adjustable="box")
    ax[:set_aspect](1, adjustable="box")
  end
  nothing
end

function makeplot!(axs, prob, diags=nothing)
  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.Z)

  sca(axs[2]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, real.(prob.vars.u))
  makesquare!(axs[1:2])

  nothing
end

prob = VerticallyCosineBoussinesq.Problem(f=f, N=N, m=m,
  nx=nx, Lx=L, nu0=nu0, nnu0=nnu0, nu1=nu1, nnu1=nnu1, mu0=mu0, nmu0=nmu0, 
  mu1=mu1, nmu1=nmu1, dt=dt, stepper="FilteredRK4")

VerticallyCosineBoussinesq.set_Z!(prob, q)
VerticallyCosineBoussinesq.set_planewave!(prob, uw, kw)

fig, axs = subplots(ncols=2, figsize=(7,4))

for i = 1:100
  @time stepforward!(prob, 100)
  makeplot!(axs, prob)
  pause(0.1)
end
