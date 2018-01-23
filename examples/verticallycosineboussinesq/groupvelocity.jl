using PyPlot, FourierFlows, JLD2, FourierFlows.VerticallyCosineBoussinesq

  n,    L =  128, 2π   # Domain
nu0, nnu0 = 1e-6,   1  # Viscosity
nu1, nnu1 = 1e-6,   1  # Viscosity
mu0, nmu0 =  0.0,   0
mu1, nmu1 =  0.0,   0
  f, N, m =  1.0, 1.0, 4.0
   uw, kw = 0.01,  16

# ε = U / σ L 
# σ = f sqrt[ 1 + (Nk/m)^2 ] = sqrt(2)-4
# tσ = 1/σ = 0.25 - 0.5
# tZ = Z^{-1}

 σ = f*sqrt(1 + (N*kw/m)^2)
tσ = 2π/σ
dt = tσ/100
ni = round(Int, tσ/dt)
ns = 100
nt = ns*ni

makesquare!(ax) = ax[:set_aspect](1, adjustable="box")
makesquare!(axs::AbstractArray) = for ax in axs; makesquare!(ax); end

ticksoff!(a) = a[:tick_params](bottom=false, left=false, labelbottom=false, 
  labelleft=false)
ticksoff!(axs::AbstractArray) = for ax in axs; ticksoff!(ax); end

e1density(u, v, p, m, N) = @. ( u^2 + v^2 + m^2*p^2/N^2 )/2
e1density(prob) = e1density(prob.vars.u, prob.vars.v, prob.vars.p, 
  prob.params.m, prob.params.N)

weightedke(u, v, σ, f) = @. u^2 + (σ/f*v)^2
weightedke(prob) = weightedke(prob.vars.u, prob.vars.v, σ, f)

function makeplot!(ax, prob)
  sca(ax); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, e1density(prob))
  makesquare!(ax)
  ticksoff!(ax)
  nothing
end

function makeplot!(ax, prob, xw, yw)
  makeplot!(ax, prob)
  plot(xw, yw, "ro")
  nothing
end

function makesanityplot!(axs, prob)
  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, weightedke(prob), vmin=0.0)
  sca(axs[2]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.p)
  makesquare!(axs)
  ticksoff!(axs)
  nothing
end


wavecentroid(prob) = (FourierFlows.xmoment(e1density(prob), prob.grid), 
  FourierFlows.ymoment(e1density(prob), prob.grid))

prob = Problem(f=f, N=N, m=m, nx=n, Lx=L, nu0=nu0, nnu0=nnu0, nu1=nu1, 
  nnu1=nnu1, mu0=mu0, nmu0=nmu0, mu1=mu1, nmu1=nmu1, dt=dt, 
  stepper="FilteredRK4")
  
r = L/10
envelope(x, y) = exp(-x^2/(2*r^2))
set_planewave!(prob, uw, kw; envelope=envelope)

t₋₁ = prob.t
xw₋₁, yw₋₁ = wavecentroid(prob)

#fig, axs = subplots()
fig, axs = subplots(ncols=2)
for i = 1:ns
  @time stepforward!(prob, ni)
  updatevars!(prob)

  xw, yw = wavecentroid(prob)

  makesanityplot!(axs, prob)
  pause(0.1)

  cg = (xw-xw₋₁) / (prob.t-t₋₁)
  cg0 = N^2*kw/(σ*m^2)

  @printf "xw/L: %.2f, yw/L: %.2f, cg/cg0: %.4f" xw/L yw/L cg/cg0

  t₋₁ = prob.t
  xw₋₁, yw₋₁ = xw, yw
end
