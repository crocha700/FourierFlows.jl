using PyPlot, FourierFlows
using VerticallyFourierBoussinesq

  n,    L  =  256, 2π
nu0, nnu0 = 1e-3,  1
nu1, nnu1 = 1e-3,  1
mu0, nmu0 = 1e-1, -1
mu1, nmu1 = 1e-1, -1
 dt,   tf = 2e-3, 1000

nt = round(Int, tf/dt)
ns = 100

# Forcing
P = 1.0
K = 16
amplitude = P*K/sqrt(dt) * n^2/2

function calcF!(F, sol, t, s, v, p, g)
  if t == s.t # not a substep
    F .= 0.0
    θk = 0.0 #2π*rand() 
    phase = σ*t

    i₁ = round(Int, abs(K*cos(θk))) + 1
    j₁ = round(Int, abs(ki*sin(θk))) + 1  # j₁ >= 1
    F[i₁, j₁] = amplitude*exp(im*phase)
  end

  nothing
end

prob = Problem(nx=n, Lx=L, nu0=nu0, nnu0=nnu0, nu1=nu1, nnu1=nnu1, 
  mu0=mu0, nmu0=nmu0, mu1=mu1, nmu1=nmu1, dt=dt, calcF=calcF!, stepper="RK4")
  
filename = @sprintf("stochastictest_ki%d.jld2", ki)
getsol(prob) = deepcopy(prob.state.sol)
out = Output(prob, filename, (:sol, getsol))

function makeplot(prob, diags)
  E, D, I, R = diags

  TwoDTurb.updatevars!(prob)  

  close("all")
  E, D, I, R = diags
  fig, axs = subplots(ncols=3, nrows=1, figsize=(13, 4))

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  xlabel(L"x")
  ylabel(L"y")

  sca(axs[2]); cla()

  i₀ = 1
  dEdt = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀+1):E.count

  # dEdt = I - D - R?
  total = I[ii] - D[ii] - R[ii]
  residual = dEdt - total

  plot(E.time[ii], I[ii], "o", markersize=0.5, label="injection (\$I\$)")
  plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  plot(E.time[ii], residual, "c-", label="residual")

  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10)

  sca(axs[3]); cla()
  plot(E.time[ii], E[ii])
  xlabel(L"t")
  ylabel(L"E")

  tight_layout()

  residual
end

# Step forward
for i = 1:ns
  tic()
  stepforward!(prob, diags, round(Int, nt/ns))
  tc = toq()

  TwoDTurb.updatevars!(prob)  
  saveoutput(out)

  cfl = prob.ts.dt*maximum(
    [maximum(prob.vars.V)/prob.grid.dx, maximum(prob.vars.U)/prob.grid.dy])

  res = makeplot(prob, diags)
  pause(0.1)

  @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s, mean(res) = %.3e\n", 
    prob.step, prob.t, cfl, tc, mean(res))
    
  savename = @sprintf("./plots/stochastictest_ki%d_%06d.png", ki, prob.step)
  savefig(savename, dpi=240)
end

savediagnostic(E, "energy", out.filename)
savediagnostic(D, "dissipation", out.filename)
savediagnostic(I, "injection", out.filename)
savediagnostic(R, "drag", out.filename)
