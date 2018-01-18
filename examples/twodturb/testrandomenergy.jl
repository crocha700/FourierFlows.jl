using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, dissipation, drag

function energytest(dt; n=128, L=2π, ν=1e-3, nν=1, tf=10, q0=nothing)

  nt = round(Int, tf/dt)
  prob = TwoDTurb.Problem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt, stepper="RK4")

  if q0 == nothing
    q00 = cos.(prob.grid.X).*cos.(prob.grid.Y)
  elseif typeof(q0) <: Function
    q00 = q0.(prob.grid.X, prob.grid.Y) 
  else
    q00 = q0
  end

  TwoDTurb.set_q!(prob, q00)

  E = Diagnostic(energy, prob, nsteps=nt) 
  D = Diagnostic(dissipation, prob, nsteps=nt)
  R = Diagnostic(drag, prob, nsteps=nt)
  diags = [E, D, R]

  stepforward!(prob, diags, nt)

  ii = 2:(E.count-1)
  ii₊₁ = 3:E.count
  dEdt = (E[ii₊₁] - E[ii]) / dt

  #?? dedtnorm = map(x->maximum(x), zip(abs.(D[ii]), abs.(R[ii]), abs.(dEdt)))
  dedtnorm = zeros(length(ii))
  for (i, j) in enumerate(ii)
    dedtnorm[i] = maximum(abs.([dEdt[i], D[j], R[j]]))
  end
  res = abs.(dEdt .+ D[ii] .+ R[ii]) ./ dedtnorm

  E.time[ii], E[ii], D[ii], R[ii], dEdt, res
end

n, L = 128, 2π
q0 = rand(n, n)
#q0(x, y) = 0.01*cos(x)*cos(y)

#dts = 10.0.^[-1, -2, -3, -4]
dts = logspace(-3, -0.5, 10)
finalres = []

fig, axs = subplots(ncols=3, figsize=(10, 3))

for (i, dt) = enumerate(dts)
  t, E, D, R, dEdt, res = energytest(dt; n=n, L=L, q0=q0, tf=1.0)
  sca(axs[1])
  plot(t, E)
  sca(axs[2])
  plot(t, res)
  push!(finalres, res[end])
end

sca(axs[3])
plot(dts, finalres ./ dts, "o")
xlabel(L"\Delta t")
ylabel("residual / \$\\Delta t\$")

sca(axs[1])
xlabel(L"t")
ylabel(L"E")

sca(axs[2])
xlabel(L"t")
ylabel(L"|E_t - R - D|/\mathrm{max}(E_t, D, R)")

axs[1][:set_yscale]("log")
axs[2][:set_yscale]("log")
axs[3][:set_xscale]("log")
#axs[3][:set_yscale]("log")

tight_layout()
show()
