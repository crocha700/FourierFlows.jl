using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

 n, L  =  256, 2π
 ν, nν = 0e-3, 1
 μ, nμ = 1e-1, 0
dt, tf = 1e-2, 2/μ

nt = round(Int, tf/dt)
ns = 20

# Forcing
kf, dkf = 15, 2
σ = 0.1*1000

g  = TwoDGrid(n, L)

force2k = exp.(-(sqrt.(g.Kr.^2+g.Lr.^2)-kf).^2/(2*dkf^2))
σ0 = sum(force2k./(2*(g.Kr.^2+g.Lr.^2)+1e-15))
force2k .= σ/σ0 * force2k

srand(1234)


function calcF!(F, sol, t, s, v, p, g)
  if t == s.t # not a substep
    eta = (randn(size(sol)) + im*randn(size(sol)))/(sqrt(2)*s.dt)
    F .= eta.*sqrt.(p.force2k)
  end

  nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt,
  stepper="RK4", calcF=calcF!, force2k=force2k)

E = Diagnostic(energy,      prob, nsteps=nt)
D = Diagnostic(dissipation, prob, nsteps=nt)
R = Diagnostic(drag,        prob, nsteps=nt)
I = Diagnostic(injection,   prob, nsteps=nt)
diags = [E, D, I, R]

filename = @sprintf("stochastictest_kf%d.jld2", kf)
getsol(prob) = deepcopy(prob.state.sol)
out = Output(prob, filename, (:sol, getsol))

function makeplot(prob, diags)
  E, D, I, R = diags

  TwoDTurb.updatevars!(prob)

  E, D, I, R = diags

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

  residual
end


fig, axs = subplots(ncols=3, nrows=1, figsize=(13, 4))

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

  savename = @sprintf("./plots/stochastictest_kf%d_%06d.png", kf, prob.step)
  savefig(savename, dpi=240)
end

savediagnostic(E, "energy", out.filename)
savediagnostic(D, "dissipation", out.filename)
savediagnostic(I, "injection", out.filename)
savediagnostic(R, "drag", out.filename)
