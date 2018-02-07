using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

 n, L  = 128, 2π
 ν, nν = 1e-6, 2
 μ, nμ = 1e-1, 0
dt, tf = 0.01, 2/μ
nt = round(Int, tf/dt)
ns = 4

# Forcing
kf, dkf = 12.0, 2.0
σ = 0.005

gr  = TwoDGrid(n, L)

force2k = exp.(-(sqrt.(gr.KKrsq)-kf).^2/(2*dkf^2))
force2k[gr.KKrsq .< 2.0^2 ] = 0
force2k[gr.KKrsq .> 20.0^2 ] = 0
force2k[gr.Kr.<2π/L] = 0
# σ0 = sum(force2k.*gr.invKKrsq/2.0)
σ0 = FourierFlows.parsevalsum(force2k.*gr.invKKrsq/2.0, gr)/(gr.Lx*gr.Ly)
force2k .= σ/σ0 * force2k

srand(1234)

function calcF!(F, sol, t, s, v, p, g)
    eta = exp.(2π*im*rand(size(sol)))/sqrt(s.dt)
    eta[1, 1] = 0
    @. F = eta .* sqrt(force2k)
    nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt,
  stepper="RK4", calcF=calcF!)

s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;

TwoDTurb.set_q!(prob, 0*g.X)
E = Diagnostic(energy,      prob, nsteps=nt)
D = Diagnostic(dissipation, prob, nsteps=nt)
R = Diagnostic(drag,        prob, nsteps=nt)
I = Diagnostic(injection,   prob, nsteps=nt)
diags = [E, D, I, R]


function makeplot(prob, diags)

  TwoDTurb.updatevars!(prob)

  E, D, I, R = diags

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  xlabel(L"$x$")
  ylabel(L"$y$")
  axis("square")

  sca(axs[2]); cla()

  i₀ = 1
  dEdt = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀):E.count-1
  ii2 = (i₀+1):E.count

  # dEdt = I - D - R?

  # If the Ito interpretation was used for the work
  # then we need to add the drift term
  total = I[ii2]+σ - D[ii] - R[ii]      # Ito
  # total = I[ii2] - D[ii] - R[ii]        # Stratonovich


  residual = dEdt - total

  # If the Ito interpretation was used for the work
  # then we need to add the drift term
  plot(E.time[ii], I[ii2] + σ, label="injection (\$I\$)")   # Ito
  # plot(E.time[ii], I[ii2] , label="injection (\$I\$)")      # Stratonovich
  plot(E.time[ii], σ+0*E.time[ii], "--")
  plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  ylabel("Energy sources and sinks")
  xlabel(L"$t$")
  legend(fontsize=10)

  sca(axs[3]); cla()
  plot(E.time[ii], total[ii], label=L"computed $dE/dt$")
  plot(E.time[ii], dEdt, "--k", label=L"numerical $dE/dt$")
  ylabel(L"$dE/dt$")
  xlabel(L"$t$")
  legend(fontsize=10)

  sca(axs[4]); cla()
  plot(E.time[ii], residual, "c-", label=L"residual = computed $-$ numerical")
  xlabel(L"$t$")
  ylabel(L"$dE/dt$")
  legend(fontsize=10)

  residual
end

fig, axs = subplots(ncols=2, nrows=2, figsize=(12, 8))

# Step forward
for i = 1:ns
  tic()

  stepforward!(prob, diags, round(Int, nt/ns))
  tc = toq()

  TwoDTurb.updatevars!(prob)
  # saveoutput(out)

  cfl = prob.ts.dt*maximum(
    [maximum(v.V)/g.dx, maximum(v.U)/g.dy])

  res = makeplot(prob, diags)
  pause(0.01)

  # @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s, mean(res) = %.3e\n",
  #   prob.step, prob.t, cfl, tc, mean(res))
  @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s\n",
    prob.step, prob.t, cfl, tc)

  # savename = @sprintf("./plots/stochastictest_kf%d_%06d.png", kf, prob.step)
  # savefig(savename, dpi=240)
end


# savediagnostic(E, "energy", out.filename)
# savediagnostic(D, "dissipation", out.filename)
# savediagnostic(I, "injection", out.filename)
# savediagnostic(R, "drag", out.filename)
