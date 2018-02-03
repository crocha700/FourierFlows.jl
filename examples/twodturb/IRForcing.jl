using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

 n, L  =  128, 2π
 ν, nν = 0e-3, 1
 μ, nμ = 0e-1, 0
dt, tf = 1e-2, 100

nt = round(Int, tf/dt)

nt = 1
ns = 1

# Forcing
kf, dkf = 12.0, 2.0
σ = 0.1

gr  = TwoDGrid(n, L)

force2k = exp.(-(sqrt.(gr.Kr.^2+gr.Lr.^2)-kf).^2/(2*dkf^2))
force2k[sqrt.(gr.Kr.^2 + gr.Lr.^2) .< 2.0 ]=0
force2k[sqrt.(gr.Kr.^2 + gr.Lr.^2) .> 20.0 ]=0
iksq = 1./(gr.Kr.^2 + gr.Lr.^2)
iksq[1, 1] = 0
σ0 = sum(force2k.*iksq/1.0)
force2k .= σ/σ0 * force2k
# if size(force2k)[1]==g.nkr
#     force2k .= 2*force2k
# end


srand(1234)

function calcF!(F, sol, t, s, v, p, g)
  if t == s.t # not a substep
    # eta = (randn(size(sol)) + im*randn(size(sol)))/(sqrt(2)*sqrt(s.dt))
    # F .= eta.*sqrt.(force2k)*(g.nx*g.ny)
    F .= 0.0
    F[2, 3] .= 1
  end
  nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt,
  stepper="ForwardEuler", calcF=calcF!)

s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;

TwoDTurb.set_q!(prob, 0*g.X)


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

  # plot(E.time[ii], I[ii], "o", markersize=0.5, label="injection (\$I\$)")
  # plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  # plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  # plot(E.time[ii], residual, "c-", label="residual")
  # plot(E.time[ii], dEdt[ii], label="dissipation (\$dE/dt\$)")
  plot(E.time[ii], I[ii], label="injection (\$I\$)")
  # plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  # plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  # plot(E.time[ii], residual, "c-", label="residual")

  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10)

  sca(axs[3]); cla()
  plot(E.time[ii], E[ii])
  # plot(E.time[ii], 0.1*E.time[ii], label="predicted (\$E\$)")
  xlabel(L"t")
  ylabel(L"E")

  residual
end


# fig, axs = subplots(ncols=3, nrows=1, figsize=(13, 4))

# Step forward
for i = 1:ns
  tic()
  stepforward!(prob, diags, round(Int, nt/ns))
  tc = toq()

  TwoDTurb.updatevars!(prob)
  # saveoutput(out)

  cfl = prob.ts.dt*maximum(
    [maximum(prob.vars.V)/prob.grid.dx, maximum(prob.vars.U)/prob.grid.dy])

  # res = makeplot(prob, diags)
  # pause(0.1)

  # @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s, mean(res) = %.3e\n",
  #   prob.step, prob.t, cfl, tc, mean(res))
  @printf("step: %04d, t: %.1f, cfl: %.3f, time: %.2f s\n",
    prob.step, prob.t, cfl, tc)

  # savename = @sprintf("./plots/stochastictest_kf%d_%06d.png", kf, prob.step)
  # savefig(savename, dpi=240)
end

println(s.sol[2, 3])
println(eq.LC[2, 3])


# savediagnostic(E, "energy", out.filename)
# savediagnostic(D, "dissipation", out.filename)
# savediagnostic(I, "injection", out.filename)
# savediagnostic(R, "drag", out.filename)
