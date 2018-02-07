using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

 n, L  =  128, 2π
 ν, nν = 1e-6, 2
 μ, nμ = 1e-1, 0
dt, tf = 0.01, 100

nt = round(Int, tf/dt)

nt = 20000
ns = 50

# Forcing
kf, dkf = 12.0, 1.0
σ = 0.005

gr  = TwoDGrid(n, L)

force2k = exp.(-(sqrt.(gr.KKrsq)-kf).^2/(2*dkf^2))
force2k[gr.KKrsq .< 2.0^2 ] = 0
force2k[gr.KKrsq .> 20.0^2 ] = 0
force2k[gr.Kr.<2π/L] = 0
# σ0 = sum(force2k.*gr.invKKrsq/2.0)
σ0 = FourierFlows.parsevalsum(force2k.*gr.invKKrsq/2.0, gr)/(gr.Lx*gr.Ly)
force2k .= σ/σ0 * force2k

force2k_fft = exp.(-(sqrt.(gr.KKsq)-kf).^2/(2*dkf^2))
force2k_fft[gr.KKsq .< 2.0^2 ] = 0
force2k_fft[gr.KKsq .> 20.0^2 ] = 0
force2k_fft[abs.(gr.K).<2π/L] = 0
σ0_fft = sum(force2k_fft.*gr.invKKsq/2.0)
force2k_fft .= σ/σ0_fft * force2k_fft * (gr.nx*gr.ny)^2

# if size(force2k)[1]==g.nkr
#     force2k .= 2*force2k
# end


srand(1234)

function calcF!(F, sol, t, s, v, p, g)
  if t == s.t # not a substep
    # eta = (randn(size(sol)) + im*randn(size(sol)))/(sqrt(2)*sqrt(s.dt))

    eta = exp.(2π*im*rand(size(sol)))/sqrt(s.dt)
    eta[1, :] .= 0
    @. F = eta .* sqrt(force2k)
    # Fphys = rfft(irfft(F, g.nx))
    # F .= Fphys

    # eta = exp.(2π*im*rand(size(v.q)))/sqrt(s.dt)
    # F1 = eta.*sqrt(force2k_fft)
    # F1[1, 1] = 0
    # Fphys = real.(ifft(F1))
    # F2 = rfft(Fphys)
    # F .= F2
  end
  nothing
end

prob = TwoDTurb.ForcedProblem(nx=n, Lx=L, ν=ν, nν=nν, μ=μ, nμ=nμ, dt=dt,
  stepper="RK4", calcF=calcF!)

s, v, p, g, eq, ts = prob.state, prob.vars, prob.params, prob.grid, prob.eqn, prob.ts;

p.calcF!(v.F, s.sol, s.t, s, v, p, g)


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

  TwoDTurb.updatevars!(prob)

  E, D, I, R = diags

  sca(axs[1]); cla()
  pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)
  xlabel(L"x")
  ylabel(L"y")

  sca(axs[2]); cla()

  i₀ = 1
  dEdt = (E[(i₀+1):E.count] - E[i₀:E.count-1])/prob.ts.dt
  ii = (i₀):E.count-1

  # dEdt = I - D - R?
  total = I[ii] + σ - D[ii] - R[ii]
  residual = dEdt - total

  # plot(E.time[ii], I[ii], "o", markersize=0.5, label="injection (\$I\$)")
  # plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  # plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  # plot(E.time[ii], residual, "c-", label="residual")
  # plot(E.time[ii], dEdt[ii], label="dissipation (\$dE/dt\$)")
  plot(E.time[ii], I[ii] + σ, label="injection (\$I\$)")
  plot(E.time[ii], σ*exp.(0*E.time[ii]), "--")
  plot(E.time[ii], -D[ii], label="dissipation (\$D\$)")
  plot(E.time[ii], -R[ii], label="drag (\$R\$)")
  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10)

  sca(axs[3]); cla()
  # plot(E.time[ii], residual, "c-", label="residual")
  plot(E.time[ii], total, label="computed")
  plot(E.time[ii], dEdt, "--k", label="numerical")

  ylabel("Energy sources and sinks")
  xlabel(L"t")
  legend(fontsize=10)

  sca(axs[4]); cla()
  plot(E.time[ii], E[ii])
  plot(E.time[ii], σ/(2*μ)*(1-exp.(-2*μ*E.time[ii])), "--")
  # plot(E.time[ii], 0.1*E.time[ii], label="predicted (\$E\$)")
  xlabel(L"t")
  ylabel(L"E")

  println(mean(I[ii]))

  residual
end

fig, axs = subplots(ncols=2, nrows=2, figsize=(10, 8))

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
