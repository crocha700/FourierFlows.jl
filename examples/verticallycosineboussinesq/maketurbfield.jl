using PyPlot, FourierFlows, JLD2
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation

 n,  L = 128, 2π   # Domain
 ν, nν = 1e-6, 1  # Viscosity
dt, nt = 1.0, 200   # Time step
ns = 10

prob = TwoDTurb.Problem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt, stepper="FilteredRK4")
TwoDTurb.set_q!(prob, 1.0*rand(n, n))
x, y = prob.grid.X, prob.grid.Y

# Step forward
fig = figure()
tic()
for i = 1:10
  stepforward!(prob, nt)
  TwoDTurb.updatevars!(prob)  

  cfl = maximum(prob.vars.U)*prob.ts.dt/prob.grid.dx
  @printf("step: %04d, t: %6.1f, cfl: %.2f, ", prob.step, prob.t, cfl)
  toc(); tic()

  clf()
  pcolormesh(x, y, prob.vars.q)
  gca()[:set_aspect](1, adjustable="box")
  colorbar()
  pause(0.01)
end

q = prob.vars.q
@save "testturb_n128_2.jld2" q
