using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation

 n,  L = 256, 2π   # Domain
 ν, nν = 1e-6, 1  # Viscosity
dt, nt = 1.0, 200   # Time step
ns = 10

prob = TwoDTurb.Problem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt, stepper="FilteredRK4")
TwoDTurb.set_q!(prob, rand(n, n))

E = Diagnostic(energy,      prob, nsteps=nt*ns) 
Z = Diagnostic(enstrophy,   prob, nsteps=nt*ns) 
D = Diagnostic(dissipation, prob, nsteps=nt*ns)
diags = [E, Z, D]

# Step forward
fig = figure()
tic()
for i = 1:10
  stepforward!(prob, diags, nt)
  TwoDTurb.updatevars!(prob)  

  cfl = maximum(prob.vars.U)*prob.ts.dt/prob.grid.dx
  @printf("step: %04d, t: %6.1f, cfl: %.2f, ", prob.step, prob.t, cfl)
  toc(); tic()

  clf(); imshow(prob.vars.q); pause(0.01)
end


# Calculate radial energy spectrum
 E = 0.5*(prob.vars.U.^2 + prob.vars.V.^2) # energy density
Eh = rfft(E)
kr, Ehr = FourierFlows.radialspectrum(Eh, prob.grid, refinement=1)

fig, axs = subplots(ncols=2, figsize=(8, 4))

sca(axs[1]); cla()
pcolormesh(prob.grid.X, prob.grid.Y, prob.vars.q)

xlabel(L"x")
ylabel(L"y")
title("Vorticity")


sca(axs[2]); cla()
plot(kr, abs.(Ehr))

xlabel(L"k_r")
ylabel(L"\int | \hat{E} | k_r \mathrm{d} k_{\theta}")
title("Radial energy spectrum")

xlim(0, n/4)
axs[2][:set_yscale]("log")

tight_layout()
show()
