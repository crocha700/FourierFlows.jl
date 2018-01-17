using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation

cfl(prob) = maximum([maximum(abs.(prob.vars.U)), maximum(abs.(prob.vars.V))]*
              prob.ts.dt/prob.grid.dx)

function lambdipoleexperiment(n, dt; L=2π, Ue=1, Re=L/20, ν=0, nν=1, 
  ti=L/Ue*0.01, nm=3)

  nt = round(Int, ti/dt)

  prob = TwoDTurb.InitialValueProblem(nx=n, Lx=L, ν=ν, nν=nν, dt=dt,
    stepper="FilteredRK4")
  x, y, q = prob.grid.X, prob.grid.Y, prob.vars.q # nicknames

  q0 = FourierFlows.lambdipole(Ue, Re, prob.grid)
  TwoDTurb.set_q!(prob, q0)

  xq = zeros(nm)   # centroid of abs(q)
  Ue_m = zeros(nm) # measured dipole speed

  # Step forward
  for i = 1:nm
    tic()
    stepforward!(prob, nt)
    TwoDTurb.updatevars!(prob)  
    xq[i] = mean(abs.(q).*x) / mean(abs.(q))

    if i > 1
      Ue_m[i] = (xq[i]-xq[i-1]) / ((nt-1)*dt)
    else
      Ue_m[i] = 0.0
    end

    @printf("     step: %04d, t: %3.1f, time: %.3f, cfl: %.2f\n",
      prob.step, prob.t, toq(), cfl(prob))
  end

  Ue, Ue_m[2:end]
end

for n = [128, 256, 512]
  for dt = [1e-3, 5e-4, 2e-4]
    @printf "n: %d:\n" n
    Ue, Ue_m = lambdipoleexperiment(n, dt)
    @printf "n: %d, dt: %.1e, mean(Ue_m/Ue): %.4f\n" n dt mean(Ue_m/Ue)
  end
end
