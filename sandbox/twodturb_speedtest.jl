include("../src/fourierflows.jl")

using PyPlot, FourierFlows
using FourierFlows.TwoDTurb

nx = 1024
dt = 1e-4   
nsub = 100
ntot = 1000

prob = TwoDTurb.InitialValueProblem(nx=nx, dt=dt, withfilter=true)
set_q!(prob, rand(nx, nx))

println("Stepping forward...")
while prob.step < ntot
  @time stepforward!(prob, nsteps=nsub)
end

updatevars!(prob)

fig, ax = subplots()
imshow(prob.vars.q)
show()
