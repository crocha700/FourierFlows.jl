using FourierFlows#, PyPlot
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

 n, L  =  128, 2π
 ν, nν = 0e-3, 1
 μ, nμ = 1e-1, 0
dt = 0.1


gr  = TwoDGrid(n, L)

println(" ")
println("Forcing ONLY component [2, 3] with constant forcing F[2, 3] = 2 + i")
println(" ")

function calcF!(F, sol, t, s, v, p, g)
  if t == s.t # not a substep
    F .= 0.0
    F[2, 3] .= 2.0 + im
  end
  nothing
end

nx, ny = n, n
Lx, Ly = L, L
_calcF = calcF!

g1  = TwoDGrid(nx, Lx, ny, Ly)
pr1 = TwoDTurb.ForcedParams(ν, nν, μ, nμ, _calcF)
vs1 = TwoDTurb.ForcedVars(g1)
eq1 = TwoDTurb.Equation(pr1, g1)
# ts1 = FourierFlows.ForwardEulerTimeStepper(dt, deepcopy(eq1.LC))    # this works
ts1 = FourierFlows.ForwardEulerTimeStepper(dt, eq1.LC)    # this DOES NOT work
prob1 = FourierFlows.Problem(g1, vs1, pr1, eq1, ts1)
s1 = prob1.state

g2  = TwoDGrid(nx, Lx, ny, Ly)
pr2 = TwoDTurb.ForcedParams(ν, nν, μ, nμ, _calcF)
vs2 = TwoDTurb.ForcedVars(g2)
eq2 = TwoDTurb.Equation(pr2, g2)
ts2 = FourierFlows.FilteredForwardEulerTimeStepper(dt, eq2.LC, g2)
prob2 = FourierFlows.Problem(g2, vs2, pr2, eq2, ts2)
s2 = prob2.state


TwoDTurb.set_q!(prob1, ones(g1.X))
TwoDTurb.set_q!(prob2, ones(g1.X))

# I found that calcN_advection! makes LC=0, then @. N += v.F makes LC=v.F
# WHY IS THAT ??

println("INITIALLY:")
println("sol1[2, 3] = ", prob1.state.sol[2, 3])
println("eq1.LC[2, 3] = ", prob1.eqn.LC[2, 3])
println("sol2[2, 3] = ", prob2.state.sol[2, 3])
println("eq2.LC[2, 3] = ", prob2.eqn.LC[2, 3])

println(" ")



println("starting")

println("eq1.LC[2, 3] = ", prob1.eqn.LC[2, 3])
eq1.calcN!(ts1.N, s1.sol, s1.t, s1, vs1, pr1, g1)
println("eq1.LC[2, 3] = ", prob1.eqn.LC[2, 3])
@. s1.sol = s1.sol + ts1.dt*(ts1.N + eq1.LC*s1.sol)
println("eq1.LC[2, 3] = ", prob1.eqn.LC[2, 3])
s1.t += ts1.dt
println("eq1.LC[2, 3] = ", prob1.eqn.LC[2, 3])
s1.step += 1
println("eq1.LC[2, 3] = ", prob1.eqn.LC[2, 3])

println(" - - - - - ")

println("eq2.LC[2, 3] = ", prob2.eqn.LC[2, 3])
eq2.calcN!(ts2.N, s2.sol, s2.t, s2, vs2, pr2, g2)
println("eq2.LC[2, 3] = ", prob2.eqn.LC[2, 3])
# @. s2.sol = ts2.filter*(s2.sol + ts2.dt*(ts2.N + eq2.LC*s2.sol))
@. s2.sol = s2.sol + ts2.dt*(ts2.N + eq2.LC*s2.sol)
println("eq2.LC[2, 3] = ", prob2.eqn.LC[2, 3])
s2.t += ts2.dt
println("eq2.LC[2, 3] = ", prob2.eqn.LC[2, 3])
s2.step += 1
println("eq2.LC[2, 3] = ", prob2.eqn.LC[2, 3])

println(" - - - - - ")

stepforward!(prob1, 1)
stepforward!(prob2, 1)

println("AFTER ONE TIME STEP:")
println("sol1[2, 3] = ", prob1.state.sol[2, 3])
println("eq1.LC[2, 3] = ", prob1.eqn.LC[2, 3])
println("sol2[2, 3] = ", prob2.state.sol[2, 3])
println("eq2.LC[2, 3] = ", prob2.eqn.LC[2, 3])

println(" ")

stepforward!(prob1, 1)
stepforward!(prob2, 1)

println("AFTER ANOTHER TIME STEP:")
println("sol1[2, 3] = ", prob1.state.sol[2, 3])
println("eq1.LC[2, 3] = ", prob1.eqn.LC[2, 3])
println("sol2[2, 3] = ", prob2.state.sol[2, 3])
println("eq2.LC[2, 3] = ", prob2.eqn.LC[2, 3])
