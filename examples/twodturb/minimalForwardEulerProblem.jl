using PyPlot, FourierFlows
import FourierFlows.TwoDTurb
import FourierFlows.TwoDTurb: energy, enstrophy, dissipation, injection, drag

 n, L  =  128, 2π
 ν, nν = 0e-3, 1
 μ, nμ = 0e-1, 0
dt = 0.1


gr  = TwoDGrid(n, L)

println(" ")
println("Forcing ONLY component [2, 3] with constant forcing F[2, 3] = 2 + i")
println(" ")

function calcF!(F, sol, t, s, v, p, g)
  if t == s.t # not a substep
    F .= 0.0
    F[2, 3] .= 2.0+im
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
ts1 = FourierFlows.autoconstructtimestepper("ForwardEuler", dt, eq1.LC, g1)
prob1 = FourierFlows.Problem(g1, vs1, pr1, eq1, ts1)

g2  = TwoDGrid(nx, Lx, ny, Ly)
pr2 = TwoDTurb.ForcedParams(ν, nν, μ, nμ, _calcF)
vs2 = TwoDTurb.ForcedVars(g2)
eq2 = TwoDTurb.Equation(pr2, g2)
ts2 = FourierFlows.autoconstructtimestepper("FilteredForwardEuler", dt, eq2.LC, g2)
prob2 = FourierFlows.Problem(g2, vs2, pr2, eq2, ts2)


TwoDTurb.set_q!(prob1, 0*g1.X)
TwoDTurb.set_q!(prob2, 0*g2.X)


println("INITIALLY:")
println("sol1[2, 3] = ", prob1.state.sol[2, 3])
println("eq1.LC[2, 3] = ", prob1.eqn.LC[2, 3])
println("sol2[2, 3] = ", prob2.state.sol[2, 3])
println("eq2.LC[2, 3] = ", prob2.eqn.LC[2, 3])

println(" ")

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
