module TwoDTurb
using FourierFlows
Grid = TwoDGrid

"""
    InitialValueProblem(; parameters...)

Construct an initial-value 2D turbulence problem.
"""
function InitialValueProblem(;
     nx = 256,
     Lx = 2π,
     ny = nx,
     Ly = Lx,
      ν = 0.0,
     nν = 1,
     dt = 0.01,
stepper = "RK4"
  )

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = TwoDTurb.Params(ν, nν)
  vs = TwoDTurb.Vars(g)
  eq = TwoDTurb.Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end

"""
    ForcedProblem(; parameters...)

Construct a forced 2D turbulence problem.
"""
function ForcedProblem(;
        nx = 256,
        Lx = 2π,
        ny = nx,
        Ly = Lx,
         ν = 0.0,
        nν = 1,
         μ = 0.0,
        nμ = 0,
        dt = 0.01,
   stepper = "RK4",
     calcF = nothing,
   force2k = force2k
  )

  if calcF == nothing; _calcF(F, sol, t, s, v, p, g) = nothing
  else;                _calcF = calcF
  end

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = TwoDTurb.ForcedParams(ν, nν, μ, nμ, _calcF, force2k)
  vs = TwoDTurb.ForcedVars(g)
  eq = TwoDTurb.Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end


"""
    Params(ν, nν)

Returns the params for unforced two-dimensional turbulence.
"""
struct Params <: AbstractParams
  ν::Float64        # Vorticity viscosity
  nν::Int           # Vorticity hyperviscous order
end

"""
  ForcedParams(ν, nν, μ, nμ, calcF!)

Returns the params for forced two-dimensional turbulence with
hyperviscosity ν and μ of order nν and nμ and forcing calculated by
calcF!.
"""
struct ForcedParams <: AbstractParams
  ν::Float64        # Vorticity viscosity
  nν::Int           # Vorticity hyperviscous order
  μ::Float64        # Bottom drag or hypoviscosity
  nμ::Float64       # Order of hypodrag
  calcF!::Function  # Function that calculates the forcing F
  force2k::Array{Complex{Float64}, 2}  # the spectrum of the spatial
                                      # correlation of the forcing
end


"""
    Equation(p, g)

Returns the equation for two-dimensional turbulence with params p and grid g.
"""
function Equation(p::Params, g::TwoDGrid)
  LC = -p.ν * g.KKrsq.^p.nν
  FourierFlows.Equation{2}(LC, calcN_advection!)
end

function Equation(p::ForcedParams, g::TwoDGrid)
  LC = -p.ν*g.KKrsq.^p.nν - p.μ*g.KKrsq.^p.nμ
  LC[1, 1] = 0
  FourierFlows.Equation{2}(LC, calcN_forced!)
end


# Construct Vars type for unforced two-dimensional turbulence
physvars = [:q, :U, :V, :psi]
transvars = [:qh, :Uh, :Vh, :psih]
expr = FourierFlows.structvarsexpr(:Vars, physvars, transvars)
eval(expr)

"""
    Vars(g)

Returns the vars for unforced two-dimensional turbulence with grid g.
"""
function Vars(g::TwoDGrid)
  @createarrays Float64 (g.nx, g.ny) q U V psi
  @createarrays Complex{Float64} (g.nkr, g.nl) sol qh Uh Vh psih
  Vars(q, U, V, psi, qh, Uh, Vh, psih)
end

# Construct Vars type for forced two-dimensional turbulence
forcedtransvars = [:qh, :Uh, :Vh, :psih, :F]
expr = FourierFlows.structvarsexpr(:ForcedVars, physvars, forcedtransvars)
eval(expr)


"""
    ForcedVars(g)

Returns the vars for unforced two-dimensional turbulence with grid g.
"""
function ForcedVars(g::TwoDGrid)
  @createarrays Float64 (g.nx, g.ny) q U V psi
  @createarrays Complex{Float64} (g.nkr, g.nl) sol qh Uh Vh psih F
  ForcedVars(q, U, V, psi, qh, Uh, Vh, psih, F)
end


# Solvers
function calcN_advection!(
  N::Array{Complex{Float64},2}, sol::Array{Complex{Float64},2},
  t::Float64, s::State, v::AbstractVars, p::AbstractParams, g::TwoDGrid)

  @. v.Uh =  im * g.l  * g.invKKrsq * sol
  @. v.Vh = -im * g.kr * g.invKKrsq * sol

  v.qh .= sol
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)
  A_mul_B!(v.q, g.irfftplan, v.qh)

  @. v.U *= v.q # U*q
  @. v.V *= v.q # V*q

  A_mul_B!(v.Uh, g.rfftplan, v.U) # \hat{U*q}
  A_mul_B!(v.Vh, g.rfftplan, v.V) # \hat{U*q}

  @. N = -im*g.kr*v.Uh - im*g.l*v.Vh
  nothing
end

function calcN_forced!(N::Array{Complex{Float64}, 2},
                sol::Array{Complex{Float64}, 2}, t::Float64,
                s::State, v::ForcedVars, p::ForcedParams, g::TwoDGrid)

  calcN_advection!(N, sol, t, s, v, p, g)
  p.calcF!(v.F, sol, t, s, v, p, g)

  @. N += v.F
  nothing
end


# Helper functions
"""
    updatevars!(v, s, g)

Update the vars in v on the grid g with the solution in s.sol.
"""
function updatevars!(v, s, g)
  v.qh .= s.sol
  @. v.psih = -g.invKKrsq * v.qh
  @. v.Uh = -im * g.l  * v.psih
  @. v.Vh =  im * g.kr * v.psih

  qh1 = deepcopy(v.qh)
  Uh1 = deepcopy(v.Uh)
  Vh1 = deepcopy(v.Vh)

  A_mul_B!(v.q, g.irfftplan, qh1)
  A_mul_B!(v.U, g.irfftplan, Uh1)
  A_mul_B!(v.V, g.irfftplan, Vh1)
  nothing
end

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.vars, prob.state, prob.grid)
end

"""
    set_q!(s, v, g, q)

Set the solution s.sol as the transform of q and update variables v
on the grid g.
"""
function set_q!(s, v, g, q)
  A_mul_B!(s.sol, g.rfftplan, q)
  s.sol[1, 1] = 0 # zero out domain average
  updatevars!(v, s, g)
end

"""
    set_q!(prob, q)

Set the solution prob.state.sol as the transform of q and update variables.
"""
function set_q!(prob::AbstractProblem, q)
  set_q!(prob.state, prob.vars, prob.grid, q)
end

"""
    energy(s, v, g)

Returns the domain-averaged kinetic energy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function energy(s, v, g)
  @. v.Uh = g.invKKrsq * abs2(s.sol)
  1/(2*g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

@inline function energy(prob)
  energy(prob.state, prob.vars, prob.grid)
end

"""
    enstrophy(s, g)

Returns the domain-averaged enstrophy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function enstrophy(s, g)
  1/(2*g.Lx*g.Ly)*FourierFlows.parsevalsum2(s.sol, g)
end

@inline function enstrophy(prob)
  enstrophy(prob.state, prob.grid)
end


"""
    dissipation(s, v, p, g)

Returns the domain-averaged dissipation rate. nν must be >= 1.
"""
@inline function dissipation(s, v, p, g)
  @. v.Uh = g.KKrsq^(p.nν-1) * abs2(s.sol)
  @. v.Uh[1, 1] = 0
  p.ν/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

@inline function dissipation(prob::AbstractProblem)
  dissipation(prob.state, prob.vars, prob.params, prob.grid)
end

"""
    injection(s, v, p, g)

Returns the domain-averaged rate of injection of energy by the forcing F.
"""
@inline function injection(s, v::ForcedVars, g)
  @. v.Uh = g.invKKrsq * s.sol * conj(v.F)
  1/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

@inline function injection(prob::AbstractProblem)
  injection(prob.state, prob.vars, prob.grid)
end

"""
    drag(s, v, p, g)

Returns the extraction of domain-averaged energy by drag μ.
"""
@inline function drag(s, v::ForcedVars, p::ForcedParams, g)
  @. v.Uh = g.KKrsq^(p.nμ-1) * abs2(s.sol)
  @. v.Uh[1, 1] = 0
  p.μ/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

@inline function drag(prob::AbstractProblem)
  drag(prob.state, prob.vars, prob.params, prob.grid)
end

end # module
