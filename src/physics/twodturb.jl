__precompile__()
module TwoDTurb
using FourierFlows

"""
    Problem(; parameters...)

Construct a forced 2D turbulence problem.
"""
function Problem(;
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
     calcF = nothing
  )

  dt = Float64(dt)
  g  = TwoDGrid(nx, Lx, ny, Ly)

  if calcF == nothing
    pr = TwoDTurb.Params(ν, nν, μ, nμ)
    vs = TwoDTurb.Vars(g)
  else
    pr = TwoDTurb.ForcedParams(ν, nν, μ, nμ, _calcF)
    vs = TwoDTurb.ForcedVars(g)
  end
  eq = TwoDTurb.Equation(pr, g)
  st = FourierFlows.State(Complex{Float64}, (g.nkr, g.nl), dt)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, st.sol, g)

  FourierFlows.Problem(g, vs, pr, eq, ts, st)
end

"""
    InitialValueProblem(; parameters...)

Construct an initial-value 2D turbulence problem.
"""
function InitialValueProblem(; kwargs...)
  Problem(; kwargs...)
end

"""
    ForcedProblem(; parameters...)

Construct a forced 2D turbulence problem.
"""
function ForcedProblem(; kwargs...)
  Problem(; kwargs...)
end

"""
    Params(ν, nν)

Returns the params for unforced two-dimensional turbulence.
"""
struct Params <: AbstractParams
  ν::Float64        # Vorticity viscosity
  nν::Int           # Vorticity hyperviscous order
  μ::Float64        # Bottom drag or hypoviscosity
  nμ::Float64       # Order of hypodrag
end

Params(ν, nν) = Params(ν, nν, 0, 0)

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
end

Params(ν, nν, μ, nμ, calcF) = ForcedParams(ν, nν, μ, nμ, calcF)

"""
    Equation(p, g)

Returns the equation for two-dimensional turbulence with params p and grid g.
"""
function Equation(p, g)
  LC = -p.ν*g.KKrsq.^p.nν - p.μ*g.KKrsq.^p.nμ
  LC[1, 1] = 0
  FourierFlows.Equation{Complex{Float64},2}(LC, calcN_advection!)
end

function Equation(p::ForcedParams, g)
  eq = Equation(p, g) 
  FourierFlows.Equation{Complex{Float64},2}(eq.LC, calcN_forced!)
end

# Properties of Vars types:
physifields = [:q, :U, :V, :psi]
transfields = [ Symbol(fld, :h) for fld in physifields ]
forcefields = [:F, :F₋₁, :sol₋₁]

fieldspecs = cat(1,
  FourierFlows.getfieldspecs(physifields, Array{Float64,2}),
  FourierFlows.getfieldspecs(transfields, Array{Complex{Float64},2}))

forcefieldspecs = cat(1, fieldspecs, 
  FourierFlows.getfieldspecs(forcefields, Array{Complex{Float64},2}))

# Define Vars type for unforced two-dimensional turbulence
eval(FourierFlows.getstructexpr(:Vars, fieldspecs; parent=:AbstractVars))

# Define Vars type for forced two-dimensional turbulence
eval(FourierFlows.getstructexpr(:ForcedVars, forcefieldspecs; 
  parent=:AbstractVars))
  
"""
    Vars(g)

Returns the vars for unforced two-dimensional turbulence with grid g.
"""
function Vars(g)
  @createarrays Float64 (g.nx, g.ny) q U V psi
  @createarrays Complex{Float64} (g.nkr, g.nl) qh Uh Vh psih
  Vars(q, U, V, psi, qh, Uh, Vh, psih)
end

"""
    ForcedVars(g)

Returns the vars for unforced two-dimensional turbulence with grid g.
"""
function ForcedVars(g)
  @createarrays Float64 (g.nx, g.ny) q U V psi
  @createarrays Complex{Float64} (g.nkr, g.nl) qh Uh Vh psih F F₋₁ sol₋₁
  ForcedVars(q, U, V, psi, qh, Uh, Vh, psih, F, F₋₁, sol₋₁)
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

# Diagnostics
"""
    energy(prob)

Returns the domain-averaged kinetic energy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function energy(s, v, g)
  @. v.Uh = g.invKKrsq * abs2(s.sol)
  1/(2*g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

"""
    enstrophy(prob)

Returns the domain-averaged enstrophy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function enstrophy(s, g)
  1/(2*g.Lx*g.Ly)*FourierFlows.parsevalsum2(s.sol, g)
end

"""
    dissipation(prob)

Returns the domain-averaged dissipation rate. nν must be >= 1.
"""
@inline function dissipation(s, v, p, g)
  @. v.Uh = g.KKrsq^(p.nν-1) * abs2(s.sol)
  @. v.Uh[1, 1] = 0
  p.ν/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

"""
    work(prob)

Returns the domain-averaged work done on the fluid by the forcing F 
defined by W = -<ψFᵢ> = -Σ ψh conj(F), where <∙> is the domain average,
Σ denotes a sum over all Fourier modes, ψ is the streamfunction and 
Fᵢ is the inverse transform of F.
"""
function work(s, v::ForcedVars, g)
  @. v.Uh = g.invKKrsq * s.sol * conj(v.F)
  1/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

"""
    drag(prob)

Returns the extraction of domain-averaged energy by drag μ.
"""
function drag(s, v::ForcedVars, p::ForcedParams, g)
  @. v.Uh = g.KKrsq^(p.nμ-1) * abs2(s.sol)
  @. v.Uh[1, 1] = 0
  p.μ/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

@inline      energy(prob) = energy(prob.state, prob.vars, prob.grid)
@inline   enstrophy(prob) = enstrophy(prob.state, prob.grid)
@inline        drag(prob) = drag(prob.state, prob.vars, prob.params, prob.grid) 
@inline        work(prob) = work(prob.state, prob.vars, prob.grid)
@inline dissipation(prob) = dissipation(prob.state, prob.vars, prob.params, 
                                        prob.grid)

end # module
