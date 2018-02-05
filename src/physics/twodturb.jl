__precompile__()
module TwoDTurb

using FourierFlows

import FourierFlows: parsevalsum, parsevalsum2

export updatevars!, set_q!, energy, enstrophy, dissipation, work, drag
       

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
     calcF = nothing,
         T = typeof(Lx),
        c0 = nothing,
  )

  g  = TwoDGrid(T, nx, Lx, ny, Ly)

  if calcF == nothing
    pr = TwoDTurb.Params(ν, nν, μ, nμ)
    vs = TwoDTurb.Vars(g)
  elseif c0 != nothing
    pr = TwoDTurb.WithTracerForcedParams(ν, nν, μ, nμ, calcF)
    vs = TwoDTurb.WithTracerForcedVars(g)
  else
    pr = TwoDTurb.ForcedParams(ν, nν, μ, nμ, calcF)
    vs = TwoDTurb.ForcedVars(g)
  end
  eq = TwoDTurb.Equation(pr, g)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LC, g, Complex{T})

  FourierFlows.Problem(g, vs, pr, eq, ts)
end

"""
    InitialValueProblem(; parameters...)

Construct an initial-value 2D turbulence problem.
"""
InitialValueProblem(; kwargs...) = Problem(; kwargs...)

"""
    ForcedProblem(; parameters...)

Construct a forced 2D turbulence problem.
"""
ForcedProblem(; kwargs...) = Problem(; kwargs...)

"""
    Params(ν, nν, μ, nμ)

Returns the params for unforced two-dimensional turbulence.
"""
struct Params <: AbstractParams
  ν::Float64  # Vorticity viscosity
  nν::Int     # Vorticity hyperviscous order
  μ::Float64  # Bottom drag or hypoviscosity
  nμ::Int     # Order of hypodrag
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
  nμ::Int           # Order of hypodrag
  calcF!::Function  # Function that calculates the forcing F
end

"""
    WithTracerForcedParams(ν, nν, μ, nμ, calcF!)

Returns the params for forced two-dimensional turbulence with
hyperviscosity ν and μ of order nν and nμ and forcing calculated by
calcF!.
"""
struct WithTracerForcedParams <: AbstractParams
  ν::Float64        # Vorticity viscosity
  nν::Int           # Vorticity hyperviscous order
  μ::Float64        # Bottom drag or hypoviscosity
  nμ::Int           # Order of hypodrag
  κ::Float64
  nκ::Int
  calcF!::Function  # Function that calculates the forcing F
end


"""
    Equation(p, g)

Returns the equation for two-dimensional turbulence with params p and grid g.
"""
function Equation(p, g)
  LC = -p.ν*g.KKrsq.^p.nν - p.μ*g.KKrsq.^p.nμ
  LC[1, 1] = 0
  FourierFlows.Equation(LC, calcN_advection!)
end

function Equation(p::ForcedParams, g)
  LC = -p.ν*g.KKrsq.^p.nν - p.μ*g.KKrsq.^p.nμ
  LC[1, 1] = 0
  FourierFlows.Equation{Complex{Float64},2}(LC, calcN_forced!)
end

function Equation(p::WithTracerForcedParams, g)
  LC = zeros(g.nkr, g.nl, 2)
  LC[:, :, 1] = -p.ν*g.KKrsq.^p.nν - p.μ*g.KKrsq.^p.nμ
  LC[:, :, 2] = -p.κ*g.KKrsq.^p.nκ
  LC[1, 1, 1] = 0
  FourierFlows.Equation{Complex{Float64},3}(LC, calcN_forced!)
end


# Properties of Vars types:
physifields = [:q, :U, :V, :psi]
transfields = [ Symbol(fld, :h) for fld in physifields ]
forcefields = [:F]

tracerphysifields = [:q, :U, :V, :psi, :c]
tracertransfields = [ Symbol(fld, :h) for fld in tracerphysifields ]

fieldspecs = cat(1,
  FourierFlows.getfieldspecs(physifields, Array{Float64,2}),
  FourierFlows.getfieldspecs(transfields, Array{Complex{Float64},2}))

forcefieldspecs = cat(1, fieldspecs, 
  FourierFlows.getfieldspecs(forcefields, Array{Complex{Float64},2}))

tracerfieldspecs = cat(1,
  FourierFlows.getfieldspecs(tracerphysifields, Array{Float64,2}),
  FourierFlows.getfieldspecs(tracertransfields, Array{Complex{Float64},2}),
  FourierFlows.getfieldspecs(forcefields, Array{Complex{Float64},2}))

# Define Vars types
eval(FourierFlows.getstructexpr(:Vars, fieldspecs; parent=:AbstractVars))
eval(FourierFlows.getstructexpr(:ForcedVars, forcefieldspecs; parent=:AbstractVars))
eval(FourierFlows.getstructexpr(:WithTracerForcedVars, tracerfieldspecs; parent=:AbstractVars))
  
  
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
  @createarrays Complex{Float64} (g.nkr, g.nl) qh Uh Vh psih F
  ForcedVars(q, U, V, psi, qh, Uh, Vh, psih, F)
end

"""
    WithTracerForcedVars(g)

Returns the vars for unforced two-dimensional turbulence with grid g.
"""
function WithTracerForcedVars(g)
  @createarrays Float64 (g.nx, g.ny) q U V psi c
  @createarrays Complex{Float64} (g.nkr, g.nl) qh Uh Vh psih ch F
  ForcedVars(q, U, V, psi, c, qh, Uh, Vh, psih, ch, F)
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

function calcN_forced!(N::Array{Complex{Float64},2}, 
                sol::Array{Complex{Float64},2}, t::Float64, 
                s::State, v::ForcedVars, p::ForcedParams, g::TwoDGrid)

  calcN_advection!(N, sol, t, s, v, p, g)
  p.calcF!(v.F, sol, t, s, v, p, g)

  @. N += v.F
  nothing
end

function calcN_tracer!(N::Array{Complex{Float64},3}, 
                sol::Array{Complex{Float64},3}, t::Float64, 
                s::State, v::WithTracerForcedVars, p::WithTracerForcedParams, g::TwoDGrid)

  @views calcN_forced!(N[:, :, 1], sol[:, :, 1], t, s, v, p, g)

  @. v.Uh =  im * g.l  * g.invKKrsq * sol[:, :, 1]
  @. v.Vh = -im * g.kr * g.invKKrsq * sol[:, :, 1]

  v.ch .= sol[:, :, 2]
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)
  A_mul_B!(v.c, g.irfftplan, v.ch)

  @. v.U *= v.c # U*q
  @. v.V *= v.c # V*q

  A_mul_B!(v.Uh, g.rfftplan, v.U) # \hat{U*q}
  A_mul_B!(v.Vh, g.rfftplan, v.V) # \hat{U*q}

  @views @. N[:, :, 2] = -im*g.kr*v.Uh - im*g.l*v.Vh

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

# Helper functions
"""
    updatevars!(v, s, g)

Update the vars in v on the grid g with the solution in s.sol.
"""
function updatevars!(v::WithTracerForcedVars, s, g)
  v.qh .= s.sol[:, :, 1]
  v.ch .= s.sol[:, :, 2]
  @. v.psih = -g.invKKrsq * v.qh
  @. v.Uh = -im * g.l  * v.psih
  @. v.Vh =  im * g.kr * v.psih

  ch1 = deepcopy(v.ch)
  qh1 = deepcopy(v.qh)
  Uh1 = deepcopy(v.Uh)
  Vh1 = deepcopy(v.Vh)

  A_mul_B!(v.q, g.irfftplan, qh1)
  A_mul_B!(v.c, g.irfftplan, ch1)
  A_mul_B!(v.U, g.irfftplan, Uh1)
  A_mul_B!(v.V, g.irfftplan, Vh1)
  nothing
end

"""
    updatevars!(v, s, g)

Update the vars in v on the grid g with the solution in s.sol.
"""
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
    set_q!(s, v, g, q)

Set the solution s.sol as the transform of q and update variables v 
on the grid g.
"""
function set_q!(s, v::WithTracerForcedVars, g, q)
  @views A_mul_B!(s.sol[:, :, 1], g.rfftplan, q)
  s.sol[1, 1, 1] = 0 # zero out domain average
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
@inline function energy(qh, v, g)
  @. v.Uh = g.invKKrsq * abs2(qh) # qh*Psih
  1/(2*g.Lx*g.Ly)*parsevalsum(v.Uh, g)
end
@inline energy(prob) = energy(prob.state, prob.vars, prob.grid)
@inline energy(s::State, v, g) = energy(s.sol, v, g)
@inline energy(s::State, v::WithTracerForcedParams, g) = energy(s.sol[:, :, 1], v, g)

"""
    enstrophy(prob)

Returns the domain-averaged enstrophy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function enstrophy(qh, g)
  1/(2*g.Lx*g.Ly)*parsevalsum2(qh, g)
end
@inline enstrophy(prob) = enstrophy(prob.state, prob.grid)
@inline enstrophy(s::State, g) = enstrophy(s.sol, g)

"""
    dissipation(prob)

Returns the domain-averaged dissipation rate. nν must be >= 1.
"""
@inline function dissipation(qh, v, p, g)
  @. v.Uh = g.KKrsq^(p.nν-1) * abs2(qh)
  @. v.Uh[1, 1] = 0
  p.ν/(g.Lx*g.Ly)*parsevalsum(v.Uh, g)
end
@inline dissipation(prob) = dissipation(prob.state, prob.vars, prob.params, prob.grid)
@inline dissipation(s::State, v, p, g) = dissipation(s.sol, v, p, g)

"""
    work(prob)

Returns the domain-averaged work done on the fluid by the forcing F 
defined by W = -<ψFᵢ> = -Σ ψh conj(F), where <∙> is the domain average,
Σ denotes a sum over all Fourier modes, ψ is the streamfunction and 
Fᵢ is the inverse transform of F.
"""
@inline function work(s, v::ForcedVars, g)
  @. v.Uh = g.invKKrsq * s.sol * conj(v.F)
  1/(g.Lx*g.Ly)*parsevalsum(v.Uh, g)
end
@inline work(prob) = work(prob.state, prob.vars, prob.grid)

"""
    drag(prob)

Returns the extraction of domain-averaged energy extraction by the drag μ.
"""
@inline function drag(qh, v, p, g)
  @. v.Uh = g.KKrsq^(p.nμ-1) * abs2(qh)
  @. v.Uh[1, 1] = 0
  p.μ/(g.Lx*g.Ly)*parsevalsum(v.Uh, g)
end
@inline drag(prob) = drag(prob.state, prob.vars, prob.params, prob.grid) 
@inline drag(s::State, v, p, g) = drag(s.sol, v, p, g)


end # module
