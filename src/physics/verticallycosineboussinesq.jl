__precompile__()
module VerticallyCosineBoussinesq
using FourierFlows

import FourierFlows: autoconstructtimestepper, parsevalsum, parsevalsum2

export updatevars!, set_Z!, set_uvp!, set_planewave!, 
       totalenergy, mode0energy, mode0enstrophy, mode1ke, mode1pe, mode1energy,
       mode0dissipation, mode1dissipation, mode0drag, mode1drag

""" 
    Problem(; parameters...)

Construct a VerticallyCosineBoussinesq problem.
"""
function Problem(;
  # Numerical parameters
  nx = 128, 
  Lx = 2π, 
  ny = nx,
  Ly = Lx,
  dt = 0.01, 
  # Drag and/or hyper-/hypo-viscosity
   nu0 = 0, # barotropic viscosity
  nnu0 = 1, 
   nu1 = 0, # baroclinic viscosity
  nnu1 = 1, 
   mu0 = 0, # barotropic drag / hypoviscosity
  nmu0 = 0, 
   mu1 = 0, # baroclinic drag / hypoviscosity
  nmu1 = 0,
  # Physical parameters
  f = 1,
  N = 1,
  m = 1,
  # Optional uniform and steady background flow
  Ub = 0,
  Vb = 0,
  # Timestepper and eqn options
  stepper = "RK4",
  linear = false,
  linearized = false,
  calcF = nothing,
  )

  g = TwoDGrid(nx, Lx, ny, Ly)

  if calcF != nothing
    pr = ForcedParams(nu0, nnu0, nu1, nnu1, mu0, nmu0, mu1, nmu1, f, N, m,
      calcF; Ub=Ub, Vb=Vb)
    vs = ForcedVars(g)
  else
    pr = Params(nu0, nnu0, nu1, nnu1, mu0, nmu0, mu1, nmu1, f, N, m; 
      Ub=Ub, Vb=Vb)
    vs = Vars(g)
  end

  if calcF != nothing; eq = ForcedEquation(pr, g)
  elseif linear;       eq = LinearEquation(pr, g)
  elseif linearized;   eq = LinearizedEquation(pr, g)
  else;                eq = Equation(pr, g)
  end

  ts = autoconstructtimestepper(stepper, dt, eq.LC, g)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end

# Params
abstract type VerticallyCosineParams <: AbstractParams end

"""
    Params(nu0, nnu0, nu1, nnu1, mu0, nmu0, mu1, nmu1, f, N, m; 
           Ub=0, Vb=0) 

Construct parameters for the Two-Fourier-mode Boussinesq problem. Suffix 0
refers to zeroth mode; 1 to first mode. f, N, m are Coriolis frequency, 
buoyancy frequency, and vertical wavenumber of the first mode, respectively.
The optional constant background velocity (Ub,Vb) is set to zero by default.
The viscosity is applied only to the first-mode horizontal velocities.
"""
struct Params <: VerticallyCosineParams
  nu0::Float64    # Mode-0 viscosity
  nnu0::Int       # Mode-0 hyperviscous order
  nu1::Float64    # Mode-1 viscosity
  nnu1::Int       # Mode-1 hyperviscous order
  mu0::Float64     # Hypoviscosity/bottom drag 
  nmu0::Float64    # Order of hypoviscosity (nmu=0 for bottom drag)
  mu1::Float64     # Hypoviscosity/bottom drag 
  nmu1::Float64    # Order of hypoviscosity (nmu=0 for bottom drag)
  f::Float64      # Planetary vorticity
  N::Float64      # Buoyancy frequency
  m::Float64      # Mode-one wavenumber
  Ub::Float64     # Steady background barotropic x-velocity
  Vb::Float64     # Steady background barotropic y-velocity
end

Params(nu0, nnu0, nu1, nnu1, mu0, nmu0, mu1, nmu1, f, N, m; Ub=0, 
  Vb=0) = Params(nu0, nnu0, nu1, nnu1, mu0, nmu0, mu1, nmu1, f, N, m, Ub, Vb)

struct ForcedParams <: VerticallyCosineParams
  nu0::Float64    # Mode-0 viscosity
  nnu0::Int       # Mode-0 hyperviscous order
  nu1::Float64    # Mode-1 viscosity
  nnu1::Int       # Mode-1 hyperviscous order
  mu0::Float64     # Hypoviscosity/bottom drag 
  nmu0::Float64    # Order of hypoviscosity (nmu=0 for bottom drag)
  mu1::Float64     # Hypoviscosity/bottom drag 
  nmu1::Float64    # Order of hypoviscosity (nmu=0 for bottom drag)
  f::Float64      # Planetary vorticity
  N::Float64      # Buoyancy frequency
  m::Float64      # Mode-one wavenumber
  Ub::Float64     # Steady background barotropic x-velocity
  Vb::Float64     # Steady background barotropic y-velocity
  calcF!::Function
end

ForcedParams(nu0, nnu0, nu1, nnu1, mu0, nmu0, mu1, nmu1, f, N, m,
   calcF::Function; Ub=0, Vb=0) = ForcedParams(nu0, nnu0, nu1, nnu1, mu0, nmu0, 
    mu1, nmu1, f, N, m, calcF, Ub, Vb)

# Equations
function Equation(p::VerticallyCosineParams, g::TwoDGrid)
  LC = zeros(Complex{Float64}, g.nkr, g.nl, 4)
  @views @. LC[:, :, 1] = -p.nu0*g.KKrsq^p.nnu0 - p.mu0*g.KKrsq^p.nmu0
  @views @. LC[:, :, 2] = -p.nu1*g.KKrsq^p.nnu1 - p.mu1*g.KKrsq^p.nmu1
  @views @. LC[:, :, 3] = -p.nu1*g.KKrsq^p.nnu1 - p.mu1*g.KKrsq^p.nmu1
  FourierFlows.Equation(LC, calcN!)
end

function LinearizedEquation(p, g)
  eq = Equation(p, g)
  FourierFlows.Equation(eq.LC, calcN_linearized!)
end

function LinearEquation(p, g)
  eq = Equation(p, g)
  FourierFlows.Equation(eq.LC, calcN_linearterms!)
end

function ForcedEquation(p, g)
  eq = Equation(p, g)
  FourierFlows.Equation(eq.LC, calcN_forced!)
end


# Vars
abstract type VerticallyCosineVars <: AbstractVars end

physifields = [:Z, :U, :V, :UZuz, :VZvz, :Ux, :Uy, :Vx, :Vy, :Psi, 
  :u, :v, :w, :p, :zeta, :Uu, :Uv, :Up, :Vu, :Vv, :Vp, :uUxvUy, :uVxvVy]
transfields = [ Symbol(var, :h) for var in physifields ]
forcefields = [:F]

fieldspecs = cat(1,
  FourierFlows.getfieldspecs(physifields, Array{Float64,2}),
  FourierFlows.getfieldspecs(transfields, Array{Complex{Float64},2}))

forcefieldspecs = cat(1, fieldspecs,
  FourierFlows.getfieldspecs(forcefields, Array{Complex{Float64},3}))

# Define Vars type for unforced problem
eval(
  FourierFlows.getstructexpr(:Vars, fieldspecs; parent=:VerticallyCosineVars))

# Define Vars type for forced problem
eval(FourierFlows.getstructexpr(:ForcedVars, forcefieldspecs; 
  parent=:VerticallyCosineVars))

"""
    Vars(g)

Returns the vars for unforced two-vertical-cosine-mode Boussinesq dynamics
on the grid g.
"""
function Vars(g)
  @createarrays Float64 (g.nx, g.ny) Z U V UZuz VZvz Ux Uy Vx Vy Psi u v w p
  @createarrays Float64 (g.nx, g.ny) zeta Uu Uv Up Vu Vv Vp uUxvUy uVxvVy
  @createarrays Complex{Float64} (g.nkr, g.nl) Zh Uh Vh UZuzh VZvzh Uxh Uyh
  @createarrays Complex{Float64} (g.nkr, g.nl) Vxh Vyh Psih uh vh wh ph zetah
  @createarrays Complex{Float64} (g.nkr, g.nl) Uuh Uvh Uph Vuh Vvh Vph
  @createarrays Complex{Float64} (g.nkr, g.nl) uUxvUyh uVxvVyh

  Vars(
    Z, U, V, UZuz, VZvz, Ux, Uy, Vx, Vy, Psi, 
    u, v, w, p, zeta, Uu, Uv, Up, Vu, Vv, Vp, uUxvUy, uVxvVy,
    Zh, Uh, Vh, UZuzh, VZvzh, Uxh, Uyh, Vxh, Vyh, Psih, 
    uh, vh, wh, ph, zetah, Uuh, Uvh, Uph, Vuh, Vvh, Vph, uUxvUyh, uVxvVyh)
end

"""
    ForcedVars(g)

Returns the vars for forced two-vertical-cosine-mode Boussinesq dynamics
on the grid g.
"""
function ForcedVars(g)
  @createarrays Float64 (g.nx, g.ny) Z U V UZuz VZvz Ux Uy Vx Vy Psi u v w p
  @createarrays Float64 (g.nx, g.ny) zeta Uu Uv Up Vu Vv Vp uUxvUy uVxvVy
  @createarrays Complex{Float64} (g.nkr, g.nl) Zh Uh Vh UZuzh VZvzh Uxh Uyh
  @createarrays Complex{Float64} (g.nkr, g.nl) Vxh Vyh Psih uh vh wh ph zetah
  @createarrays Complex{Float64} (g.nkr, g.nl) Uuh Uvh Uph Vuh Vvh Vph
  @createarrays Complex{Float64} (g.nkr, g.nl) uUxvUyh uVxvVyh
  @createarrays Complex{Float64} (g.nkr, g.nl, 4) F

  ForcedVars(
    Z, U, V, UZuz, VZvz, Ux, Uy, Vx, Vy, Psi, 
    u, v, w, p, zeta, Uu, Uv, Up, Vu, Vv, Vp, uUxvUy, uVxvVy,
    Zh, Uh, Vh, UZuzh, VZvzh, Uxh, Uyh, Vxh, Vyh, Psih, 
    uh, vh, wh, ph, zetah, Uuh, Uvh, Uph, Vuh, Vvh, Vph, uUxvUyh, uVxvVyh, 
    F)
end



# Solvers
function calcN_linearterms!(
  N::Array{Complex{Float64},3}, sol::Array{Complex{Float64},3},
  t::Float64, s::State, v::Vars, p::VerticallyCosineParams, g::TwoDGrid)

  # Zeroth-mode:
  @views @. N[:, :, 1] = 0
  # First-mode linear terms:
  # u
  @views @. N[:, :, 2] =  p.f*sol[:, :, 3] - im*g.kr*sol[:, :, 4]
  # v
  @views @. N[:, :, 3] = -p.f*sol[:, :, 2] - im*g.l *sol[:, :, 4]
  # p
  @views @. N[:, :, 4] = -im*p.N^2/p.m^2*(g.kr*sol[:, :, 2] + g.l*sol[:, :, 3])

  nothing
end


function calcN!(
  N::Array{Complex{Float64},3}, sol::Array{Complex{Float64},3},
  t::Float64, s::State, v::Vars, p::VerticallyCosineParams, g::TwoDGrid)

  @views v.Zh .= sol[:, :, 1]
  @views v.uh .= sol[:, :, 2]
  @views v.vh .= sol[:, :, 3]
  @views v.ph .= sol[:, :, 4]

  @. v.wh = -im/p.m*(g.kr*v.uh + g.l*v.vh)

  # Spectral-space calculations
  @. v.Psih = -g.invKKrsq*v.Zh

  @. v.Uh = -im*g.l*v.Psih
  @. v.Vh =  im*g.kr*v.Psih

  @. v.Uxh = im*g.kr*v.Uh
  @. v.Vxh = im*g.kr*v.Vh

  @. v.Uyh = im*g.l*v.Uh
  @. v.Vyh = im*g.l*v.Vh

  @. v.zetah = im*g.kr*v.vh - im*g.l*v.uh

  v.Uh[1, 1] += p.Ub*g.nx*g.ny
  v.Vh[1, 1] += p.Vb*g.nx*g.ny

  # Inverse transforms
  A_mul_B!(v.u, g.irfftplan, v.uh)
  A_mul_B!(v.v, g.irfftplan, v.vh)
  A_mul_B!(v.p, g.irfftplan, v.ph)
  A_mul_B!(v.w, g.irfftplan, v.wh)
  A_mul_B!(v.zeta, g.irfftplan, v.zetah)

  A_mul_B!(v.Z, g.irfftplan, v.Zh)
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)

  # Multiplies
  @. v.UZuz = v.U*v.Z + 0.5*v.u*v.zeta - 0.5*p.m*v.v*v.w 
  @. v.VZvz = v.V*v.Z + 0.5*v.v*v.zeta + 0.5*p.m*v.u*v.w 

  @. v.Uu = v.U * v.u
  @. v.Vu = v.V * v.u
  @. v.Uv = v.U * v.v
  @. v.Vv = v.V * v.v
  @. v.Up = v.U * v.p
  @. v.Vp = v.V * v.p

  @. v.uUxvUy = v.u*v.Ux + v.v*v.Uy
  @. v.uVxvVy = v.u*v.Vx + v.v*v.Vy

  # Forward transforms
  A_mul_B!(v.UZuzh, g.rfftplan, v.UZuz)
  A_mul_B!(v.VZvzh, g.rfftplan, v.VZvz)

  A_mul_B!(v.Uuh, g.rfftplan, v.Uu)
  A_mul_B!(v.Uvh, g.rfftplan, v.Uv)
  A_mul_B!(v.Vuh, g.rfftplan, v.Vu)
  A_mul_B!(v.Vvh, g.rfftplan, v.Vv)
  A_mul_B!(v.Uph, g.rfftplan, v.Up)
  A_mul_B!(v.Vph, g.rfftplan, v.Vp)

  A_mul_B!(v.uUxvUyh, g.rfftplan, v.uUxvUy)
  A_mul_B!(v.uVxvVyh, g.rfftplan, v.uVxvVy)

  # Linear terms
  calcN_linearterms!(N, sol, t, s, v, p, g)

  # Zeroth-mode nonlinear term
  @views @. N[:, :, 1] = - im*g.kr*v.UZuzh - im*g.l*v.VZvzh
  # First-mode nonlinear terms:
  # u
  @views @. N[:, :, 2] += -im*g.kr*v.Uuh - im*g.l*v.Vuh - v.uUxvUyh
  # v
  @views @. N[:, :, 3] += -im*g.kr*v.Uvh - im*g.l*v.Vvh - v.uVxvVyh
  # p
  @views @. N[:, :, 4] += -im*g.kr*v.Uph - im*g.l*v.Vph

  nothing
end

function calcN_linearized!(N, sol, t, s, v, p, g)
  calcN!(N, sol, t, s, v, p, g)

  # "Linearize" the vorticity term by recalculating
  @. v.UZuz = v.U*v.Z
  @. v.VZvz = v.V*v.Z
  A_mul_B!(v.UZuzh, g.rfftplan, v.UZuz)
  A_mul_B!(v.VZvzh, g.rfftplan, v.VZvz)

  @views @. N[:, :, 1] = - im*g.kr*v.UZuzh - im*g.l*v.VZvzh
  nothing
end

function calcN_forced!(N, sol, t, s, v, p, g)
  calcN!(N, sol, t, s, v, p, g)
  p.calcF!(v.F, sol, t, s, v, p, g)
  @. N += v.F
  nothing
end



# Helper functions
"""
    updatevars!(prob)

Update variables to correspond to the solution in s.sol or prob.state.sol.
"""

function updatevars!(v, s, p, g)
  @views @. v.Zh = s.sol[:, :, 1]
  @views @. v.uh = s.sol[:, :, 2]
  @views @. v.vh = s.sol[:, :, 3]
  @views @. v.ph = s.sol[:, :, 4]

  @. v.Psih = -g.invKKrsq*v.Zh
  @. v.Uh   = -im*g.l*v.Psih
  @. v.Vh   =  im*g.kr*v.Psih
  @. v.wh   = -im/p.m*(g.kr*v.uh + g.l*v.vh)

  Psih = deepcopy(v.Psih)
  Uh = deepcopy(v.Uh)
  Vh = deepcopy(v.Vh)
  Zh = deepcopy(v.Zh)

  uh = deepcopy(v.uh)
  vh = deepcopy(v.vh)
  ph = deepcopy(v.ph)
  wh = deepcopy(v.wh)

  A_mul_B!(v.Psi, g.irfftplan, Psih)
  A_mul_B!(v.U, g.irfftplan, Uh)
  A_mul_B!(v.V, g.irfftplan, Vh)
  A_mul_B!(v.Z, g.irfftplan, Zh)
  A_mul_B!(v.u, g.irfftplan, uh)
  A_mul_B!(v.v, g.irfftplan, vh)
  A_mul_B!(v.p, g.irfftplan, ph)
  A_mul_B!(v.w, g.irfftplan, wh)
  nothing
end
updatevars!(prob) = updatevars!(prob.vars, prob.state, prob.params, prob.grid)

"""
    set_Z!(prob, Z)

Set zeroth mode vorticity and update vars. 
"""
function set_Z!(s, v, p, g, Z)
  @views A_mul_B!(s.sol[:, :, 1], g.rfftplan, Z)
  updatevars!(v, s, p, g)
  nothing
end
set_Z!(prob, Z) = set_Z!(prob.state, prob.vars, prob.params, prob.grid, Z)

""" 
    set_uvp!(prob)

Set first mode u, v, and p and update vars.
"""
function set_uvp!(s, vs, pr, g, u, v, p)
  @views A_mul_B!(s.sol[:, :, 2], g.rfftplan, u)
  @views A_mul_B!(s.sol[:, :, 3], g.rfftplan, v)
  @views A_mul_B!(s.sol[:, :, 4], g.rfftplan, p)
  updatevars!(vs, s, pr, g)
  nothing
end
set_uvp!(prob, u, v, p) = set_uvp!(prob.state, prob.vars, prob.params, 
                                   prob.grid, u, v, p)

""" 
    set_planewave!(prob, u₀, κ, θ=0)

Set a plane wave solution with initial speed u₀, non-dimensional wave
number κ, and angle θ with the horizontal. The non-dimensional wavenumber 
vector is (k, l) = (κ cos θ, κ sin θ), is normalized by 2π/Lx, and is rounded
to the nearest integer.
"""
function set_planewave!(s, vs, pr, g, u₀, κ, θ=0; envelope=nothing)
  k = 2π/g.Lx*round(Int, κ*cos(θ))
  l = 2π/g.Lx*round(Int, κ*sin(θ))
  x, y = g.X, g.Y

  # Wave parameters
  f, N, m = pr.f, pr.N, pr.m
  σ = sqrt( f^2 + N^2*(k^2 + l^2)/m^2 )

  v₀ = u₀ * (f*k - σ*l)/(σ*k - f*l)
  p₀ = u₀ * (σ^2 - f^2)/(σ*k - f*l)

  Φ = k*x + l*y
  u = u₀ * cos.(Φ)
  v = v₀ * sin.(Φ)
  p = p₀ * cos.(Φ)

  if envelope != nothing
    @. u *= envelope(x, y)
    @. v *= envelope(x, y)
    @. p *= envelope(x, y)
  end
  
  set_uvp!(s, vs, pr, g, u, v, p)
  nothing
end
set_planewave!(prob, uw, nkw, θ=0; kwargs...) = set_planewave!(
  prob.state, prob.vars, prob.params, prob.grid, uw, nkw, θ; kwargs...)

# Diagnostics
""" 
    mode0energy(prob)

Returns the domain-averaged energy in the zeroth mode.
"""
@inline function mode0energy(s, v, g)
  @views @. v.Uh = g.invKKrsq * abs2(s.sol[:, :, 1])
  1/(2*g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

"""
    mode0enstrophy(prob)

Returns the domain-averaged enstrophy in the Fourier-transformed vorticity
solution s.sol.
"""
@inline function mode0enstrophy(s, g)
  @views 1/(2*g.Lx*g.Ly)*FourierFlows.parsevalsum2(s.sol[:, :, 1], g)
end

"""
    mode0dissipation(prob)

Returns the domain-averaged barotropic dissipation rate. nnu0 must be >= 1.
"""
@inline function mode0dissipation(s, v, p, g)
  @views @. v.Uh = g.KKrsq^(p.nnu-1) * abs2(s.sol[:, :, 1])
  p.nu/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

"""
    mode0drag(prob)

Returns the extraction of domain-averaged barotropic energy by drag μ.
"""
@inline function mode0drag(s, v, p, g)
  @. v.Uh = g.KKrsq^(p.nmu-1) * abs2(s.sol)
  @. v.Uh[1, 1] = 0
  p.mu/(g.Lx*g.Ly)*FourierFlows.parsevalsum(v.Uh, g)
end

"""
    mode1ke(prob)

Returns the domain-averaged kinetic energy in the first mode.
"""
@inline function mode1ke(s, g) 
  @views 1/(4*g.Lx*g.Ly)*(parsevalsum2(s.sol[:, :, 2], g) 
    + parsevalsum2(s.sol[:, :, 3], g))
end 

"""
    mode1pe(prob)

Returns the domain-averaged potential energy in the first mode.
"""
@inline mode1pe(s, p, g) = @views p.m^2/(4*g.Lx*g.Ly*p.N^2)*parsevalsum2(
  s.sol[:, :, 4], g)

"""
    mode1energy(prob)

Returns the domain-averaged total energy in the first mode.
"""
@inline mode1energy(s, p, g) = mode1ke(s, g) + mode1pe(s, p, g)

"""
    mode1dissipation(prob)

Returns the domain-averaged kinetic energy dissipation of the first mode 
by horizontal viscosity.
"""
@inline function mode1dissipation(s, v, p, g)
  @views @. v.Uuh = g.kr^p.nnu1*s.sol[:, :, 2]
  @views @. v.Vuh =  g.l^p.nnu1*s.sol[:, :, 3]
  p.nu1/(2*g.Lx*g.Ly)*(parsevalsum2(v.Uuh, g) + parsevalsum2(v.Vuh, g))    
end

"""
    mode1drag(prob)

Returns the domain-averaged kinetic energy dissipation of the first mode 
by horizontal viscosity.
"""
@inline function mode1drag(s, v, p, g)
  @views @. v.Uuh = g.kr^p.nmu1*s.sol[:, :, 2]
  @views @. v.Vuh =  g.l^p.nmu1*s.sol[:, :, 3]
  if p.nmu1 != 0 # zero out zeroth mode
    @views @. v.Uuh[1, :] = 0 
    @views @. v.Vuh[:, 1] = 0
  end
  p.mu1/(2*g.Lx*g.Ly)*(parsevalsum2(v.Uuh, g) + parsevalsum2(v.Vuh, g))
end
                                                  
"""
    totalenergy(prob)

Returns the total energy projected onto the zeroth mode.
"""
@inline totalenergy(s, v, p, g) = mode0energy(s, v, g) + mode1energy(s, p, g)

mode0energy(pb) = mode0energy(pb.state, pb.vars, pb.grid)
mode0enstrophy(pb) = mode0enstrophy(pb.state, pb.grid)
mode0dissipation(pb) = mode0dissipation(pb.state, pb.vars, pb.params, pb.grid)
mode0drag(pb) = mode0drag(pb.state, pb.vars, pb.params, pb.grid) 
mode1ke(pb) = mode1ke(pb.state, pb.grid)
mode1pe(pb) = mode1pe(pb.state, pb.params, pb.grid)
mode1energy(pb) = mode1energy(pb.state, pb.params, pb.grid)
mode1dissipation(pb) = mode1dissipation(pb.state, pb.vars, pb.params, pb.grid)
mode1drag(pb) = mode1drag(pb.state, pb.vars, pb.params, pb.grid)
totalenergy(pb) = totalenergy(pb.state, pb.vars, pb.params, pb.grid)
                                        
end # module
