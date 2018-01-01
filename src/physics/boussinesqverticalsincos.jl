__precompile__()

# A module for solving a two-vertical-mode truncation of the Boussinesq 
# equations.
module RigidLidBoussinesq

using FourierFlows


""" 
Construct a RigidLidBoussinesq initial value problem.
"""
function InitialValueProblem(;
  # Numerical parameters
  nx   = 128, 
  Lx   = 2π, 
  ny   = nothing,
  Ly   = nothing,
  dt   = 0.01, 
  # Viscosity and damping
  mu   = nothing,
  nu0  = nothing, 
  nnu0 = 2, 
  nu1  = nothing,
  nnu1 = 2, 
  # Physical parameters
  f    = 1.0,
  N    = 10.0,
  m    = 40.0,
  # Optional uniform and steady background flow
  Ub   = 0.0,
  Vb   = 0.0
  )

  if Ly == nothing; Ly = Lx; end
  if ny == nothing; ny = nx; end
  if nu0 == nothing; nu0 = 1e-1/(dt*(0.65π*nx/Lx)^nnu0); end
  if nu1 == nothing; nu1 = 1e-1/(dt*(0.65π*nx/Lx)^nnu1); end

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = RigidLidBoussinesq.Params(mu, nu0, nnu0, nu1, nnu1, f, N, m)
  vs = RigidLidBoussinesq.Vars(g)
  eq = RigidLidBoussinesq.Equation(pr, g)
  ts = ETDRK4TimeStepper(dt, eq.LC)

  FourierFlows.Problem(g, vs, pr, eq, ts)
end


# Params
abstract type RigidLidParams <: AbstractParams end

struct Params <: RigidLidParams
  mu::Float64     # Bottom drag 
  nu0::Float64    # Mode-0 viscosity
  nnu0::Int       # Mode-0 hyperviscous order
  nu1::Float64    # Mode-1 viscosity
  nnu1::Int       # Mode-1 hyperviscous order
  f::Float64      # Planetary vorticity
  N::Float64      # Buoyancy frequency
  m::Float64      # Mode-one wavenumber
  Ub::Float64     # Steady background barotropic x-velocity
  Vb::Float64     # Steady background barotropic y-velocity
end

function Params(mu, nu0, nnu0::Int, nu1, nnu1::Int, f, N, m; Ub=0.0, Vb=0.0)
  Params(mu, nu0, nnu0, nu1, nnu1, f, N, m, Ub, Vb)
end


# Equations
struct Equation <: AbstractEquation
  LC::Array{Complex{Float64}, 3}
  calcNL!::Function            
end

function Equation(p::RigidLidParams, g::TwoDGrid)
  LC = zeros(g.nx, g.ny, 4)
  LC[:, :, 1] = -p.nu0 * g.KKrsq.^(0.5*p.nnu0) - p.mu
  LC[:, :, 2] = -p.nu1 * g.KKrsq.^(0.5*p.nnu1)
  LC[:, :, 3] = -p.nu1 * g.KKrsq.^(0.5*p.nnu1)
  Equation(LC, calcNL!)
end


# Vars
abstract type TwoModeVars <: AbstractVars end

varlist = [
  :Z, :U, :V, :UZuz, :VZvz, :Ux, :Uy, :Vx, :Vy, :Psi,
  :u, :v, :p, :w, :zeta, :Uu, :Uv, :Up, :Vu, :Vv, :Vp, :uUxvUy, :uVxvVy
]

struct Vars <: TwoModeVars
  t::Float64
  sol::Array{Complex{Float64}, 3}

  # Auxiliary zeroth-mode vars
  Z::Array{Float64, 2}
  U::Array{Float64, 2}
  V::Array{Float64, 2}
  UZuz::Array{Float64, 2}
  VZvz::Array{Float64, 2}
  Ux::Array{Float64, 2}
  Uy::Array{Float64, 2}
  Vx::Array{Float64, 2}
  Vy::Array{Float64, 2}
  Psi::Array{Float64, 2}

  # Auxiliary first-mode vars
  u::Array{Float64, 2}
  v::Array{Float64, 2}
  w::Array{Float64, 2}
  p::Array{Float64, 2}
  zeta::Array{Float64, 2}

  # Multiplies
  Uu::Array{Float64, 2}
  Uv::Array{Float64, 2}
  Up::Array{Float64, 2}
  Vu::Array{Float64, 2}
  Vv::Array{Float64, 2}
  Vp::Array{Float64, 2}
  uUxvUy::Array{Float64, 2}
  uVxvVy::Array{Float64, 2}

  # Zeroth-mode transforms
  Zh::Array{Complex{Float64}, 2}
  Uh::Array{Complex{Float64}, 2}
  Vh::Array{Complex{Float64}, 2}
  UZuzh::Array{Complex{Float64}, 2}
  VZvzh::Array{Complex{Float64}, 2}
  Uxh::Array{Complex{Float64}, 2}
  Uyh::Array{Complex{Float64}, 2}
  Vxh::Array{Complex{Float64}, 2}
  Vyh::Array{Complex{Float64}, 2}
  Psih::Array{Complex{Float64}, 2}

  # First-mode transforms
  uh::Array{Complex{Float64}, 2}
  vh::Array{Complex{Float64}, 2}
  wh::Array{Complex{Float64}, 2}
  ph::Array{Complex{Float64}, 2}
  zetah::Array{Complex{Float64}, 2}

  # Multiply transforms
  Uuh::Array{Complex{Float64}, 2}
  Uvh::Array{Complex{Float64}, 2}
  Uph::Array{Complex{Float64}, 2}
  Vuh::Array{Complex{Float64}, 2}
  Vvh::Array{Complex{Float64}, 2}
  Vph::Array{Complex{Float64}, 2}
  uUxvUyh::Array{Complex{Float64}, 2}
  uVxvVyh::Array{Complex{Float64}, 2}
end

function Vars(g::TwoDGrid)
  # Initialize with t=0
  t = 0.0
  sol = zeros(Complex{Float64}, g.nkr, g.nl, 4)

  for var in varlist
    varh = Symbol(var, :h)
    eval(:($var = zeros(Float64, g.nx, g.ny)))
    eval(:($varh = zeros(Complex{Float64}, g.nkr, g.nl)))
  end

  Vars(t, sol,
    Z, U, V, UZuz, VZvz, Ux, Uy, Vx, Vy, Psi, 
    u, v, w, p, zeta, Uu, Uv, Up, Vu, Vv, Vp, uUxvUy, uVxvVy,
    Zh, Uh, Vh, UZuzh, VZvzh, Uxh, Uyh, Vxh, Vyh, Psih, 
    uh, vh, wh, ph, zetah, Uuh, Uvh, Uph, Vuh, Vvh, Vph, uUxvUyh, uVxvVyh)
end




# Solvers --------------------------------------------------------------------- 
function calcNL!(
  NL::Array{Complex{Float64}, 3}, sol::Array{Complex{Float64}, 3},
  t::Float64, v::Vars, p::RigidLidParams, g::TwoDGrid)

  @views v.Zh .= sol[:, :, 1]
  @views v.uh .= sol[:, :, 2]
  @views v.vh .= sol[:, :, 3]
  @views v.ph .= sol[:, :, 4]

  # Spectral-space calculations
  @. v.Psih = -g.invKKrsq*v.Zh

  @. v.Uh = -im*g.l*v.Psih
  @. v.Vh =  im*g.kr*v.Psih

  @. v.Uxh = im*g.kr*v.Uh
  @. v.Vxh = im*g.kr*v.Vh

  @. v.Uyh = im*g.l*v.Uh
  @. v.Vyh = im*g.l*v.Vh

  @. v.wh = -im/p.m*(g.kr*v.uh + g.l*v.vh)
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


  # Zeroth-mode nonlinear term
  @views @. NL[:, :, 1] = - im*g.kr*v.UZuzh - im*g.l*v.VZvzh

  # First-mode nonlinear terms:
  # u
  @views @. NL[:, :, 2] = ( p.f*sol[:, :, 3] - im*g.k*sol[:, :, 4]
    - im*g.kr*v.Uuh - im*g.l*v.Vuh - v.uUxvUyh
  )

  # v
  @views @. NL[:, :, 3] = ( -p.f*sol[:, :, 2] - im*g.l*solc[:, :, 4]
    - im*g.kr*v.Uvh - im*g.l*v.Vvh - v.uVxvVyh
  )

  # p
  @views @. NL[:, :, 4] = ( p.N^2.0/p.m*v.wh
    - im*g.kr*v.Uph - im*g.l*v.Vph
  )

  #dealias!(NL, g)

  nothing
end




# Helper functions ------------------------------------------------------------ 
function updatevars!(v::TwoModeVars, p::RigidLidParams, g::TwoDGrid, 
  Zh::AbstractArray)

  @views v.Zh .= v.sol[:, :, 1]
  @views v.uh .= v.sol[:, :, 2]
  @views v.vh .= v.sol[:, :, 3]
  @views v.ph .= v.sol[:, :, 4]

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

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.vars, prob.params, prob.grid)
end




"""
Set zeroth mode vorticity and update vars. 
"""
function set_Z!(v::Vars, p::RigidLidParams, g::TwoDGrid, Z)
  @views A_mul_B!(v.sol[:, :, 1], g.rfftplan, Z)
  updatevars!(v, p, g)
  nothing
end

function set_Z!(prob::AbstractProblem, Z)
  set_Z!(prob.vars, prob.params, prob.grid, Z)
end




""" 
Set first mode u, v, and p and update vars.
"""
function set_uvp!(vs::TwoModeVars, pr::RigidLidParams, g::TwoDGrid, u, v, p)
  @views A_mul_B!(vs.sol[:, :, 2], g.rfftplan, u)
  @views A_mul_B!(vs.sol[:, :, 3], g.rfftplan, v)
  @views A_mul_B!(vs.sol[:, :, 4], g.rfftplan, p)
  updatevars!(vs, pr, g)
  nothing
end

function set_uvp!(prob::AbstractProblem, u, v, p)
  set_uvp!(prob.vars, prob.params, prob.grid, u, v, p)
end




""" 
Set a plane wave solution with initial speed uw and non-dimensional wave
number nkw. The dimensional wavenumber will be 2π*nkw/Lx. 
"""
function set_planewave!(vs::TwoModeVars, pr::RigidLidParams, g::TwoDGrid,
  uw::Real, nkw::Int)

  x, y = g.X, g.Y

  # Wave parameters
  kw = 2π*nkw/g.Lx
  σ = sqrt(pr.f^2 + pr.N^2*kw^2/pr.m^2)
  alpha = pr.N^2*kw^2/(pr.f^2*pr.m^2) # also (sig^2-f^2)/f^2

  u0 = uw/2
  v0 = -uw * im*pr.f/2σ
  p0 = uw * kw*pr.N^2/(2σ*pr.m^2)

  u = u0 * exp.(im*kw*x)
  v = v0 * exp.(im*kw*x)
  p = p0 * exp.(im*kw*x)
  
  set_uvp!(vs, pr, g, u, v, p)
  nothing
end

function set_planewave!(prob::AbstractProblem, uw::Real, nkw::Int)
  set_planewave!(prob.vars, prob.params, prob.grid, uw::Real, nkw::Int)
end




""" 
Generate an isotropic spectrum of waves.
"""
function set_isotropicwavefield!(
  vs::TwoModeVars, pr::RigidLidParams, g::TwoDGrid, amplitude::Function; KE=1.0,
  maxspeed=nothing)

  # For clarity
  f, N, m = pr.f, pr.N, pr.m

  # Initialize
  phase = zeros(Float64, g.nx, g.ny)
  u0 = zeros(Float64, g.nx, g.ny)
  v0 = zeros(Float64, g.nx, g.ny)
  p0 = zeros(Float64, g.nx, g.ny)

  # Sum Fourier components
  for k in g.k, l in g.l
    if amplitude(k, l) > 1e-15
     
      # Dispersion relation
      σ = sqrt(f^2 + N^2/m^2*(k^2 + l^2))

      # Random phases
      phase .= k*g.X .+ l*g.Y .+ 2π*rand()

      # Polarization
      u0 .+= amplitude(k, l)*cos.(phase)
      v0 .+= -u0*(im*f/σ - k*l*N^2/(σ*m)^2)/(1 - (l*N)^2/(σ*m)^2)
      p0 .+= N^2/(σ*m^2) * (k*u0 .+ l*v0)
    end
  end

  
  if maxspeed == nothing # Normalize by kinetic energy
    uh, vh = fft(u0), fft(v0)
    ke = mode1ke(uh, vh, g)
    norm = sqrt(KE)/sqrt(ke/(g.Lx*g.Ly))
  else
    norm = maxspeed / maximum(sqrt.(u0.^2 + v0.^2))
  end

  u0 .*= norm
  v0 .*= norm
  p0 .*= norm
   
  set_uvp!(vs, pr, g, u0, v0, p0)

  nothing
end

function set_isotropicwavefield!(
  vs::TwoModeVars, pr::RigidLidParams, g::TwoDGrid; kwargs...)
  amplitude(k, l) = 1.0
  set_isotropicwavefield!(vs, pr, g, amplitude; kwargs...)
end

function set_isotropicwavefield!(prob::AbstractProblem, amplitude::Function;
  kwargs...)
  set_isotropicwavefield!(prob.vars, prob.params, prob.grid, amplitude; 
    kwargs...)
end


 

""" 
Returns the integrated energy in the zeroth mode energy. 
"""
function mode0energy(v::Vars, p::RigidLidParams, g::TwoDGrid)
  @views 0.5*FourierFlows.parsevalsum(g.invKKrsq.*abs2.(v.sol[:, :, 1]), g)
end

function mode0energy(prob::AbstractProblem)
  mode0energy(prob.vars, prob.params, prob.grid)
end




""" 
Returns the projection of the integrated first mode kinetic energy
onto the zeroth mode.
"""
function mode1ke(uh, vh, g)
  FourierFlows.parsevalsum2(uh, g) + FourierFlows.parsevalsum2(vh, g)
end

function mode1ke(v::TwoModeVars, p::RigidLidParams, g::TwoDGrid)
  @views mode1ke(v.sol[:, :, 2], v.sol[:, :, 3], g)
end

function mode1ke(prob::AbstractProblem)
  mode1ke(prob.vars, prob.params, prob.grid)
end




""" 
Returns the projection of the integrated first mode potential energy onto the
zeroth mode. 
"""
function mode1pe(v::TwoModeVars, p::RigidLidParams, g::TwoDGrid)
  p.m^2/p.N^2*FourierFlows.parsevalsum2(v.sol[:, :, 4], g)
end

function mode1pe(prob::AbstractProblem)
  mode1pe(prob.vars, prob.params, prob.grid)
end




""" 
Returns the projection of the total integrated first mode energy onto the
zeroth mode.
"""
function mode1energy(v::TwoModeVars, p::RigidLidParams, g::TwoDGrid)
  mode1ke(v, p, g) + mode1pe(v, p, g)
end

function mode1energy(prob::AbstractProblem)
  mode1energy(prob.vars, prob.params, prob.grid)
end




""" 
Returns the total energy projected onto the zeroth mode.
"""
function totalenergy(v::TwoModeVars, p::RigidLidParams, g::TwoDGrid)
  mode0energy(v, p, g) + mode1energy(v, p, g)
end

function totalenergy(prob::AbstractProblem)
  totalenergy(prob.vars, prob.params, prob.grid)
end




""" 
Returns kinetic energy dissipation of the zeroth mode. 
"""
function mode0dissipation(v::TwoModeVars, p::RigidLidParams, g::TwoDGrid)
  delzeta = irfft(
    (-1.0)^(p.nnu0/2) .* g.KKrsq.^(p.nnu0/2) .* vs.sol[:, :, 1], g.nx)
  -p.nu*g.dx*g.dy*sum(vs.Psi.*delzeta)
end


"""
Return the speed associated with mode-0.
"""
mode1speed(prob) = sqrt.(prob.vars.u.^2.0 .+ prob.vars.v.^2.0)


"""
Return the speed associated with mode-0.
"""
mode0speed(prob) = sqrt.(prob.vars.U.^2.0 .+ prob.vars.V.^2.0)


""" 
Returns the Courant-Freidrichs-Lewy number. 
"""
function CFL(prob, dt)
  dx = minimum([prob.grid.dx, prob.grid.dy])
  U = maximum([
    maximum(abs.(prob.vars.U)), 
    maximum(abs.(prob.vars.V)),
    maximum(abs.(prob.vars.u)),
    maximum(abs.(prob.vars.v))])

  U*dt/dx
end

 


# End module
end
