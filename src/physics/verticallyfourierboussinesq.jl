module VerticallyFourierBoussinesq
using FourierFlows

# Problem
"""
    InitialValueProblem(; parameters...)

Construct a VerticallyFourierBoussinesq initial value problem.
"""
function InitialValueProblem(;
   nx = 128,
   Lx = 2π,
   ny = nothing,
   Ly = nothing,
   nu0 = 0,
  nnu0 = 1,
   nu1 = 0,
  nnu1 = 1,
    f = 1,
    N = 10,
    m = 40,
   Ub = 0,
   Vb = 0,
   dt = 0.01,
  stepper = "RK4"
  )

  if Ly == nothing; Ly = Lx; end
  if ny == nothing; ny = nx; end

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = Params(nu0, nnu0, nu1, nnu1, f, N, m)
  vs = Vars(g)
  eq = Equation(pr, g)
  st = DualState(Complex{Float64}, (g.nk, g.nl, 3), 
                 Complex{Float64}, (g.nkr, g.nl), dt)
  ts = FourierFlows.autoconstructtimestepper(stepper, dt, eq.LCc, eq.LCr, g)

  FourierFlows.Problem(g, vs, pr, eq, ts, st)
end


# Params
abstract type TwoModeParams <: AbstractParams end

"""
    Params(nu0, nnu0, nu1, nnu1, f, N, m)
    Params(nu0, nnu0, nu1, nnu1, f, N, m, Ub, Vb) 

Construct parameters for the Two-Fourier-mode Boussinesq problem. Suffix 0
refers to zeroth mode; 1 to first mode. f, N, m are Coriolis frequency, 
buoyancy frequency, and vertical wavenumber of the first mode, respectively.
The optional constant background velocity (Ub,Vb) is set to zero by default.
The viscosity is applied only to the first-mode horizontal velocities.
"""
struct Params <: TwoModeParams
  nu0::Float64                 # Mode-0 viscosity
  nnu0::Int                    # Mode-0 hyperviscous order
  nu1::Float64                 # Mode-1 viscosity
  nnu1::Int                    # Mode-1 hyperviscous order
  f::Float64                  # Planetary vorticity
  N::Float64                  # Buoyancy frequency
  m::Float64                  # Mode-one wavenumber
  Ub::Float64                 # Steady mode-0 mean x-velocity
  Vb::Float64                 # Steady mode-0 mean y-velocity
end

Params(nu0, nnu0, nu1, nnu1, f, N, m, Ub=0) = Params(nu0, nnu0, nu1, nnu1, f, N, m, 
  Ub, 0)

# Equations
function Equation(p::TwoModeParams, g::TwoDGrid)
  LCc, LCr = getlinearcoefficients(p, g)
  FourierFlows.DualEquation{3,2}(LCc, LCr, calcN!)
end

function getlinearcoefficients(p::TwoModeParams, g::TwoDGrid)
  LCr = -p.nu0 * g.KKrsq.^p.nnu0
  LCc = zeros(g.nk, g.nl, 3)
  LCc[:, :, 1] = -p.nu1 * g.KKsq.^p.nnu1
  LCc[:, :, 2] = -p.nu1 * g.KKsq.^p.nnu1
  LCc, LCr
end

# Vars
abstract type TwoModeVars <: AbstractVars end

struct Vars <: TwoModeVars
  # Zeroth-mode
  Z::Array{Float64,2}
  U::Array{Float64,2}
  V::Array{Float64,2}
  UZuzvw::Array{Float64,2}
  VZvzuw::Array{Float64,2}
  Ux::Array{Float64,2}
  Uy::Array{Float64,2}
  Vx::Array{Float64,2}
  Vy::Array{Float64,2}
  Psi::Array{Float64,2}
  # First-mode
  u::Array{Complex{Float64},2}
  v::Array{Complex{Float64},2}
  w::Array{Complex{Float64},2}
  p::Array{Complex{Float64},2}
  zeta::Array{Complex{Float64},2}
  # Multiplies
  Uu::Array{Complex{Float64},2}
  Uv::Array{Complex{Float64},2}
  Up::Array{Complex{Float64},2}
  Vu::Array{Complex{Float64},2}
  Vv::Array{Complex{Float64},2}
  Vp::Array{Complex{Float64},2}
  uUxvUy::Array{Complex{Float64},2}
  uVxvVy::Array{Complex{Float64},2}
  # Zeroth-mode transforms
  Zh::Array{Complex{Float64},2}
  Uh::Array{Complex{Float64},2}
  Vh::Array{Complex{Float64},2}
  UZuzvwh::Array{Complex{Float64},2}
  VZvzuwh::Array{Complex{Float64},2}
  Uxh::Array{Complex{Float64},2}
  Uyh::Array{Complex{Float64},2}
  Vxh::Array{Complex{Float64},2}
  Vyh::Array{Complex{Float64},2}
  Psih::Array{Complex{Float64},2}
  # First-mode transforms
  uh::Array{Complex{Float64},2}
  vh::Array{Complex{Float64},2}
  wh::Array{Complex{Float64},2}
  ph::Array{Complex{Float64},2}
  zetah::Array{Complex{Float64},2}
  # Multiply transforms
  Uuh::Array{Complex{Float64},2}
  Uvh::Array{Complex{Float64},2}
  Uph::Array{Complex{Float64},2}
  Vuh::Array{Complex{Float64},2}
  Vvh::Array{Complex{Float64},2}
  Vph::Array{Complex{Float64},2}
  uUxvUyh::Array{Complex{Float64},2}
  uVxvVyh::Array{Complex{Float64},2}
end

function Vars(g)
  @createarrays Float64 (g.nx, g.ny) Z U V UZuzvw VZvzuw Ux Uy Vx Vy Psi
  @createarrays Complex{Float64} (g.nx, g.ny) u v w p zeta Uu Uv Up Vu Vv Vp
  @createarrays Complex{Float64} (g.nx, g.ny) uUxvUy uVxvVy
  @createarrays Complex{Float64} (g.nkr, g.nl) Zh Uh Vh UZuzvwh VZvzuwh
  @createarrays Complex{Float64} (g.nkr, g.nl) Uxh Uyh Vxh Vyh Psih
  @createarrays Complex{Float64} (g.nk, g.nl) uh vh wh ph zetah
  @createarrays Complex{Float64} (g.nk, g.nl) Uuh Uvh Uph Vuh Vvh Vph
  @createarrays Complex{Float64} (g.nk, g.nl) uUxvUyh uVxvVyh
  return Vars(
    Z, U, V, UZuzvw, VZvzuw, Ux, Uy, Vx, Vy, Psi,
    u, v, w, p, zeta, Uu, Uv, Up, Vu, Vv, Vp, uUxvUy, uVxvVy,
    Zh, Uh, Vh, UZuzvwh, VZvzuwh, Uxh, Uyh, Vxh, Vyh, Psih,
    uh, vh, wh, ph, zetah, Uuh, Uvh, Uph, Vuh, Vvh, Vph, uUxvUyh, uVxvVyh,
    )
end


# Solvers
function calcN!(
  Nc::Array{Complex{Float64},3}, Nr::Array{Complex{Float64},2},
  solc::Array{Complex{Float64},3}, solr::Array{Complex{Float64},2},
  t::Float64, s::DualState, v::Vars, p::TwoModeParams, g::TwoDGrid)

  v.Zh .= solr

  @. v.Psih = -g.invKKrsq*v.Zh

  @. v.Uh  = -im*g.l  * v.Psih
  @. v.Vh  =  im*g.kr * v.Psih
  @. v.Uxh =  im*g.kr * v.Uh
  @. v.Uyh =  im*g.l  * v.Uh
  @. v.Vxh =  im*g.kr * v.Vh
  @. v.Vyh =  im*g.l  * v.Vh

  v.Uh[1, 1] += p.Ub*g.nx*g.ny
  v.Vh[1, 1] += p.Vb*g.nx*g.ny

  # Inverse transforms
  A_mul_B!(v.Z,  g.irfftplan, v.Zh)
  A_mul_B!(v.U,  g.irfftplan, v.Uh)
  A_mul_B!(v.V,  g.irfftplan, v.Vh)
  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)

  @views A_mul_B!(v.u, g.ifftplan, solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, solc[:, :, 3])

  @views @. v.zetah = im*g.k*solc[:, :, 2] - im*g.l*solc[:, :, 1]
  @views @. v.wh = -(g.k*solc[:, :, 1] + g.l*solc[:, :, 2]) / p.m

  A_mul_B!(v.w, g.ifftplan, v.wh)
  A_mul_B!(v.zeta, g.ifftplan, v.zetah)

  # Multiplies
  @. v.UZuzvw = (v.U * v.Z
    + real(   v.u*conj(v.zeta)  +  im*p.m*v.v*conj(v.w)
            + conj(v.u)*v.zeta  -  im*p.m*conj(v.v)*v.w   ))

  @. v.VZvzuw = (v.V * v.Z
    + real(   v.v*conj(v.zeta)  -  im*p.m*v.u*conj(v.w)
            + conj(v.v)*v.zeta  +  im*p.m*conj(v.u)*v.w   ))

  @. v.uUxvUy = v.u*v.Ux + v.v*v.Uy
  @. v.uVxvVy = v.u*v.Vx + v.v*v.Vy

  @. v.Uu = v.U * v.u
  @. v.Vu = v.V * v.u
  @. v.Uv = v.U * v.v
  @. v.Vv = v.V * v.v
  @. v.Up = v.U * v.p
  @. v.Vp = v.V * v.p

  # Forward transforms
  A_mul_B!(v.UZuzvwh, g.rfftplan, v.UZuzvw)
  A_mul_B!(v.VZvzuwh, g.rfftplan, v.VZvzuw)
  A_mul_B!(v.uUxvUyh, g.fftplan, v.uUxvUy)
  A_mul_B!(v.uVxvVyh, g.fftplan, v.uVxvVy)

  A_mul_B!(v.Uuh, g.fftplan, v.Uu)
  A_mul_B!(v.Uvh, g.fftplan, v.Uv)
  A_mul_B!(v.Vuh, g.fftplan, v.Vu)
  A_mul_B!(v.Vvh, g.fftplan, v.Vv)
  A_mul_B!(v.Uph, g.fftplan, v.Up)
  A_mul_B!(v.Vph, g.fftplan, v.Vp)

  # Zeroth-mode nonlinear term
  @. Nr = - im*g.kr*v.UZuzvwh - im*g.l*v.VZvzuwh

  # First-mode nonlinear terms:
  # u
  @views @. Nc[:, :, 1] = (  p.f*solc[:, :, 2] - im*g.k*(solc[:, :, 3] + v.Uuh)
    - im*g.l*v.Vuh - v.uUxvUyh )

  # v
  @views @. Nc[:, :, 2] = ( -p.f*solc[:, :, 1] - im*g.l*(solc[:, :, 3] + v.Vvh)
    - im*g.k*v.Uvh - v.uVxvVyh )

  # p
  @views @. Nc[:, :, 3] = ( im*p.N^2.0/p.m*v.wh
    - im*g.k*v.Uph - im*g.l*v.Vph )

  nothing
end


# Helper functions
"""
    updatevars!(prob)
    updatevars!(v, s, p, g)

Update variables to correspond to the solution in s.sol or prob.state.sol.
"""
function updatevars!(v, s, p, g)
  v.Zh .= s.solr
  @. v.Psih = -g.invKKrsq*v.Zh
  @. v.Uh   = -im*g.l*v.Psih
  @. v.Vh   =  im*g.kr*v.Psih

  Psih1 = deepcopy(v.Psih)
  Zh1 = deepcopy(v.Zh)
  Uh1 = deepcopy(v.Uh)
  Vh1 = deepcopy(v.Vh)

  A_mul_B!(v.Psi, g.irfftplan, Psih1)
  A_mul_B!(v.Z, g.irfftplan, Zh1)
  A_mul_B!(v.U, g.irfftplan, Uh1)
  A_mul_B!(v.V, g.irfftplan, Vh1)

  @views v.uh .= s.solc[:, :, 1]
  @views v.vh .= s.solc[:, :, 2]
  @views v.ph .= s.solc[:, :, 3]

  @. v.wh = -1.0/p.m*(g.k*v.uh + g.l*v.vh)

  A_mul_B!(v.u, g.ifftplan, v.uh)
  A_mul_B!(v.v, g.ifftplan, v.vh)
  A_mul_B!(v.p, g.ifftplan, v.ph)
  A_mul_B!(v.w, g.ifftplan, v.wh)
  nothing
end

"""
    set_Z!(prob, Z)
    set_Z!(s, v, p, g, Z)

Set zeroth mode vorticity and update variables.
"""
function set_Z!(s, v, p, g, Z)
  A_mul_B!(s.solr, g.rfftplan, Z)
  updatevars!(v, s, p, g)
  nothing
end

"""
    set_uvp!(prob, u, v, p)
    set_uvp!(st, vs, pr, g, u, v, p)

Set first mode u, v, and p and update variables.
"""
function set_uvp!(s, vs, pr, g, u, v, p)
  uh = fft(u)
  vh = fft(v)
  ph = fft(p)

  @. s.solc[:, :, 1] = uh
  @. s.solc[:, :, 2] = vh
  @. s.solc[:, :, 3] = ph

  updatevars!(vs, s, pr, g)
  nothing
end

set_uvp!(prob, u, v, p) = set_uvp!(prob.state, prob.vars, prob.params, 
  prob.grid, u, v, p)
updatevars!(prob) = updatevars!(prob.vars, prob.state, prob.params, prob.grid)
set_Z!(prob, Z) = set_Z!(prob.state, prob.vars, prob.params, prob.grid, Z)


"""
    set_planewave!(prob, uw, nkw)
    set_planewave!(s, vs, pr, g, uw, nkw)

Set a plane wave solution with initial speed uw and non-dimensional wave
number nkw. The dimensional wavenumber will be 2π*nkw/Lx.
"""
function set_planewave!(s, vs, pr, g, uw, nkw)
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

  set_uvp!(s, vs, pr, g, u, v, p)
  nothing
end


"""
    set_isotropicwavefield!(prob, amplitude; KE=1.0, maxspeed=nothing)

Generate an isotropic spectrum of waves with an envelope given by the function
amplitude(k, l), and either total kinetic energy KE or maximum speed maxspeed.
"""
function set_isotropicwavefield!(s, vs, pr, g, amplitude; 
                                 KE=1.0, maxspeed=nothing)
  f, N, m = pr.f, pr.N, pr.m # for clarity
  @createarrays Complex{Float64} (g.nx, g.ny) phase u0 v0 p0
  
  # Sum Fourier components
  for k in real.(g.k), l in real.(g.l)
    if amplitude(k, l) > 1e-15
      σ = sqrt(f^2 + N^2/m^2*(k^2 + l^2))   # dispersion relation
      phase .= k*g.X .+ l*g.Y .+ 2π*rand()  # random phases
      # Polarization relations
      u0 .+= amplitude(k, l)*exp.(im*phase)
      v0 .+= -u0*(im*f/σ - k*l*N^2/(σ*m)^2)/(1 - (l*N)^2/(σ*m)^2)
      p0 .+= N^2/(σ*m^2) * (k*u0 .+ l*v0)
    end
  end

  if maxspeed == nothing # Normalize by kinetic energy
    uh, vh = fft(u0), fft(v0)
    ke = mode1ke(uh, vh, g)
    norm = sqrt(KE)/sqrt(ke/(g.Lx*g.Ly))
  else
    norm = maxspeed / maximum(
      sqrt.(real.(u0+conj.(u0)).^2 + real.(v0+conj.(v0)).^2))
  end

  u0 .*= norm
  v0 .*= norm
  p0 .*= norm

  set_uvp!(s, vs, pr, g, u0, v0, p0)
  nothing
end

set_planewave!(prob, uw, nkw) = set_planewave!(prob.state, prob.vars, 
  prob.params, prob.grid, uw::Real, nkw::Int)

function set_isotropicwavefield!(prob::AbstractProblem, amplitude; kwargs...)
  set_isotropicwavefield!(prob.state, prob.vars, prob.params, prob.grid, 
    amplitude; kwargs...)
end


"""
    mode0energy(prob)

Returns the integrated energy in the zeroth mode.
"""
function mode0energy(s::DualState, g::TwoDGrid)
  0.5*FourierFlows.parsevalsum(real.(g.invKKrsq).*abs2.(s.solr), g)
end

function mode0energy(prob::AbstractProblem)
  mode0energy(prob.state, prob.grid)
end


"""
    mode1ke(prob)

Returns the projection of the integrated first mode kinetic energy
onto the zeroth mode.
"""
function mode1ke(uh, vh, g)
  FourierFlows.parsevalsum2(uh, g) + FourierFlows.parsevalsum2(vh, g)
end

mode1ke(s, g) = @views mode1ke(s.solc[:, :, 1], s.solc[:, :, 2], g)

"""
    mode1pe(prob)

Returns the projection of the integrated first mode potential energy onto the
zeroth mode.
"""
function mode1pe(s, p, g)
  p.m^2/p.N^2*FourierFlows.parsevalsum2(s.solc[:, :, 3], g)
end

"""
    mode1energy(s, p, g)
Returns the projection of the total integrated first mode energy onto the
zeroth mode.
"""
mode1energy(s, p, g) = mode1ke(s, g) + mode1pe(s, p, g)

    mode1pe(prob) = mode1pe(prob.state, prob.params, prob.grid)
    mode1ke(prob) = mode1ke(prob.state, prob.grid)
mode1energy(prob) = mode1energy(prob.state, prob.params, prob.grid)

"""
    totalenergy(s, p, g)

Returns the total energy projected onto the zeroth mode.
"""
totalenergy(s, p, g) = mode0energy(s, g) + mode1energy(s, p, g)
totalenergy(prob) = totalenergy(prob.state, prob.params, prob.grid)

"""
    mode0dissipation(s, v, p, g)

Returns kinetic energy dissipation of the zeroth mode.
"""
function mode0dissipation(s, v, p, g)
  v.Zh .= s.solr
  @. v.Psih = -g.invKKrsq*v.Zh
  @. v.UZuzvwh = -1.0^(p.nnu0/2) * g.KKrsq.^(p.nnu0/2) * s.solr
  A_mul_B!(v.UZuzvw, g.irfftplan, v.UZuzvwh) 
  A_mul_B!(v.Psi, g.irfftplan, v.Psih) 
  @. v.VZvzuw = v.Psi*v.UZuzvw

  -p.nu*g.dx*g.dy*sum(v.VZvzuw)
end

"""
    shearproduction(s, v, p, g)

Returns the domain-integrated shear production.
"""
function shearproduction(s, v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  v.Zh .= s.solr
  @. v.Psih = -g.invKKrsq*v.Zh
  @. v.Uh  = -im*g.l  * v.Psih
  @. v.Vh  =  im*g.kr * v.Psih
  @. v.Uxh =  im*g.kr * v.Uh
  @. v.Uyh =  im*g.l  * v.Uh
  @. v.Vxh =  im*g.kr * v.Vh
  @. v.Vyh =  im*g.l  * v.Vh

  A_mul_B!(v.Ux, g.irfftplan, v.Uxh)
  A_mul_B!(v.Uy, g.irfftplan, v.Uyh)
  A_mul_B!(v.Vx, g.irfftplan, v.Vxh)
  A_mul_B!(v.Vy, g.irfftplan, v.Vyh)

  @views A_mul_B!(v.u, g.ifftplan, s.solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, s.solc[:, :, 2])

  @. v.UZuzvw = real( 2.0*abs2(v.u)*v.Ux + 2.0*abs2(v.v)*v.Vy
                       + (conj(v.u)*v.v + v.u*conj(v.v))*(v.Uy + v.Vx) )

  g.dx*g.dy*sum(v.UZuzvw)
end

shearproduction(prob) = shearproduction(prob.state, prob.vars, prob.params, 
                                        prob.grid)

"""
Return the domain-integrated conversion from potential to kinetic energy.
"""
function energyconversion(s, v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  @views @. v.wh = -(g.k*s.solc[:, :, 1] + g.l*s.solc[:, :, 2]) / p.m
  @views A_mul_B!(v.p, g.ifftplan, s.solc[:, :, 3])
  A_mul_B!(v.w, g.ifftplan, v.wh)
  # b = i*m*p
  @. v.UZuzvw = real(im*p.m*conj(v.w)*v.p - im*p.m*v.w*conj(v.p))
  g.dx*g.dy*sum(v.UZuzvw)
end

function energyconversion(prob::AbstractProblem)
  energyconversion(prob.state, prob.vars, prob.params, prob.grid)
end


"""
Returns the projection of available potential vorticity onto the
zeroth mode.
"""
function mode0apv(Z, u, v, p, pr::TwoModeParams, g::TwoDGrid)
  (Z .+ irfft( pr.m^2.0./pr.N^2.0 .* (
      im.*g.l.*rfft(real.(u.*conj.(p) .+ conj.(u).*p))
    - im.*g.kr.*rfft(real.(v.*conj.(p) .+ conj.(v).*p))
  ), g.nx))
end

function mode0apv(s, v::Vars, p::TwoModeParams, g::TwoDGrid)
  v.Z = irfft(s.solr, g.nx)
  @views A_mul_B!(v.u, g.ifftplan, s.solc[:, :, 1])
  @views A_mul_B!(v.v, g.ifftplan, s.solc[:, :, 2])
  @views A_mul_B!(v.p, g.ifftplan, s.solc[:, :, 3])
  mode0apv(v.Z, v.u, v.v, v.p, p, g)
end

function mode0apv(prob::AbstractProblem)
  mode0apv(prob.state, prob.vars, prob.params, prob.grid)
end


"""
Returns the projection of available potential energy onto the first mode.
"""
function mode1apv(Z, zeta, p, pr::TwoModeParams, g::TwoDGrid)
  zeta .- pr.m.^2.0./pr.N.^2.0 .* (pr.f .+ Z) .* p
end

function mode1apv(Z, s, v::TwoModeVars, p::TwoModeParams, g::TwoDGrid)
  @views @. v.ph = s.solc[:, :, 3]
  @views @. v.zetah = im*g.k*s.solc[:, :, 2] - im*g.l*s.solc[:, :, 1]

  A_mul_B!(v.p,  g.ifftplan, v.ph)
  A_mul_B!(v.zeta,  g.ifftplan, v.zetah)

  mode1apv(Z, v.zeta, v.p, p, g)
end


"""
Return the apv associated with mode-1.
"""

function mode1apv(s, v::Vars, p::TwoModeParams, g::TwoDGrid)
  v.Z = irfft(s.solr, g.nx)
  mode1apv(v.Z, v, p, g)
end

mode1apv(prob::AbstractProblem) = mode1apv(prob.state, prob.vars, prob.params, 
                                           prob.grid)

"""
    mode1u(prob)

Return the x-velocity associated with mode-1 at z=0.
"""
mode1u(v::AbstractVars) = real.(v.u .+ conj.(v.u))
mode1u(prob::AbstractProblem) = mode1u(prob.vars)

"""
    mode1v(prob)

Return the y-velocity associated with mode-1 at z=0.
"""
mode1v(v::AbstractVars) = real.(v.v .+ conj.(v.v))
mode1v(prob::AbstractProblem) = mode1v(prob.vars)

"""
    mode1w(prob)

Return the z-velocity associated with mode-1 at z=0.
"""
mode1w(v::AbstractVars) = real.(v.w .+ conj.(v.w))
mode1w(prob::AbstractProblem) = mode1w(prob.vars)

"""
    mode1p(prob)

Return the pressure associated with mode-1 at z=0.
"""
mode1p(v::AbstractVars) = real.(v.p .+ conj.(v.p))
mode1p(prob::AbstractProblem) = mode1p(prob.vars)

"""
    mode1buoyancy(prob)
Return the buoyancy associated with mode-1 at z=0.
"""
mode1buoyancy(v, p) = real.(im.*p.m.*v.p .- im.*p.m.*conj.(v.p))
mode1buoyancy(prob::AbstractProblem) = mode1buoyancy(prob.vars, prob.params)

"""
    mode1speed(prob)

Return the speed associated with mode-1 at z=0.
"""
mode1speed(v::AbstractVars) = sqrt.(mode1u(v).^2.0 .+ mode1v(v).^2.0)
mode1speed(prob::AbstractProblem) = mode1speed(prob.vars)

"""
    mode0speed(prob)

Return the speed associated with mode-0.
"""
mode0speed(prob) = sqrt.(prob.vars.U.^2.0 .+ prob.vars.V.^2.0)

# End module
end
