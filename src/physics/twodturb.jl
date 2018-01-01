__precompile__()

module TwoDTurb

using FourierFlows,
      PyPlot

export Params,
       Vars,
       Equation

export set_q!, updatevars!


# Problem
function InitialValueProblem(;
  nx = 256,
  Lx = 2π,
  ny = nothing,
  Ly = nothing,
  nu = nothing,
  nnu = 2,
  dt = 0.01,
  withfilter = false
  )

  if Ly == nothing; Ly = Lx; end
  if ny == nothing; ny = nx; end
  if nu == nothing 
    if !withfilter; nu = 0.0
    else            nu = 1e-1/(dt*(0.65π*nx/Lx)^nnu)
    end
  end

  g  = TwoDGrid(nx, Lx, ny, Ly)
  pr = TwoDTurb.Params(nu, nnu)
  vs = TwoDTurb.Vars(g)
  eq = TwoDTurb.Equation(pr, g)
  if withfilter
    ts = FilteredETDRK4TimeStepper(dt, eq.LC, g)
  else
    ts = ETDRK4TimeStepper(dt, eq.LC)
  end

  FourierFlows.Problem(g, vs, pr, eq, ts)
end

function InitialValueProblem(n, L, nu, nnu, dt, withfilter)
  InitialValueProblem(nx=n, Lx=L, nu=nu, nnu=nnu, dt=dt, withfilter=withfilter)
end


# P A R A M S
struct Params <: AbstractParams
  nu::Float64  # Vorticity viscosity
  nnu::Int     # Vorticity hyperviscous order
end


# E Q U A T I O N S
struct Equation <: AbstractEquation
  LC::Array{Complex{Float64},2}  # Element-wise coeff of the eqn's linear part
  calcNL!::Function               # Function to calculate eqn's nonlinear part
end

function Equation(p::Params, g::TwoDGrid)
  LC = -p.nu*g.KKrsq.^(0.5*p.nnu)
  Equation(LC, calcNL!)
end




# V A R S
mutable struct Vars <: AbstractVars
  t::Float64
  sol::Array{Complex{Float64},2}

  q::Array{Float64,2}
  U::Array{Float64,2}
  V::Array{Float64,2}
  Uq::Array{Float64,2}
  Vq::Array{Float64,2}
  psi::Array{Float64,2}

  qh::Array{Complex{Float64},2}
  Uh::Array{Complex{Float64},2}
  Vh::Array{Complex{Float64},2}
  Uqh::Array{Complex{Float64},2}
  Vqh::Array{Complex{Float64},2}
  psih::Array{Complex{Float64},2}
end

function Vars(g::TwoDGrid)
  # Initialize with t=0
  t = 0.0
  sol  = zeros(Complex{Float64}, g.nkr, g.nl)

  # Vorticity auxiliary vars
  q    = zeros(Float64, g.nx, g.ny)
  U    = zeros(Float64, g.nx, g.ny)
  V    = zeros(Float64, g.nx, g.ny)
  Uq   = zeros(Float64, g.nx, g.ny)
  Vq   = zeros(Float64, g.nx, g.ny)
  psi  = zeros(Float64, g.nx, g.ny)

  qh   = zeros(Complex{Float64}, g.nkr, g.nl)
  Uh   = zeros(Complex{Float64}, g.nkr, g.nl)
  Vh   = zeros(Complex{Float64}, g.nkr, g.nl)
  Uqh  = zeros(Complex{Float64}, g.nkr, g.nl)
  Vqh  = zeros(Complex{Float64}, g.nkr, g.nl)
  psih = zeros(Complex{Float64}, g.nkr, g.nl)

  # Random initial condition
  sol = exp.( 2.0*pi*im*rand(g.nkr, g.nl) )

  Vars(t, sol, q, U, V, Uq, Vq, psi, qh, Uh, Vh, Uqh, Vqh, psih)
end




# S O L V E R S
function calcNL!(NL::Array{Complex{Float64},2}, sol::Array{Complex{Float64},2},
  t::Float64, v::Vars, p::Params, g::TwoDGrid)

  @. v.qh = sol
  @. v.Uh =  im*g.l *g.invKKrsq*sol
  @. v.Vh = -im*g.kr*g.invKKrsq*sol

  A_mul_B!(v.q, g.irfftplan, v.qh)
  A_mul_B!(v.U, g.irfftplan, v.Uh)
  A_mul_B!(v.V, g.irfftplan, v.Vh)

  @. v.Uq = v.U*v.q
  @. v.Vq = v.V*v.q

  A_mul_B!(v.Uqh, g.rfftplan, v.Uq)
  A_mul_B!(v.Vqh, g.rfftplan, v.Vq)

  @. NL = -im*g.kr*v.Uqh - im*g.l*v.Vqh

  nothing
end




# H E L P E R   F U N C T I O N S
function updatevars!(v::Vars, g::TwoDGrid)
  v.qh .= v.sol
  @. v.psih = -g.invKKrsq*v.qh
  @. v.Uh = -im*g.l  * v.psih
  @. v.Vh =  im*g.kr * v.psih

  qh = deepcopy(v.qh)
  Uh = deepcopy(v.Uh)
  Vh = deepcopy(v.Vh)

  A_mul_B!(v.q, g.irfftplan, qh)
  A_mul_B!(v.U, g.irfftplan, Uh)
  A_mul_B!(v.V, g.irfftplan, Vh)

  nothing
end

function updatevars!(v::Vars, p::Params, g::TwoDGrid)
  updatevars!(v, g)
end

function updatevars!(prob::AbstractProblem)
  updatevars!(prob.vars, prob.grid)
end



""" Set the vorticity field. """
function set_q!(v::Vars, g::TwoDGrid, q::Array{Float64, 2})
  A_mul_B!(v.sol, g.rfftplan, q)
  updatevars!(v, g)
end

function set_q!(v::Vars, p::Params, g::TwoDGrid, q::Array{Float64, 2})
  set_q!(v, g, q)
end

function set_q!(prob::AbstractProblem, q)
  set_q!(prob.vars, prob.grid, q)
end




""" Calculate the domain integrated kinetic energy. """
function energy(v::Vars, g::TwoDGrid)
  0.5*(FourierFlows.parsevalsum2(g.Kr.*g.invKKrsq.*v.sol, g)
        + FourierFlows.parsevalsum2(g.Lr.*g.invKKrsq.*v.sol, g))
end

function energy(prob::AbstractProblem)
  energy(prob.vars, prob.grid)
end




"""
Returns the domain-integrated enstrophy.
"""
function enstrophy(v::Vars, g::TwoDGrid)
  0.5*FourierFlows.parsevalsum2(v.sol, g)
end

function enstrophy(prob)
  enstrophy(prob.vars, prob.grid)
end




""" Make a field of mature turbulence on a square grid.

  Args:
    nx: grid resolution
    Lx: grid extent
    qf: final maximum vorticity
    q0: initial maximum vorticity
    nnu: order of hyperviscosity
    maxsteps: maximum number of steps to take
    dt: time step
    nu: hyperviscosity
    k0: initial wavenumber
    E0: initial energy
    tf: final time
    plots: whether or not to plot field evolution

  Returns
    q: The vorticity field
"""
function makematureturb(nx::Int, Lx::Real; qf=0.1, q0=0.2, nnu=4,
  maxsteps=10000, dt=nothing, nu=nothing, k0=nx/2,
  E0=nothing, tf=nothing, plots=false, loginterval=5)

  g  = TwoDGrid(nx, Lx)
  vs = TwoDTurb.Vars(g)

  if E0 != nothing # set initial energy rather than vorticity

    # Closely following the formulation in Rocha, Wagner, Young
    modk = sqrt(g.KKsq)

    psik = zeros(g.nk, g.nl)
    psik =  (modk .* (1 + (modk/k0).^4)).^(-0.5)
    psik[1, 1] = 0.0
    C = real(sqrt(E0/sum(g.KKsq.*abs2.(psik))))

    psi = zeros(g.nx, g.ny)
    for i = 1:128
      for j = 1:128
        psi .+= real.(C*psik[i, j]*cos.(
          g.k[i]*g.X + g.l[j]*g.Y + 2*pi*rand(1)[1]))
      end
    end

    psih = rfft(psi)
    qi = -irfft(g.KKrsq.*psih, g.nx)
    set_q!(vs, g, qi)
    E0 = FourierFlows.parsevalsum(g.KKrsq.*abs2.(psih), g)

  else
    qi = FourierFlows.peaked_isotropic_spectrum(nx, k0; maxval=q0)
    set_q!(vs, g, qi)
    E0 = energy(vs, g)
  end

  maxq = q0 = maximum(abs.(vs.q))

  # Defaults
  if dt == nothing; dt = 0.1*g.dx/maximum([vs.U; vs.V]);  end
  if nu == nothing; nu = 0.1/(dt*(0.65*nx/Lx)^nnu);          end
  if tf != nothing; maxsteps = ceil(Int, tf/dt); qf=0.0   end

  # Number of substeps between vorticity-checking
  substeps = ceil(Int, loginterval/(maxq*dt))

  pr = TwoDTurb.Params(nu, nnu)
  eq = TwoDTurb.Equation(pr, g)
  ts = ETDRK4TimeStepper(dt, eq.LC)

  if plots
    fig, axs = subplots()
    imshow(vs.q)
    pause(0.01)
  end

  @printf("\nMaking a mature turbulence field...\n")
  starttime = time()
  while maxq > qf && ts.step < maxsteps

    stepforward!(vs, ts, eq, pr, g; nsteps=substeps)
    TwoDTurb.updatevars!(vs, g)
    maxq = maximum(abs.(vs.q))

    if plots
      imshow(vs.q)
      pause(0.01)
    end

    log1 = @sprintf("τ: %.3f s, step: %d, t*q0: %.1e, max q: %.3e, ",
      time()-starttime, ts.step, vs.t*q0, maxq)

    log2 = @sprintf("ΔE: %.3f, CFL: %.3f",
      energy(vs, g)/E0, maximum([vs.U; vs.V])*ts.dt/g.dx)

    println(log1*log2)

  end

  @printf("... done.")

  return vs.q
end







end
# E N D   T W O D T U R B >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
