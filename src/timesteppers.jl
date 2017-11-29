__precompile__()




export ForwardEulerTimeStepper, FilteredForwardEulerTimeStepper,
       AB3TimeStepper,
       RK4TimeStepper,
       ETDRK4TimeStepper, FilteredETDRK4TimeStepper

export stepforward!




# Looping stepforward function ------------------------------------------------
function stepforward!(v::AbstractVars, ts::AbstractTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid; nsteps=1)
  for step = 1:nsteps
    stepforward!(v, ts, eq, p, g)
  end
  nothing
end

function stepforward!(prob::Problem; nsteps=1)
  for step = 1:nsteps
    stepforward!(prob.vars, prob.ts, prob.eqn, prob.params, prob.grid)
    prob.t = prob.vars.t
    prob.step = prob.ts.step
  end
  nothing
end

function stepforward!(prob::Problem, diags::AbstractArray;
  nsteps=1)

  # Initialize diagnostics for speed
  for diag in diags
    newnum = ceil(Int, (diag.count+nsteps)/diag.freq)
    if newnum > diag.num
      resize!(diag, newnum)
    end
  end

  for step = 1:nsteps
    stepforward!(prob.vars, prob.ts, prob.eqn, prob.params, prob.grid)
    prob.t = prob.vars.t
    prob.step = prob.ts.step
    for diag in diags
      if (prob.step+1) % diag.freq == 0.0
        increment!(diag)
      end
    end

  end
  nothing
end

function stepforward!(prob::Problem, diag::AbstractDiagnostic; nsteps=1)
  stepforward!(prob, [diag]; nsteps=nsteps)
end




# Utilities -------------------------------------------------------------------
"""
Returns an filter with an exponentially-decaying profile that, when multiplied,
removes high-wavenumber content from a spectrum on the grid g. The filter is
constructed to transition from unity inside a circle in wavenumber space with
radius innerK to zero outside a circle with radius outerK. The filter has 
size=(nk, nl, nvars), where nk may be equal to nx or nx/2+1 depending on 
whether the associated variable is real or complex.
"""
function makefilter(filtersize, g::TwoDGrid; 
  order=4.0, innerK=0.65, outerK=1.0)

  if filtersize[1] == g.nkr
    realvars = true
  else
    realvars = false
  end

  # Get decay rate for filter
  decay = 15.0*log(10.0) / (outerK-innerK)^order

  # Non-dimensional square wavenumbers
  if realvars
    KK = sqrt.( (g.Kr*g.dx/π).^2 + (g.Lr*g.dy/π).^2 )
  else
    KK = sqrt.( (g.K*g.dx/π).^2  + (g.L*g.dy/π).^2  )
  end

  filter = exp.( -decay*(KK-innerK).^order ) # filter is Complex{Float64} array
  filter[ real.(KK) .< innerK ] = 1.0

  # Return filter broadcasted to the specified shape
  ones(Complex{Float64}, filtersize).*filter
end




"""
Return the ETDRK4 coefficients that correspond to the (diagonal) linear 
coefficients of the equation, LC, and the time step dt. The coefficients are 
calculated by contour integration using ncirc points over a circle of radius
rcirc in complex space.
"""
function get_etdrk4coeffs(dt::Float64, LC; ncirc=32, rcirc=1.0)

  # Find the circle's shape, which is (1, 1, ..., ncirc)
  circshape = [1 for dim=1:ndims(LC)]
  push!(circshape, ncirc)

  # Complex values on a circle in complex space with dimension (1, ..., ncirc)
  circ = reshape(rcirc*exp.(2π*im*(0.5:1.0:(ncirc-0.5))/ncirc), 
    Tuple(circshape))

  # Construct intermediate vars
  zc = broadcast(+, dt*LC, circ)
  zdims = ndims(zc)

  # Four coefficients: ζ, α, β, Γ
  ζ = dt*squeeze(mean( (exp.(zc/2)-1)./zc, zdims), zdims)

  α = dt*squeeze(mean(
    ( -4  .- zc .+ exp.(zc).*(4 .- 3zc .+ zc.^2) )./zc.^3, zdims), zdims)
  β = dt*squeeze(mean(
    (  2  .+ zc .+ exp.(zc).*(-2 .+ zc)          )./zc.^3, zdims), zdims)
  Γ = dt*squeeze(mean(
    ( -4 .- 3zc .- zc.^2 .+ exp.(zc).*(4 .- zc)  )./zc.^3, zdims), zdims)

  ζ, α, β, Γ
end




# Forward Euler ---------------------------------------------------------------
"""
The Forward Euler method is the simplest time-stepping method in the books. 
It's fully explicit and 1st-order accurate.
"""
type ForwardEulerTimeStepper{dim} <: AbstractTimeStepper
  step::Int
  dt::Float64
  NL::Array{Complex{Float64}, dim}    # Nonlinear term
end

function ForwardEulerTimeStepper(dt::Float64, sol::Array{Complex{Float64}, 2})
  ForwardEulerTimeStepper{ndims(LC)}(0, dt, zeros(sol))
end

function ForwardEulerTimeStepper(dt::Float64, v::AbstractVars)
  ForwardEulerTimeStepper{ndims(NL)}(dt, v.sol)
end




function stepforward!(v::AbstractVars, ts::ForwardEulerTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  eq.calcNL!(ts.NL, v.sol, v.t, v, p, g)

  @. v.sol += ts.dt*(ts.NL + eq.LC*v.sol)
  v.t += ts.dt
  ts.step += 1

  nothing
end




# Filtered Forward Euler ------------------------------------------------------
"""
The 'filtered Forward Euler' scheme uses a simple linear filter to remove 
high-wavenumber variance from the solution in each time-step. The filter is 
designed is unity in low wavenumber regions and falls off to zero at
high wavenumbers.
"""
type FilteredForwardEulerTimeStepper{dim} <: AbstractTimeStepper
  step::Int
  dt::Float64
  NL::Array{Complex{Float64}, dim}        # Nonlinear term
  filter::Array{Complex{Float64}, dim}    # Filter for solution
end

function FilteredForwardEulerTimeStepper(dt::Float64, g::AbstractGrid,
  sol::AbstractArray; filterorder=4.0, innerfilterK=0.65, outerfilterK=0.95)

  filter = makefilter(size(sol), g; order=filterorder, innerK=innerfilterK,
    outerK=outerfilterK)

  FilteredForwardEulerTimeStepper{ndims(NL)}(0, dt, zeros(sol), filter)
end

function FilteredForwardEulerTimeStepper(dt::Float64, g::AbstractGrid,
  v::AbstractVars; filterorder=4.0, innerfilterK=0.65, outerfilterK=0.95)

  FilteredForwardEulerTimeStepper(dt, g, v.sol; filterorder=filterorder,
    innerfilterK=innerfilterK, outerfilterK=outerfilterK)
end



function stepforward!(v::AbstractVars, ts::FilteredForwardEulerTimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  eq.calcNL!(ts.NL, v.sol, v.t, v, p, g)

  @. v.sol = ts.filter*(v.sol + ts.dt*(ts.NL + eq.LC.*v.sol))
  v.t += ts.dt
  ts.step += 1

  nothing
end








# ETDRK4 ----------------------------------------------------------------------
"""
The Exponential Time Differencing RK4 is the Rolls-Royce of time-stepping. 
It solves the linear part of the equation exactly over each time-step, 
while integrating the nonlinear part with 4th-order RK4-like scheme.
"""
type ETDRK4TimeStepper{dim} <: AbstractTimeStepper
  step::Int
  dt::Float64

  LC::Array{Complex{Float64}, dim}          # Linear coefficient

  # ETDRK4 coefficents
  ζ::Array{Complex{Float64}, dim}
  α::Array{Complex{Float64}, dim}
  β::Array{Complex{Float64}, dim}
  Γ::Array{Complex{Float64}, dim}
  expLCdt::Array{Complex{Float64}, dim}     # Precomputed exp(LC*dt)
  expLC½dt::Array{Complex{Float64}, dim}    # Precomputed exp(LC*dt/2)

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64
  sol₁::Array{Complex{Float64}, dim}
  sol₂::Array{Complex{Float64}, dim}
  NL₁::Array{Complex{Float64}, dim}
  NL₂::Array{Complex{Float64}, dim}
  NL₃::Array{Complex{Float64}, dim}
  NL₄::Array{Complex{Float64}, dim}
end

"""
Construct an ETDRK4TimeStepper for timestep dt and array of linear coefficients
LC. The ETDRK4 coefficients are calculated with contour integration over ncirc
points along a circle with radius rcirc.
"""
function ETDRK4TimeStepper(dt::Float64, LC::AbstractArray; ncirc=32, rcirc=1.0)

   expLCdt = exp.(dt*LC)
  expLC½dt = exp.(dt*LC/2.0)

  ζ, α, β, Γ = get_etdrk4coeffs(dt, LC; ncirc=ncirc, rcirc=rcirc)

  ETDRK4TimeStepper{ndims(LC)}(0, dt, LC, ζ, α, β, Γ, expLCdt, expLC½dt, 0.0,
    zeros(LC), zeros(LC), zeros(LC), zeros(LC), zeros(LC), zeros(LC))
end




function stepforward!(v::AbstractVars, ts::ETDRK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  # Substep 1
  eq.calcNL!(ts.NL₁, v.sol, v.t, v, p, g)
  @. ts.sol₁ = ts.expLC½dt*v.sol + ts.ζ*ts.NL₁

  # Substep 2
  ts.ti = v.t + 0.5*ts.dt
  eq.calcNL!(ts.NL₂, ts.sol₁, ts.ti, v, p, g)
  @. ts.sol₂ = ts.expLC½dt*v.sol + ts.ζ*ts.NL₂

  # Substep 3
  eq.calcNL!(ts.NL₃, ts.sol₂, ts.ti, v, p, g)
  @. ts.sol₂ = ts.expLC½dt*ts.sol₁ + ts.ζ*(2.0*ts.NL₃ - ts.NL₁)

  # Substep 4
  ts.ti = v.t + ts.dt
  eq.calcNL!(ts.NL₄, ts.sol₂, ts.ti, v, p, g)

  # Update
  @. v.sol = (ts.expLCdt*v.sol +     ts.α * ts.NL₁
                               + 2.0*ts.β * (ts.NL₂ + ts.NL₃)
                               +     ts.Γ * ts.NL₄ )

  v.t += ts.dt
  ts.step += 1

  nothing
end




# Filtered ETDRK4 --------------------------------------------------------------
"""
The Exponential Time Differencing RK4 is the Rolls-Royce of time-stepping. 
It solves the linear part of the equation exactly over each time-step, 
while integrating the nonlinear part with 4th-order RK4-like scheme.
"""
type FilteredETDRK4TimeStepper{dim} <: AbstractTimeStepper
  step::Int
  dt::Float64

  LC::Array{Complex{Float64}, dim}          # Linear coefficient

  # ETDRK4 coefficents
  ζ::Array{Complex{Float64}, dim}
  α::Array{Complex{Float64}, dim}
  β::Array{Complex{Float64}, dim}
  Γ::Array{Complex{Float64}, dim}
  expLCdt::Array{Complex{Float64}, dim}     # Precomputed exp(LC*dt)
  expLC½dt::Array{Complex{Float64}, dim}    # Precomputed exp(LC*dt/2)

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64
  sol₁::Array{Complex{Float64}, dim}
  sol₂::Array{Complex{Float64}, dim}
  NL₁::Array{Complex{Float64}, dim}
  NL₂::Array{Complex{Float64}, dim}
  NL₃::Array{Complex{Float64}, dim}
  NL₄::Array{Complex{Float64}, dim}
  filter::Array{Complex{Float64}, dim}    # Filter for solution
end


function FilteredETDRK4TimeStepper(dt::Float64, LC::AbstractArray,
    g::AbstractGrid; filterorder=4.0, innerfilterK=0.65, outerfilterK=0.95,
    ncirc=32, rcirc=1.0)

  expLCdt  = exp.(dt*LC)
  expLC½dt = exp.(dt*LC/2.0)

  ζ, α, β, Γ = get_etdrk4coeffs(dt, LC; ncirc=ncirc, rcirc=rcirc)

  filter = makefilter(size(LC), g; order=filterorder, innerK=innerfilterK,
    outerK=outerfilterK)

  FilteredETDRK4TimeStepper{ndims(LC)}(
    0, dt, LC, ζ, α, β, Γ, expLCdt, expLC½dt, 0.0,
    zeros(LC), zeros(LC), zeros(LC), zeros(LC), zeros(LC), zeros(LC), filter)
end




function stepforward!(v::AbstractVars, ts::FilteredETDRK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  # Substep 1
  eq.calcNL!(ts.NL₁, v.sol, v.t, v, p, g)
  @. ts.sol₁ = ts.expLC½dt*v.sol + ts.ζ*ts.NL₁

  # Substep 2
  ts.ti = v.t + 0.5*ts.dt
  eq.calcNL!(ts.NL₂, ts.sol₁, ts.ti, v, p, g)
  @. ts.sol₂ = ts.expLC½dt*v.sol + ts.ζ*ts.NL₂

  # Substep 3
  eq.calcNL!(ts.NL₃, ts.sol₂, ts.ti, v, p, g)
  @. ts.sol₂ = ts.expLC½dt*ts.sol₁ + ts.ζ*(2.0*ts.NL₃ - ts.NL₁)

  # Substep 4
  ts.ti = v.t + ts.dt
  eq.calcNL!(ts.NL₄, ts.sol₂, ts.ti, v, p, g)

  # Update
  @. v.sol = ts.filter*(ts.expLCdt.*v.sol +     ts.α * ts.NL₁
                                          + 2.0*ts.β * (ts.NL₂ + ts.NL₃)
                                          +     ts.Γ * ts.NL₄ )
  v.t += ts.dt
  ts.step += 1

end




# Dual ETDRK4 -----------------------------------------------------------------
"""
The DualETDRK4TimeStepper is a convenience type that wraps two individual
ETDRK4TimeSteppers into a single object, for use when an equation has both 
real variables and imaginary variables that require the separate
ETDRK4TimeStepper objects.
"""
type DualETDRK4TimeStepper{dimc, dimr} <: AbstractTimeStepper
  step::Int
  dt::Float64
  c::ETDRK4TimeStepper{dimc}
  r::ETDRK4TimeStepper{dimr}
end

function ETDRK4TimeStepper(dt::Float64, LCc::AbstractArray, LCr::AbstractArray)
  c = ETDRK4TimeStepper(dt::Float64, LCc)
  r = ETDRK4TimeStepper(dt::Float64, LCr)
  DualETDRK4TimeStepper{ndims(LCc), ndims(LCr)}(0, dt, c, r)
end




function stepforward!(v::AbstractVars, ts::DualETDRK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  # Substep 1
  eq.calcNL!(ts.c.NL₁, ts.r.NL₁, v.solc, v.solr, v.t, v, p, g)
  @. ts.c.sol₁ = ts.c.expLC½dt*v.solc + ts.c.ζ*ts.c.NL₁
  @. ts.r.sol₁ = ts.r.expLC½dt*v.solr + ts.r.ζ*ts.r.NL₁

  # Substep 2
  ts.c.ti = v.t + 0.5*ts.c.dt
  eq.calcNL!(ts.c.NL₂, ts.r.NL₂, ts.c.sol₁, ts.r.sol₁, ts.c.ti, v, p, g)
  @. ts.c.sol₂ = ts.c.expLC½dt*v.solc + ts.c.ζ*ts.c.NL₂
  @. ts.r.sol₂ = ts.r.expLC½dt*v.solr + ts.r.ζ*ts.r.NL₂

  # Substep 3
  eq.calcNL!(ts.c.NL₃, ts.r.NL₃, ts.c.sol₂, ts.r.sol₂, ts.c.ti, v, p, g)
  @. ts.c.sol₂ = (ts.c.expLC½dt*ts.c.sol₁ + ts.c.ζ*(2.0*ts.c.NL₃ - ts.c.NL₁))
  @. ts.r.sol₂ = (ts.r.expLC½dt*ts.r.sol₁ + ts.r.ζ*(2.0*ts.r.NL₃ - ts.r.NL₁))
    

  # Substep 4
  ts.c.ti = v.t + ts.c.dt
  eq.calcNL!(ts.c.NL₄, ts.r.NL₄, ts.c.sol₂, ts.r.sol₂, ts.c.ti, v, p, g)

  # Update
  @. v.solc = (ts.c.expLCdt*v.solc +     ts.c.α * ts.c.NL₁
                                   + 2.0*ts.c.β * (ts.c.NL₂ + ts.c.NL₃)
                                   +     ts.c.Γ * ts.c.NL₄ )

  @. v.solr = (ts.r.expLCdt*v.solr +     ts.r.α * ts.r.NL₁
                                   + 2.0*ts.r.β * (ts.r.NL₂ + ts.r.NL₃)
                                   +     ts.r.Γ * ts.r.NL₄ )

  v.t += ts.dt
  ts.step += 1
  ts.c.step += 1
  ts.r.step += 1

  nothing
end




# RK4 -------------------------------------------------------------------------
"""
RK4 is the classical explicit 4th-order Runge-Kutta time-stepping
method. It uses a series of substeps/estimators to achieve 4th-order
accuracy over each individual time-step, at the cost of requiring
relatively more evaluations of the nonlinear right hand side.
It is described, among other places, in Bewley's Numerical
Renaissance.
"""
type RK4TimeStepper <: AbstractTimeStepper
  step::Int
  dt::Float64

  # Intermediate times, solutions, and nonlinear evaluations
  ti::Float64
  sol₁::Array{Complex{Float64}, 2}
  RHS₁::Array{Complex{Float64}, 2}
  RHS₂::Array{Complex{Float64}, 2}
  RHS₃::Array{Complex{Float64}, 2}
  RHS₄::Array{Complex{Float64}, 2}
end

function RK4TimeStepper(dt::Float64, v::AbstractVars)
  RK4TimeStepper(dt, v.sol)
end

function RK4TimeStepper(dt::Float64, LC::Array{Complex{Float64}, 2})
  RK4TimeStepper(0, dt, 0.0, 
    zeros(LC), zeros(LC), zeros(LC), zeros(LC), zeros(LC))
end




function stepforward!(v::AbstractVars, ts::RK4TimeStepper,
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid)

  eq.calcNL!(ts.RHS₁, v.sol, v.t, v, p, g)
  @. ts.RHS₁ += eq.LC*v.sol

  # Substep 1
  ts.ti = v.t + 0.5*ts.dt
  @. ts.sol₁ = v.sol + 0.5*ts.dt*ts.RHS₁
  eq.calcNL!(ts.RHS₂, ts.sol₁, v.t, v, p, g)
  @. ts.RHS₂ += eq.LC*ts.sol₁

  # Substep 2
  @. ts.sol₁ = v.sol + 0.5*ts.dt*ts.RHS₂
  eq.calcNL!(ts.RHS₃, ts.sol₁, v.t, v, p, g)
  @. ts.RHS₃ += eq.LC*ts.sol₁

  # Substep 3
  ts.ti = v.t + ts.dt
  @. ts.sol₁ = v.sol + ts.dt*ts.RHS₃
  eq.calcNL!(ts.RHS₄, ts.sol₁, v.t, v, p, g)
  @. ts.RHS₄ += eq.LC*ts.sol₁

  # Substep 4 and final step
  @. v.sol += ts.dt*( ts.RHS₁/6.0 + ts.RHS₂/3.0
                    + ts.RHS₃/3.0 + ts.RHS₄/6.0 )

  v.t += ts.dt
  ts.step += 1
end








# AB3 -------------------------------------------------------------------------
"""
3rd order Adams-Bashforth time stepping is an explicit scheme that uses
solutions from two previous time-steps to achieve 3rd order accuracy.
"""
type AB3TimeStepper <: AbstractTimeStepper
  step::Int
  dt::Float64
  RHS::Array{Complex{Float64}, 2}     # RHS at current step
  RHS₋₁::Array{Complex{Float64}, 2}   # RHS at one step previous, step-1
  RHS₋₂::Array{Complex{Float64}, 2}   # RHS two steps previous at step-2
  _initstep::Int                      # Counter for initialization steps
end

function AB3TimeStepper(dt::Float64, sol::Array{Complex{Float64}, 2})
  AB3TimeStepper(0, dt, zeros(sol), zeros(sol), zeros(sol), 0)
end

function AB3TimeStepper(dt::Float64, v::AbstractVars)
  AB3TimeStepper(dt, v.sol)
end




"""
Special stepforward function for Adams-Bashforth that accounts for
initialization.
"""
function stepforward!(v::AbstractVars, ts::AB3TimeStepper, 
  eq::AbstractEquation, p::AbstractParams, g::AbstractGrid;
  nsteps=1, initialize=false)

  # Reset _initsteps if called for 
  ts._initstep = initialize ? 0 : ts._initstep
  step₀ = 1

  while ts._initstep < 2 # take forward Euler steps to initialize AB3

    eq.calcNL!(ts.RHS, v.sol, v.t, v, p, g)
    @. ts.RHS += eq.LC*v.sol

    @. v.sol += ts.dt*ts.RHS
    v.t += ts.dt
    ts.step += 1

    ts._initstep += 1

    # Update stored RHS's
    @. ts.RHS₋₂ = ts.RHS₋₁
    @. ts.RHS₋₁ = ts.RHS
    
    step₀ += 1
  end

  # Loop
  for step = step₀:nsteps
    eq.calcNL!(ts.RHS, v.sol, v.t, v, p, g)
    @. ts.RHS += eq.LC*v.sol

    @. v.sol += ts.dt*(
      23.0/12.0 * ts.RHS - 16.0/12.0 * ts.RHS₋₁ + 5.0/12.0 * ts.RHS₋₂ )

    v.t += ts.dt
    ts.step += 1

    # Update stored RHS's
    @. ts.RHS₋₂ = ts.RHS₋₁
    @. ts.RHS₋₁ = ts.RHS
  end

end
