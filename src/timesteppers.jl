export ForwardEulerTimeStepper, FilteredForwardEulerTimeStepper,
       RK4TimeStepper, FilteredRK4TimeStepper, 
       ETDRK4TimeStepper, FilteredETDRK4TimeStepper,
       AB3TimeStepper

export stepforward!


# Looping stepforward function
"""
    stepforward!(prob)

Step forward the Problem prob for one timestep.
"""
function stepforward!(prob::Problem)
    stepforward!(prob.state, prob.ts, prob.eqn, prob.vars, prob.params, prob.grid)
    prob.t = prob.state.t
    prob.step = prob.state.step
end

"""
    stepforward!(prob, nsteps)

Step forward the problem 'prob' for 'nsteps'.
"""
stepforward!(prob::Problem, nsteps) = for step = 1:nsteps; stepforward!(prob); end

"""
    stepforward!(prob, diags, nsteps)

Step forward the problem prob for nsteps while calculating the 
diagnostics in diags.
"""
function stepforward!(prob::Problem, diags::AbstractArray, nsteps)
  # Initialize diagnostics for speed
  for diag in diags
    newnum = ceil(Int, (diag.count+nsteps)/diag.freq)
    if newnum > diag.num
      warn("Resizing diags before stepping forward...")
      resize!(diag, newnum)
    end
  end

  for step = 1:nsteps
    stepforward!(prob)
    for diag in diags
      if (prob.step+1) % diag.freq == 0.0
        increment!(diag)
      end
    end
  end
  nothing
end


# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
# Timestepper utilities -------------------------------------------------------
# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
"""
    getetdcoeffs(dt, LC; ncirc=32, rcirc=1)

Calculate ETDRK4 coefficients associated with the (diagonal) linear coefficient
LC by integrating over a small circle in complex space.

Note: arbitrary-precision arithmetic might provide a more robust method for 
calculating these coefficients.
"""
function getetdcoeffs(dt, LC; ncirc=32, rcirc=eltype(LC)(1))

  shape = Tuple(cat(1, ncirc, ones(Int, ndims(LC))))

  circ = zeros(cxeltype(LC), shape)
  circ .= rcirc * exp.(2π*im/ncirc*(0.5:1:(ncirc-0.5)))
  circ = permutedims(circ, ndims(circ):-1:1)

  zc = dt*LC .+ circ
  M = ndims(LC)+1

  # Four coefficients: ζ, α, β, Γ
  ζc = @.          ( exp(zc/2)-1 ) / zc
  αc = @. ( -4 - zc + exp(zc)*(4 - 3zc + zc^2) ) / zc^3
  βc = @.    ( 2  + zc + exp(zc)*(-2 + zc) ) / zc^3
  Γc = @. ( -4 - 3zc - zc^2 + exp(zc)*(4 - zc) ) / zc^3

  if eltype(LC) <: Real
    ζ = dt*real.(squeeze(mean(ζc, M), M))
    α = dt*real.(squeeze(mean(αc, M), M))
    β = dt*real.(squeeze(mean(βc, M), M))
    Γ = dt*real.(squeeze(mean(Γc, M), M))
  else
    ζ = dt*squeeze(mean(ζc, M), M)
    α = dt*squeeze(mean(αc, M), M)
    β = dt*squeeze(mean(βc, M), M)
    Γ = dt*squeeze(mean(Γc, M), M)
  end

  ζ, α, β, Γ
end


# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 
# Timesteppers
# ----------------------------------------------------------------------------- 
# ----------------------------------------------------------------------------- 


# Forward Euler ---------------------------------------------------------------
# The simplest time-stepping method in the books. Explicit and 1st-order
# accurate.

struct ForwardEulerTimeStepper{Tdt,Tsol,dim} <: AbstractTimeStepper
  dt::Tdt
  N::Array{Tsol,dim}      # Explicit linear and nonlinear terms
end

struct FilteredForwardEulerTimeStepper{Tdt,Tsol,Tfilt,dim} <: AbstractTimeStepper
  dt::Tdt
  N::Array{Tsol,dim}      # Explicit linear and nonlinear terms
  filter::Array{Tfilt,dim}  # Filter for solution
end

function ForwardEulerTimeStepper(dt, LC, soltype::DataType=cxeltype(LC))
  ForwardEulerTimeStepper(dt, zeros(soltype, size(LC)))
end

function FilteredForwardEulerTimeStepper(dt, LC, g, soltype::DataType=cxeltype(LC); filterkwargs...) 
  filter = makefilter(g, real(soltype), size(LC); filterkwargs...)
  FilteredForwardEulerTimeStepper(dt, zeros(soltype, size(LC)), filter)
end

function stepforward!(s, ts::ForwardEulerTimeStepper, eq, v, p, g)
  eq.calcN!(ts.N, s.sol, s.t, s, v, p, g)
  @. s.sol += ts.dt*(ts.N + eq.LC*s.sol)
  s.t += ts.dt
  s.step += 1
  nothing
end

function stepforward!(s, ts::FilteredForwardEulerTimeStepper, eq, v, p, g)
  eq.calcN!(ts.N, s.sol, s.t, s, v, p, g)
  @. s.sol = ts.filter*(s.sol + ts.dt*(ts.N + eq.LC.*s.sol) )
  s.t += ts.dt
  s.step += 1
  nothing
end


# RK4 -------------------------------------------------------------------------
# RK4 is the classical explicit 4th-order Runge-Kutta time-stepping
# method. It uses a series of substeps/estimators to achieve 4th-order
# accuracy over each individual time-step, at the cost of requiring
# relatively more evaluations of the nonlinear right hand side.
# It is described, among other places, in Bewley's Numerical
# Renaissance.

abstract type AbstractRK4TimeStepper <: AbstractTimeStepper end
abstract type AbstractDualRK4TimeStepper <: AbstractTimeStepper end

struct RK4TimeStepper{Tdt,Tsol,dim} <: AbstractRK4TimeStepper
  dt::Tdt
  sol₁::Array{Tsol,dim}
  RHS₁::Array{Tsol,dim}
  RHS₂::Array{Tsol,dim}
  RHS₃::Array{Tsol,dim}
  RHS₄::Array{Tsol,dim}
end

struct FilteredRK4TimeStepper{Tdt,Tsol,Tfilt,dim} <: AbstractRK4TimeStepper
  dt::Tdt
  sol₁::Array{Tsol,dim}
  RHS₁::Array{Tsol,dim}
  RHS₂::Array{Tsol,dim}
  RHS₃::Array{Tsol,dim}
  RHS₄::Array{Tsol,dim}
  filter::Array{Tfilt,dim}
end

struct DualRK4TimeStepper{Tdt,Tc,Tr,dimc,dimr} <: AbstractDualRK4TimeStepper
  dt::Tdt
  c::RK4TimeStepper{Tdt,Tc,dimc}
  r::RK4TimeStepper{Tdt,Tr,dimr}
end

struct DualFilteredRK4TimeStepper{Tdt,Tc,Tr,dimc,dimr} <: AbstractDualRK4TimeStepper
  dt::Tdt
  c::FilteredRK4TimeStepper{Tdt,Tc,dimc}
  r::FilteredRK4TimeStepper{Tdt,Tr,dimr}
end

function RK4TimeStepper(dt, LC, soltype::DataType=cxeltype(LC))
  @createarrays soltype size(LC) sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  RK4TimeStepper(dt, sol₁, RHS₁, RHS₂, RHS₃, RHS₄)
end

function FilteredRK4TimeStepper(dt, LC, g::AbstractGrid, soltype::DataType=cxeltype(LC); filterkwargs...)
  @createarrays soltype size(LC) sol₁ RHS₁ RHS₂ RHS₃ RHS₄
  filter = makefilter(g, real(soltype), size(LC); filterkwargs...)
  FilteredRK4TimeStepper(dt, sol₁, RHS₁, RHS₂, RHS₃, RHS₄, filter)
end

function RK4TimeStepper(dt, LCc, LCr, solctype::DataType=cxeltype(LCc), solrtype::DataType=cxeltype(LCr))
  c = RK4TimeStepper(dt, LCc, solctype)
  r = RK4TimeStepper(dt, LCr, solrtype)
  DualRK4TimeStepper(dt, c, r)
end

function FilteredRK4TimeStepper(dt, LCc, LCr, g::AbstractGrid, solctype::DataType=cxeltype(LCc), 
  solrtype::DataType=cxeltype(LCr); filterkwargs...)
  c = FilteredRK4TimeStepper(dt, LCc, g, solctype)
  r = FilteredRK4TimeStepper(dt, LCr, g, solrtype)
  DualFilteredRK4TimeStepper(dt, c, r)
end

function stepRK4!(s, ts::AbstractRK4TimeStepper, eq, v, p, g)
  eq.calcN!(ts.RHS₁, s.sol, s.t, s, v, p, g)
  @. ts.RHS₁ += eq.LC*s.sol
  # Substep 1
  t2 = s.t + 0.5*ts.dt
  @. ts.sol₁ = s.sol + 0.5*ts.dt*ts.RHS₁
  eq.calcN!(ts.RHS₂, ts.sol₁, t2, s, v, p, g)
  @. ts.RHS₂ += eq.LC*ts.sol₁
  # Substep 2
  @. ts.sol₁ = s.sol + 0.5*ts.dt*ts.RHS₂
  eq.calcN!(ts.RHS₃, ts.sol₁, t2, s, v, p, g)
  @. ts.RHS₃ += eq.LC*ts.sol₁
  # Substep 3
  t3 = s.t + ts.dt
  @. ts.sol₁ = s.sol + ts.dt*ts.RHS₃
  eq.calcN!(ts.RHS₄, ts.sol₁, t3, s, v, p, g)
  @. ts.RHS₄ += eq.LC*ts.sol₁

  # Substep 4 and final step
  @. s.sol += ts.dt*(   1.0/6.0*ts.RHS₁ + 1.0/3.0*ts.RHS₂
                      + 1.0/3.0*ts.RHS₃ + 1.0/6.0*ts.RHS₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end

function stepRK4!(s, ts::AbstractDualRK4TimeStepper, eq, v, p, g)
  eq.calcN!(ts.c.RHS₁, ts.r.RHS₁, s.solc, s.solr, s.t, s, v, p, g)
  @. ts.c.RHS₁ += eq.LCc*s.solc
  @. ts.r.RHS₁ += eq.LCr*s.solr
  # Substep 1
  t2 = s.t + 0.5*ts.dt
  @. ts.c.sol₁ = s.solc + 0.5*ts.dt*ts.c.RHS₁
  @. ts.r.sol₁ = s.solr + 0.5*ts.r.dt*ts.r.RHS₁
  eq.calcN!(ts.c.RHS₂, ts.r.RHS₂, ts.c.sol₁, ts.r.sol₁, t2, s, v, p, g)
  @. ts.c.RHS₂ += eq.LCc*ts.c.sol₁
  @. ts.r.RHS₂ += eq.LCr*ts.r.sol₁
  # Substep 2
  @. ts.c.sol₁ = s.solc + 0.5*ts.dt*ts.c.RHS₂
  @. ts.r.sol₁ = s.solr + 0.5*ts.dt*ts.r.RHS₂
  eq.calcN!(ts.c.RHS₃, ts.r.RHS₃, ts.c.sol₁, ts.r.sol₁, t2, s, v, p, g)
  @. ts.c.RHS₃ += eq.LCc*ts.c.sol₁
  @. ts.r.RHS₃ += eq.LCr*ts.r.sol₁
  # Substep 3
  t3 = s.t + ts.dt
  @. ts.c.sol₁ = s.solc + ts.dt*ts.c.RHS₃
  @. ts.r.sol₁ = s.solr + ts.dt*ts.r.RHS₃
  eq.calcN!(ts.c.RHS₄, ts.r.RHS₄, ts.c.sol₁, ts.r.sol₁, t3, s, v, p, g)
  @. ts.c.RHS₄ += eq.LCc*ts.c.sol₁
  @. ts.r.RHS₄ += eq.LCr*ts.r.sol₁

  # Substep 4 and final step
  @. s.solc += ts.dt*(   1/6*ts.c.RHS₁ + 1/3*ts.c.RHS₂
                       + 1/3*ts.c.RHS₃ + 1/6*ts.c.RHS₄ )

  @. s.solr += ts.dt*(   1/6*ts.r.RHS₁ + 1/3*ts.r.RHS₂
                       + 1/3*ts.r.RHS₃ + 1/6*ts.r.RHS₄ )

  s.t += ts.dt
  s.step += 1

  nothing
end

stepforward!(s, ts::RK4TimeStepper, eq, v, p, g) = stepRK4!(s, ts, eq, v, p, g) 
stepforward!(s, ts::DualRK4TimeStepper, eq, v, p, g) = stepRK4!(s, ts, eq, v, p, g)  

function stepforward!(s, ts::FilteredRK4TimeStepper, eq, v, p, g)
  stepRK4!(s, ts, eq, v, p, g)
  @. s.sol *= ts.filter
  nothing
end

function stepforward!(s, ts::DualFilteredRK4TimeStepper, eq, v, p, g)
  stepRK4!(s, ts, eq, v, p, g) 
  @. s.solc *= ts.c.filter
  @. s.solr *= ts.r.filter
  nothing
end


# ETDRK4 ----------------------------------------------------------------------
# The Rolls-Royce of time-stepping. Exact treatment of the implicit linear part
# of the equation, explicit and 4th-order accurate integration of nonlinear
# parts of equation.

abstract type AbstractETDRK4TimeStepper <: AbstractTimeStepper end
abstract type AbstractDualETDRK4TimeStepper <: AbstractTimeStepper end

struct ETDRK4TimeStepper{Tdt,TLC,Tsol,dim} <: AbstractETDRK4TimeStepper
  dt::Tdt
  LC::Array{TLC,dim}          # Linear coefficient
  ζ::Array{TLC,dim}           # ETDRK4 coefficents
  α::Array{TLC,dim}           # ...
  β::Array{TLC,dim}           # ...
  Γ::Array{TLC,dim}
  expLCdt::Array{TLC,dim}     # Precomputed exp(LC*dt)
  expLCdt2::Array{TLC,dim}    # Precomputed exp(LC*dt/2)
  sol₁::Array{Tsol,dim}
  sol₂::Array{Tsol,dim}
  N₁::Array{Tsol,dim}
  N₂::Array{Tsol,dim}
  N₃::Array{Tsol,dim}
  N₄::Array{Tsol,dim}
end

struct FilteredETDRK4TimeStepper{Tdt,TLC,Tsol,Tfilt,dim} <: AbstractETDRK4TimeStepper
  dt::Tdt
  LC::Array{TLC,dim}          # Linear coefficient
  ζ::Array{TLC,dim}
  α::Array{TLC,dim}
  β::Array{TLC,dim}
  Γ::Array{TLC,dim}
  expLCdt::Array{TLC,dim}     # Precomputed exp(LC*dt)
  expLCdt2::Array{TLC,dim}    # Precomputed exp(LC*dt/2)
  sol₁::Array{Tsol,dim}
  sol₂::Array{Tsol,dim}
  N₁::Array{Tsol,dim}
  N₂::Array{Tsol,dim}
  N₃::Array{Tsol,dim}
  N₄::Array{Tsol,dim}
  filter::Array{Tfilt,dim}    # Filter for solution
end

struct DualETDRK4TimeStepper{Tdt,TLCc,TLCr,Tc,Tr,dimc,dimr} <: AbstractDualETDRK4TimeStepper
  dt::Tdt
  c::ETDRK4TimeStepper{Tdt,TLCc,Tc,dimc}
  r::ETDRK4TimeStepper{Tdt,TLCr,Tr,dimr}
end

struct DualFilteredETDRK4TimeStepper{Tdt,TLCc,TLCr,Tc,Tr,dimc,dimr} <: AbstractTimeStepper
  dt::Tdt
  c::FilteredETDRK4TimeStepper{Tdt,TLCc,Tc,dimc}
  r::FilteredETDRK4TimeStepper{Tdt,TLCr,Tr,dimr}
end

function ETDRK4TimeStepper(dt, LC, soltype::DataType=cxeltype(LC))
   expLCdt = exp.(dt*LC)
  expLCdt2 = exp.(dt*LC/2)
  ζ, α, β, Γ = getetdcoeffs(dt, LC)
  @createarrays soltype size(LC) sol₁ sol₂ N₁ N₂ N₃ N₄
  ETDRK4TimeStepper(dt, LC, ζ, α, β, Γ, expLCdt, expLCdt2, sol₁, sol₂, N₁, N₂, N₃, N₄)
end

function FilteredETDRK4TimeStepper(dt, LC, g::AbstractGrid, soltype::DataType=cxeltype(LC); filterkwargs...)
  expLCdt = exp.(dt*LC)
  expLCdt2 = exp.(0.5*dt*LC)
  ζ, α, β, Γ = getetdcoeffs(dt, LC)
  @createarrays soltype size(LC) sol₁ sol₂ N₁ N₂ N₃ N₄
  filter = makefilter(g, real(soltype), size(LC); filterkwargs...)
  FilteredETDRK4TimeStepper(dt, LC, ζ, α, β, Γ, expLCdt, expLCdt2, sol₁, sol₂, N₁, N₂, N₃, N₄, filter)
end

function ETDRK4TimeStepper(dt, LCc, LCr, solctype::DataType=cxeltype(LCc), solrtype::DataType=cxeltype(LCr))
  c = ETDRK4TimeStepper(dt, LCc, solctype)
  r = ETDRK4TimeStepper(dt, LCr, solrtype)
  DualETDRK4TimeStepper(dt, c, r)
end

function FilteredETDRK4TimeStepper(dt, LCc, LCr, g::AbstractGrid,
  solctype::DataType=cxeltype(LCc), solrtype::DataType=cxeltype(LCr); filterkwargs...)
  c = FilteredETDRK4TimeStepper(dt, LCc, g, solctype; filterkwargs...)
  r = FilteredETDRK4TimeStepper(dt, LCr, g, solrtype; filterkwargs...)
  DualFilteredETDRK4TimeStepper(dt, c, r)
end

function stepETDRK4!(s, ts::AbstractETDRK4TimeStepper, eq, v, p, g)
  # Substep 1
  eq.calcN!(ts.N₁, s.sol, s.t, s, v, p, g)
  @. ts.sol₁ = ts.expLCdt2*s.sol + ts.ζ*ts.N₁
  # Substep 2
  t2 = s.t + 0.5*ts.dt
  eq.calcN!(ts.N₂, ts.sol₁, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*s.sol + ts.ζ*ts.N₂
  # Substep 3
  eq.calcN!(ts.N₃, ts.sol₂, t2, s, v, p, g)
  @. ts.sol₂ = ts.expLCdt2*ts.sol₁ + ts.ζ*(2*ts.N₃ - ts.N₁)
  # Substep 4
  t3 = s.t + ts.dt
  eq.calcN!(ts.N₄, ts.sol₂, t3, s, v, p, g)

  # Update
  @. s.sol = (ts.expLCdt.*s.sol +   ts.α * ts.N₁
                                + 2*ts.β * (ts.N₂ + ts.N₃)
                                +   ts.Γ * ts.N₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end

function stepETDRK4!(s, ts::AbstractDualETDRK4TimeStepper, eq, v, p, g)
  # Substep 1
  eq.calcN!(ts.c.N₁, ts.r.N₁, s.solc, s.solr, s.t, s, v, p, g)
  @. ts.c.sol₁ = ts.c.expLCdt2*s.solc + ts.c.ζ*ts.c.N₁
  @. ts.r.sol₁ = ts.r.expLCdt2*s.solr + ts.r.ζ*ts.r.N₁
  # Substep 2
  t2 = s.t + 0.5*ts.c.dt
  eq.calcN!(ts.c.N₂, ts.r.N₂, ts.c.sol₁, ts.r.sol₁, t2, s, v, p, g)
  @. ts.c.sol₂ = ts.c.expLCdt2*s.solc + ts.c.ζ*ts.c.N₂
  @. ts.r.sol₂ = ts.r.expLCdt2*s.solr + ts.r.ζ*ts.r.N₂
  # Substep 3
  eq.calcN!(ts.c.N₃, ts.r.N₃, ts.c.sol₂, ts.r.sol₂, t2, s, v, p, g)
  @. ts.c.sol₂ = (ts.c.expLCdt2*ts.c.sol₁
    + ts.c.ζ*(2.0*ts.c.N₃ - ts.c.N₁))
  @. ts.r.sol₂ = (ts.r.expLCdt2*ts.r.sol₁
    + ts.r.ζ*(2.0*ts.r.N₃ - ts.r.N₁))
  # Substep 4
  t3 = s.t + ts.c.dt
  eq.calcN!(ts.c.N₄, ts.r.N₄, ts.c.sol₂, ts.r.sol₂, t3, s, v, p, g)

  # Update
  @. s.solc = (ts.c.expLCdt.*s.solc +   ts.c.α * ts.c.N₁
                                    + 2*ts.c.β * (ts.c.N₂ + ts.c.N₃)
                                    +   ts.c.Γ * ts.c.N₄ )
  @. s.solr = (ts.r.expLCdt.*s.solr +   ts.r.α * ts.r.N₁
                                    + 2*ts.r.β * (ts.r.N₂ + ts.r.N₃)
                                    +   ts.r.Γ * ts.r.N₄ )
  s.t += ts.dt
  s.step += 1

  nothing
end

stepforward!(s, ts::ETDRK4TimeStepper, eq, v, p, g) = stepETDRK4!(s, ts, eq, v, p, g)
stepforward!(s, ts::DualETDRK4TimeStepper, eq, v, p, g) = stepETDRK4!(s, ts, eq, v, p, g)

function stepforward!(s, ts::FilteredETDRK4TimeStepper, eq, v, p, g)
  stepETDRK4!(s, ts, eq, v, p, g)
  @. s.sol *= ts.filter
  nothing
end

function stepforward!(s::DualState, ts::DualFilteredETDRK4TimeStepper, eq::AbstractEquation, 
                      v::AbstractVars, p::AbstractParams, g::AbstractGrid)
  stepETDRK4!(s, ts, eq, v, p, g)
  @. s.solc *= ts.c.filter
  @. s.solr *= ts.r.filter
  nothing
end

# AB3 -------------------------------------------------------------------------
# 3rd order Adams-Bashforth time stepping is an explicit scheme that uses
# solutions from two previous time-steps to achieve 3rd order accuracy.

struct AB3TimeStepper{Tdt,Tsol,dim} <: AbstractTimeStepper
  dt::Tdt
  RHS::Array{Tsol,dim}
  RHS₋₁::Array{Tsol,dim}
  RHS₋₂::Array{Tsol,dim}
end

function AB3TimeStepper(dt, LC, soltype::DataType=cxeltype(LC))
  @createarrays soltype size(LC) RHS RHS₋₁ RHS₋₂
  AB3TimeStepper(dt, RHS, RHS₋₁, RHS₋₂)
end

function stepforward!(s, ts::AB3TimeStepper, eq, v, p, g)
  if s.step < 2                 # forward Euler steps to initialize AB3
    eq.calcN!(ts.RHS, s.sol, s.t, s, v, p, g)
    @. ts.RHS += eq.LC.*s.sol   # Add linear term to RHS

    @. s.sol += ts.dt*ts.RHS    # Update
    s.t += ts.dt                # ...
    s.step += 1                 # ...

    ts.RHS₋₂ .= ts.RHS₋₁        # Store
    ts.RHS₋₁ .= ts.RHS          # ... previous values of RHS.
  end

  # Otherwise, stepforward with 3rd order Adams Bashforth:
  eq.calcN!(ts.RHS, s.sol, s.t, s, v, p, g)
  @. ts.RHS += eq.LC*s.sol      # Add linear term to RHS

  @. s.sol += ts.dt*(23/12*ts.RHS - 16/12*ts.RHS₋₁ + 5/12*ts.RHS₋₂)
  s.t += ts.dt
  s.step += 1

  ts.RHS₋₂ .= ts.RHS₋₁          # Store
  ts.RHS₋₁ .= ts.RHS            # ... previous values of RHS

  nothing
end
