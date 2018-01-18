export Equation, Problem, State, DualState

# Equation Composite Type
"""
This type defines the linear implicit and explicit components of an equation.
The linear implicit part of an equation is defined by an array of coefficients
which multiply the solution. The explicit part of an equation is calculated
by a function that may define linear and nonlinear parts.
"""
struct Equation{T,dim} <: AbstractEquation
  LC::Array{T,dim} # Coeffs of the eqn's implicit linear part
  calcN!::Function # Function that calcs linear & nonlinear parts
end

struct DualEquation{Tc,Tr,dimc,dimr} <: AbstractEquation
  LCc::Array{Tc,dimc}
  LCr::Array{Tr,dimr}
  calcN!::Function
end


# Problem state 
mutable struct State{Tdt,Tsol,dim} <: AbstractState
  t::Tdt
  step::Int
  dt::Tdt
  sol::Array{Tsol,dim}
end

function State(Tsol::DataType, sz::Tuple, dt)
  sol = zeros(Tsol, sz)
  State(0.0, 0, dt, sol)
end

mutable struct DualState{Tdt,Tc,Tr,dimc,dimr} <: AbstractState
  t::Tdt
  step::Int
  dt::Tdt
  solc::Array{Tc,dimc}
  solr::Array{Tr,dimr}
end

function DualState(Tc::DataType, sizec, Tr::DataType, sizer, dt)
  solc = zeros(Tc, sizec)
  solr = zeros(Tr, sizer)
  DualState(0.0, 0, dt, solc, solr)
end

# Problem type and associated functions
mutable struct Problem <: AbstractProblem
  grid::AbstractGrid
  vars::AbstractVars
  params::AbstractParams
  eqn::AbstractEquation
  ts::AbstractTimeStepper
  state::AbstractState
  t::Float64
  step::Int
end

#=
For v1.0 release:
  The `t` and `step` properties can be removed from Problem, instead
  overloading the `getfield` function to intercept attempts to access
  `t` and `step`, which will be redirected to prob.state.t and 
  prob.state.step (and perhaps sol as well). This'll make lots of things
  just a little bit nicer.

import Base: getfield

function getfield(prob::AbstractProblem, name)
  if name ∈ fieldnames(prob.state)
    return getfield(prob.state, name)
  else
    return getfield(prob, name)
  end
end
=#

Problem(g, v, p, eq, ts, st) = Problem(g, v, p, eq, ts, st, st.t, st.step)

function Problem(g, v, p, eq::Equation, ts; soltype=cxeltype(eq.LC))
  st = State(0.0, 0, ts.dt, zeros(soltype, size(eq.LC)))
  Problem(g, v, p, eq, ts, st)
end

function Problem(g, v, p, eq::DualEquation, ts; 
                 solctype=cxeltype(eq.LCc), solrtype=cxeltype(eq.LCr))
  st = DualState(solctype, size(eq.LCc), solrtype, size(eq.LCr), ts.dt)
  Problem(g, v, p, eq, ts, st)
end
