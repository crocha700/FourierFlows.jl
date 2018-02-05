import Base: getindex
export AbstractDiagnostic, Diagnostic, update!, increment!

abstract type AbstractDiagnostic end

""" A diagnostic type associated with FourierFlows.Problem types """
mutable struct Diagnostic{Tdt<:AbstractFloat,T} <: AbstractDiagnostic
  calc::Function
  prob::Problem
  num::Int
  data::Array{T,1}
  time::Array{Tdt,1}
  step::Array{Int,1}
  value::T
  count::Int
  freq::Int
end

""" Constructor for the ProblemDiagnostic type. """
function Diagnostic(calc::Function, prob::FourierFlows.Problem; freq=1,
  nsteps=1, num=ceil(Int, (nsteps+1)/freq))

  value = calc(prob)
  T = typeof(value)

  data = Array{T}(num)
  time = zeros(typeof(prob.ts.dt), num)
  step = zeros(Int, num)

  data[1] = value
  time[1] = prob.t
  step[1] = prob.step

  Diagnostic(calc, prob, num, data, time, step, value, 1, freq)
end

getindex(d::Diagnostic, inds...) = getindex(d.data, inds...)
function getindex(d::Diagnostic, ::Colon)
  d.count < d.num ? d.data[1:d.count] : (
    cat(1, d.data[1:mod1(d.count, d.num)], d.data[mod1(d.count, d.num)+1:end]))   
end

function getindex(d::Diagnostic, time::Symbol)
  if d.count < d.num
    return getfield(d, :time)[1:d.count]
  else
    i = mod1(d.count, d.num)
    diagfld = getfield(d, :time)
    return cat(1, diagfld[1:i], diagfld[i+1:end])
  end
end

""" 
    update!(diag)

Update diag with its current value.
"""
function update!(diag::AbstractDiagnostic)
  i = mod1(diag.count, diag.num)
  diag.value = diag.data[i] = diag.calc(diag.prob)
  diag.time[i] = diag.prob.t
  diag.step[i] = diag.prob.step
  nothing
end

update!(diags::AbstractArray) = for diag in diags; update!(diag); end
  
""" 
    increment!(diag)

Increment the Diagnostic diag.
"""
function increment!(diag::AbstractDiagnostic)
  diag.count += 1
  i = mod1(diag.count, diag.num)
  diag.value = diag.data[i] = diag.calc(diag.prob)
  diag.time[i] = diag.prob.t
  diag.step[i] = diag.prob.step
  nothing
end

increment!(diags::AbstractArray) = for d in diags; increment!(d); end
