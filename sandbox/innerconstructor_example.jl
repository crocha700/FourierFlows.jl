abstract type AbstractTest end

mutable struct teststruct{dim} <: AbstractTest
    dt::Float64
    N::Array{Complex{Float64},dim}
    teststruct(dt::Float64, N::Array{Complex{Float64},dim}) where dim = new{dim}(dt, zeros(size(N)))
end
#teststruct(dt::Float64, N) = teststruct{ndims(N)}(dt, N)


dt=0.1;a2=im*ones(2, 2); a3 = im*ones(3, 3, 2);

a = teststruct(dt, a3);

a.N === a3
