#!/usr/bin/env julia

# Start Test Script
using FourierFlows
using Base.Test

# Run tests
tic()

@testset "Grids" begin
  include("test_grid.jl")
end

@testset "FFT" begin
  include("test_fft.jl")
end

@testset "IFFT" begin
  include("test_ifft.jl")
end

@testset "Timesteppers" begin
  include("test_timesteppers.jl")
end

@testset "Utils" begin
  include("test_utils.jl")
end

@testset "Physics: TwoDTurb           " begin
  include("test_twodturb.jl")    
end

@testset "Physics: Cosine Boussinesq  " begin
  include("test_verticallycosineboussinesq.jl")    
end

@testset "Physics: Fourier Boussinesq " begin
  include("test_verticallyfourierboussinesq.jl")    
end

#@time @testset "BarotropicQG and Timestepper tests" begin
#    include("test_BarotropicQG_timestep.jl")
#end

println("Total test time: ", toq())
