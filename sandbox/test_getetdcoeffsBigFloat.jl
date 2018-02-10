using FourierFlows

n, L = 128, 2π
g = TwoDGrid(n, L)
ν = 1e-3
LC = -ν*g.KKrsq.^2  + im*g.Kr
LC[1, 1]=0
dt = 0.1




function getetdcoeffs(dt, LC; ncirc=32, rcirc=eltype(LC)(14))

  shape = Tuple(cat(1, ncirc, ones(Int, ndims(LC))))

  circ = zeros(FourierFlows.cxeltype(LC), shape)
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

function getetdcoeffsBigFloat(dt, LC)
  L = dt*LC
  M = ndims(LC)+1

  expL  = big.(exp.(L))
  expL2 = big.(exp.(L/2))
  L = big.(L)

  ζc = @. Complex{Float64}( ( expL2-1 ) / L )
  αc = @. Complex{Float64}( (-4 - L + expL*(4-3L + L^2))/L^3 )
  βc = @. Complex{Float64}( ( 2 + L + expL*(-2 + L))/L^3 )
  Γc = @. Complex{Float64}( (-4 - 3L -L^2 + expL*(4 - L))/L^3 )


  if eltype(LC) <: Real
    ζ = dt*real.(ζc)
    α = dt*real.(αc)
    β = dt*real.(βc)
    Γ = dt*real.(Γc)
  else
    ζ = dt*ζc
    α = dt*αc
    β = dt*βc
    Γ = dt*Γc
  end

  ζ, α, β, Γ
end



@time ζ1, α1, β1, Γ1 = getetdcoeffs(dt, LC, ncirc=32)
ζ1[1, 1], α1[1, 1], β1[1, 1], Γ1[1, 1] = 0, 0, 0, 0
@time ζ2, α2, β2, Γ2 = getetdcoeffsBigFloat(dt, LC)
ζ2[1, 1], α2[1, 1], β2[1, 1], Γ2[1, 1] = 0, 0, 0, 0


println("diff error ζ: ", norm(ζ1-ζ2, Inf))
println("diff error α: ", norm(α1-α2, Inf))
println("diff error β: ", norm(β1-β2, Inf))
println("diff error Γ: ", norm(Γ1-Γ2, Inf))
