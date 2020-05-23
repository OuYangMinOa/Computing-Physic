using PyPlot
using Distributed
pygui(true)

function DFT(array,Δ,ran)
    N = length(array)
    out_array = Array{ComplexF64,1}(undef,N)
    omaga = Array{Float64,1}(undef,N)
    for j in 0:N-1
        this = 0#::ComplexF64
        omaga[j+1] = (j)*2π/N/Δ
        for k in 0:N-1
            this += (array[k+1]*exp(-im*2*π*(j-N/2)*k/N))::ComplexF64
        end # end for          
        out_array[j+1] = this/(sqrt(N)/(sqrt(2pi)/sqrt(N)/Δ))
    end # end for
    return out_array,omaga
end


function Fourth(N::Int64,a::Int64,b::Int64, initial= [0,2])
    y = [0.0]::Array{Float64}    # my y
    τ = (abs(b-a)/(N))::Float64  # step
    for t in a:τ:b-τ             
        c₁ = τ .* generation(t  ,initial)::Array{Float64}
        c₂ = τ .* generation(t+τ/2 ,initial + c₁./2)::Array{Float64}
        c₃ = τ .* generation(t+τ/2,initial + c₂./2)::Array{Float64}
        c₄ = τ .* generation(t+τ,initial + c₃)::Array{Float64}
        initial += ((c₁+2c₂+2c₃+c₄)/6)::Array{Float64}
        append!(y,initial[1])
    end # end for
    return y,τ
end # end function

function generation(t::Float64, last_point)
    fff = (0.9*cos(t*2/3) - 0.5*last_point[2] - sin(last_point[1]))::Float64   
    return [last_point[2],fff]
end

function main()
    for i in 6:16
        N = 2^i
        println(N)
        a = 1::Int64
        b = 2^8::Int64
        Δ = (b-a)/(N-1)
        @time y,step = Fourth(N,a,b)
        out,ww = DFT(y,Δ,b-a)
        arr = collect(range(a,b,step=step))
        subplots()
        ion()
        title("N=$N range(0~$b) original")
        plot(arr,y)
        grid(true)
        subplots()
        ion()
        title("N=$N range(0~$b) log-scale")
        plot(ww,log10.(abs.(out).^2),"--")
        grid(true)
    end # end for 
end # end function

@time main()