using PyPlot
using FFTW
using Printf
pygui(true)

# -N/2 ... N/2-1
func(t) = exp(-t^2)*exp(2im*t)
func_ft(w) = 1/√2*exp(-(2-w)^2/4)

function DFT(array)
    out_array = []
    N = length(array)
        for j in 0:N-1
            this = 0#::ComplexF64
            for k in 0:N-1
                this += (array[k+1]*exp(-im*2*π*(j-N/2)*k/N))::ComplexF64
            end # end for
            append!(out_array,1/(sqrt(N))*this)
            #println(1/(sqrt(N))*this)
        end # end for
    return out_array
end


function fourier_draw(N,A)
    subplots()
    Δ = 2A/(N-1)
    fc = 1/(2Δ)
    w_max = pi/Δ
    arr = collect(range(-A,A,length=N))
    in_array = func.(arr)
    arr_ft = collect(range(-w_max,w_max,length=N+1))
    out = DFT(in_array)
    new ,ww = [],[] 
    for (num,i) in enumerate(out)
        wₖ = ((num-1-N/2)*2π/N/Δ)
        append!(ww,wₖ)
        append!(new,i*exp(im*wₖ*A))
    end
    new = new./ ( sqrt(2pi)/sqrt(N)/(2A/(N-1)) )
    plot(ww,real.(new),"x")
    #plot(ww,real.(new),"--",label="N=$N, range=( $(-A) ~ $A )")
    #plot(ww,func_ft.(ww),label="exact")
    plot(ww,func_ft.(ww),"+")
    grid(true)
    legend()
    return ww,real.(new),func_ft.(ww)
end
function main() 
    ##################
    A = 7::Int64
    N = 2048::Int64
    fourier_draw(N,A)
    ################
end


function write_data(N,out)
    file_name1 = "D:\\計算物理\\ex6\\K-DFT(RE)-DFt(Im).txt"
    file_name2 = "D:\\計算物理\\ex6\\w-DFT(RE)-exact(RE).txt"
    open(file_name1,"w") do io
        write(io,"   k                            DFT(re)                         DFT(im)\n")
        for (num,DFT) in enumerate(out)
            k = num-1-N/2
            write(io,"$(@sprintf("%10.5f     %25.20f     %20.25f\n",k,real(DFT),imag(DFT)))")
        end # end for 
    end # end open
    open(file_name2,"w") do io
        write(io,"   w                            DFT(re)                         exact(re)\n")
        for (num,DFT) in enumerate(out)
            wₖ = ((num-1-N/2)*2π/N/Δ)
            write(io,"$(@sprintf("%10.5f     %25.20f     %20.25f\n",wₖ,real(DFT),real(func_ft(wₖ))))")
        end # end for 
    end # end open
end # end function

@time main()