using PyPlot
pygui(true)
global f(y,t) = exp(t)
global ans_func(t) =exp(t)
global f0 = 1.0
global down = 0
global up = 10

function Fourth(N,a,b,pp = false)
    y::Array{Float64} = [f0]
    step = abs(b-a)/N
    print("N = $N")
    x = collect(a:step:b)
    for i in 1:N
        # println(x[i],y[i])
        c1 = step*f(y[i],x[i])
        c2 = step*f(y[i]+c1/2,x[i]+step/2)
        c3 = step*f(y[i]+c2/2,x[i]+step/2)
        c4 = step*f(y[i]+c3,x[i]+step)
        append!(y,y[i] + (c1+(2*c2)+(2*c3)+c4)/6)
        if (pp)
            println("y $(y[i])")
        end
    end # end for
    return x,y
end

function main1()
    
    N = 500
    x,y = Fourth(N,down,up)
    println(size(x),"\n",size(y))
    # println(x)
    title("Fourth order Runge-Kutta N=$N")
    plot(x,y)
    subplots()
    title("Real answer N=$N")
    real = map(ans_func,x)
    plot(x,real)
    subplots()
    title("log(error) N=$N")

    error = []
    for i in 1:length(x)
        this = abs(y[i] - real[i])/real[i]
        append!(error,this)
    end # end for 
    plot(log.(x),log.(error))
    println(error[end])
end

function main2()
    ans = ans_func(up)
    error_X = []
    error2 = []
    for i in 0:7
        N = 10^i
        x,y = Fourth(N,down,up)
        append!(error_X,log10(x[2]-x[1]))
        append!(error2,abs(y[end]-ans)/ans)
        println(" , error = $(abs(y[end]-ans)) , log10(error) = $(log10(abs(y[end]-ans)/ans))")
    end # end for 
    plot(error_X,log10.(error2))
end
@time main2()

# @time Fourth(1000,down,up,true)
∙
