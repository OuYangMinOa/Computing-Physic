using LinearAlgebra
using PyPlot
using SymPy
using Distributed
pygui(true)
######### function #############
global f(x) = exp(x) + exp(x)*sin(x) # formula
################################
function Definite_integral()
    @vars x
    println(integrate(exp(x) + exp(x)*sin(x),x))
    println()
    return N(
        integrate(exp(x) + exp(x)*sin(x),(x,0,5))
        )  # real answer 
end

function pp()
    arr = collect(0:0.001:5)
    plot(arr,f.(arr))
end
function simpson(n::Int64,up::Int64,down::Int64,Func = f)
    step = (up - down) / n
    ans = 0
    for i in 0:n-1
        a = down + step*i
        b = down + step*(i+1)
        ans += (b-a)/6*(Func(a)+Func(b)+4*Func((a+b)/2))
    end # for end 
    return ans
end # end function

function read_file()
    real_ans::Float64 = Definite_integral()
    out_ans = []
    out_len = []
    cut_to_time = 7
    file_name = "D:\\計算物理\\ex3\\Legendre(1).dat"
    open(file_name,"r") do io
        while true
        get = readline(io)
        if (get == "")
            break
        else
            get = parse(Int64,get)
            node = zeros(get)
            weight = zeros(get)
            
            for i in 1:get
                line = lstrip(readline(io))
                if (get >= 64)
                    line = replace(line,"D"=>"E")
                end # end if
                this_line = []
                for num in split((line)," ")
                    try 
                        new = parse(Float64,num)
                        append!(this_line,new)
                    catch err
                        #println(err)
                    end
                end # end for
                #println(this_line)
                node[i] = this_line[1]
                weight[i] = this_line[2]
            end # end for
            println(repeat("=",25))
            println(get)
            println(repeat("=",25))
            println("cut         answer         error(%)")
            for cut in 0:cut_to_time
                hhh = 10^cut   # real_ans
                gettt = Gaussian_quadrature(get, node, weight,hhh)
                # append!(out_len,get)
                # append!(out_ans,Gaussian_quadrature(get, node, weight,hhh))
                println("$hhh   $(gettt)  error = $(abs(real_ans-gettt)/real_ans)%")
            end # end for
        end # end if
        end # end while
    end; # end open
    return out_len,out_ans

end # end function

function Gaussian_quadrature(n, node_array, weigth_array,cut=1)
    #println(n," ",length(node_array)," ",length(weigth_array))
    #print(node_array,weigth_array)
    a = 0
    b = 5
    step = (b - a)/cut
    ans = 0
    for j in 0:cut-1
        up_a = a + step*j
        up_b = a + step*(j+1)
        for i in 1:n
            ans += f(
                ((up_b-up_a)*node_array[i]+(up_b+up_a))/2
                )*step/2*weigth_array[i]
        end #end for
    end # end for 
    return ans
end

function main()
    println(repeat("=",60))
    real_ans::Float64 = Definite_integral()
    print("Definite integral answer = $real_ans")
    println("\n",repeat("=",60))
    println("Simpson's Integral:\n")
    for xx in 0:6
        cut_n = 10^xx
        simp_ans::Float64 = simpson(cut_n,5,0)
        print("$cut_n  $simp_ans  ")
        println("error = $(abs(real_ans-simp_ans)/real_ans*100)%")
    end # end for
    each_n,guass_ans = read_file()
    println("\n",repeat("=",60))
    println("Gaussian Quadrature: (歐陽出品)\n")
    # for i in 1:length(guass_ans)
    #     println("$(each_n[i])  $(guass_ans[i])  error = $(abs(real_ans-guass_ans[i])/guass_ans[i]*100)%")
    # end # end for
end

@time main()
# println()
# println( Definite_integral())











