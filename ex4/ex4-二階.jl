using PyPlot
using Distributed
using SymPy
pygui(true)
println("import ok")
-
function Fourth(N::Int64,a::Int64,b::Int64, initial= [0,0])
    y = [0.0]::Array{Float64}    # 儲存我要的y
    y_1 = [0.0]::Array{Float64}    # 儲存我要的y
    τ = (abs(b-a)/(N))::Float64  # 間距 = (上限-下限)/N 
    for t in a:τ:b-τ             # 迴圈
        ################## k₁ k₂ k₃ k₄ ##################
        c₁ = τ .* generation(t  ,initial)::Array{Float64}
        c₂ = τ .* generation(t+τ/2 ,initial + c₁./2)::Array{Float64}
        c₃ = τ .* generation(t+τ/2,initial + c₂./2)::Array{Float64}
        c₄ = τ .* generation(t+τ,initial + c₃)::Array{Float64}
        ################## yₙ = yₙ₊₁ + slope*h ##################
        initial += ((c₁+2c₂+2c₃+c₄)/6)::Array{Float64}  # 前一個 加上 斜率*dx
        append!(y,initial[1])    # 把 y 存起來
        append!(y_1,initial[2])
    end # end for
    return y,τ,y_1  # 返回 y 跟 間距
end # end function

function generation(t::Float64, last_point)
    fff = (exp(2t)*sin(t) + 2*last_point[2] - 2*last_point[1])::Float64
    #y''=  exp(2t)*sin(t) +      2y'        - 2y
    return [last_point[2],fff]
end

function exact(x)
    return 0.1*exp(x)*(4*cos(x) - 5*exp(x)*cos(x) + exp(x)*cos(x)*cos(2x)+2*sin(x)-2*exp(x)*cos(2x)*sin(x)+2*exp(x)*cos(x)*sin(2x)+exp(x)*sin(x)*sin(2x))
end

function draw_excat(up = 2)
    subplot()
    title("exact function 0~2")
    arr = collect(0:0.02:up)
    plot(arr,exact.(arr),label="exact function")
end

function main()
    down = 0::Int64
    up = 4::Int64
    ans = exact(2)::Float64
    println("real ans = $ans")

    N = 1000
    y,h,gg =  Fourth(N,down,up,[0,0])
    plot(collect(down:h:up),y,"-",label="y'")
    plot(collect(down:h:up),gg,"--",label="y")
    legend()
    #plot(collect(down:h:up),get,"*--",label="N = $N")
    #################### for ########################
    # for i in 0:0.1:5               # for 迴圈
    #     N = trunc(Int64,2*(10^i))  # 轉整數
    #     get,h,get_2 =  Fourth(N,down,up) # 去做runge kutta
    #     y = get[end]               # 最後一個點是答案
    #     error = log10(abs(ans-y))  # 誤差
    # ###################### 畫圖 ######################
    #     subplot()                  # 新圖
    #     title("error")             # 標題
    #     xlabel("log10(step)")      # x座標
    #     ylabel("log10(error)")     # y座標
    #     println(i,". step: ",h," ans: ",y," error: ",abs(ans-y)) # print
    #     scatter(log10(h),error,label="N = $N") # 畫點
    # end # end for
end # end function
#draw_excat()
# @time main()
# @time [Fourth(10000,0,2) for i in 1:1000]
# Fourth(100,0,2)

function real_ans()
    @syms t
    y = SymFunction("y")
    functiob = y''(t) + 0.9*cos(t*2/3) - 0.5*y'(t) + sin(y(t)) 
    get = dsolve(functiob,y(t), ics=((y, 0, 0),(y', 0, 2)))
    #print(N(simplify(get.rhs)))
end
real_ans()