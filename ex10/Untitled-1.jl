using PyPlot
pygui(true)

function Fourth(N::Int64,a::Int64,b::Int64, initial= [0,1,0,1])
    y = [0.0]::Array{Float64}    # 儲存我要的y
    y_1 = [1.0]::Array{Float64}    # 儲存我要的y
    y_2 = [0.0]::Array{Float64}    # 儲存我要的y
    y_3 = [-1.0]::Array{Float64}    # 儲存我要的y
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
        append!(y_2,initial[3])
        append!(y_3,initial[4])
    end # end for
    return y,y_1,y_2,y_3,τ  # 返回 y 跟 間距
end # end function


function generation(t::Float64, last_point)
    fff = sin(t)::Float64
    return [last_point[2],last_point[3],last_point[4],fff]
end

function loss(y,y♡)
    return log10(sum((y.-y♡).^2)/length(y))
end

function main()
    y_ans(t) = sin(t)
    subplot()
    for i in 1:0.1:5
        down = 0
        N = trunc(Int64,2*(10^i))
        up = 5
        y,y1,y2,y3,h =  Fourth(N,down,up,[0,1,0,-1])
        arr = collect(down:h:up)
        losss = loss(y,y_ans.(arr))
        losss = log10(abs(y[end]-y_ans(5)))
        title("error")             # 標題
        xlabel("log10(step)")      # x座標
        ylabel("log10(error)")     # y座標
        scatter(log10(h),losss,label = "loss")
    end # end for
end 

main()