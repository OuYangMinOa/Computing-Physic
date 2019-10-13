using PyPlot
using Distributed
using Printf
pygui(true)

function integer_part(num)
    out = ""
    exp = 0
    while num > 1
        exp += 1
        num = num ÷ 2
        if (num%2==1)
            out = string('1',out)
        else
            out = string('0',out)
        end
    end
    return out
end

function float_part(num)
    two = 1
    out = ""
    exp = 0
    if (num >1)
        return "number must smaller than 1"
    end
    while num>10^(-5)
        exp += 1
        two /= 2
        if (num>=two)
            out = string(out,"1")
            num -= two
        else
            out = string(out,"0")
        end
    end
    return out 
end

function exercise_1()
    #while true
    #print("input your decimal number: ")
    n =  "7.156" #readline()
    integer_num = floor(parse(Float64,n))
    float_num = parse(Float64,n) - integer_num
    #print(integer_num,float_num)
    println("($(n))₂ = ($(integer_part(integer_num)).$(float_part(float_num)))₁₀")
    #end
end

#exercise_1()
################################################
function up(n)
    out::Float64 =  0.0
    for i::Float64 in 1:floor(n)
        out::Float64 += (1.0/i)
    end
    return out
end

function down(n)
    out::Float64 =  0.0
    for i::Float64 in floor(n):-1:1
        out::Float64 += 1.0/i
    end
    return out
end

function exercise_2()
    #subplots()
    pre  = 1
    nu = 10
    arr = []
    subplots()
    for i in 1:pre:nu
        u = up(10^i)
        d = down(10^i)
        p = abs(u-d)/d
        println("10^$i: $u, $d,誤差= $(p)")#, 誤差(log)= $(log10(abs(u-d)/d))")
        append!(arr,p)
        scatter(i,log10(p))
    end
    savefig("D:/計算物理/32log10 1~$(nu).png")
    subplots()
    plot(arr,"ro")
    plot(arr,"b--")
    savefig("D:/計算物理/32 linear 1~$(nu).png")
end


@time exercise_2()