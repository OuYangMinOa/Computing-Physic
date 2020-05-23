using PyPlot
pygui(true)
global const ħ = 1
global const m = 1 
global v(x) = 0.5x^2
global q(ϵ,x) = 2m/(ħ^2) * (ϵ-v(x))
global const initial = 0
global s(x) = 0
function Numerov(ϵ,up,down,N,initial,pp=false,factor=1,next = 10.0^-10)
    u = [initial,next]
    all_length = N-2
    h = abs(up-down)/(N)
    for (num,i) in enumerate(range(down+2h,step=h,length=all_length))
        c_plus_i = (1 + h^2/12*q(ϵ,i+h)) 
        c_minus_i = 1 + h^2/12*q(ϵ,i-h)
        c_i = 2 -5h^2/6*q(ϵ,i)
        di = h^2/12*(s(i+h) + 10*s(i) + s(i-h))
        next_u = (c_i*u[num+1]+di-c_minus_i*u[num])/c_plus_i
        append!(u,next_u)
    end #for
    u_r = [initial,next*factor]
    for (num,i) in enumerate(range(-down-2h,step=-h,length=all_length))
        c_plus_i = 1 + h^2/12*q(ϵ,i+h)
        c_minus_i = 1 + h^2/12*q(ϵ,i-h)
        c_i = 2 -5h^2/6*q(ϵ,i)
        di = h^2/12*(s(i+h) + 10*s(i) + s(i-h))
        next_u = (c_i*u_r[num+1]+di-c_minus_i*u_r[num])/c_plus_i
        append!(u_r,next_u)
    end #for
    if (pp)
        #scatter(range(down+2h,step=h,length=all_length+2),u)
        plot(range(down+2h,step=h,length=all_length+2),u)
        plot(range(-down-2h,step=-h,length=all_length+2),u_r)
        #scatter(range(-down-2h,step=-h,length=all_length+2),u_r)
    end # end if
    return u,u_r,h
end

function positiveenergy_answer(ϵ,pp=false)
    left, right, h = Numerov(ϵ,0,-5,100000,initial,pp,1)
    return ((left[end] - left[end-2])/(2*h*left[end-1]) - ( right[end]-right[end-2] )/(2*h*right[end-1]))
end

function ne_energy_answer(ϵ,pp=false)
    left, right, h = Numerov(ϵ,0,-5,100000,initial,pp,-1)
    #println( ((left[end] - left[end-2])/(2*left[end-1]) ,"\n",  (right[end-2] - right[end])/(2*right[end-1])))
    return ((left[end] - left[end-2])/(2*left[end-1]) -  (right[end-2] - right[end])/(2*right[end-1]))
end

function newtowns(x,func = positiveenergy_answer)
    if (abs(func(x))< 10.0^-6)
        return x
    end # end if
    h = 10^-10
    f2 = (func(x+h)-func(x-h))/(2h)
    newtowns(x-func(x)/f2,func)
end # end function

function main()
    for i in 1:9
    if (i%2==1)
    ans = newtowns(i-0.5)
    println("answer = ",ans)
    subplots()
    title("ϵ = $ans")
    positiveenergy_answer(ans,true)
    else
    ANS_2 = newtowns(i-0.5,ne_energy_answer)
    println("answer = ",ANS_2)
    subplots()
    title("ϵ = $ANS_2")
    ne_energy_answer(ANS_2,true)
    end
    end # end for
end # end function
subplots()
arr = collect(0:0.1:5)
title("ϵ v.s continuity")
scatter(arr,positiveenergy_answer.(arr)+ne_energy_answer.(arr))

# main()
