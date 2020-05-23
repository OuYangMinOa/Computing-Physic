using PyPlot
using FFTW
pygui(true)

global func(t) = exp(-t^2)*exp(2im*t)
global func_ft(w) = 1/√2*exp(-(2-w)^2/4)
A = 7

function DFT(array)
    N = length(array) # size of input array
    Δ = 2A/(N-1) # step = range/N-1
    sq = (sqrt(2pi)/Δ)
    DFT_out = Array{ComplexF64,1}(undef,N)
    DFT_w = Array{Float64,1}(undef,N)
    for j in 0:N-1
        wₖ = ((j-N/2)*2π/N/Δ)
        this = 0
        for k in 0:N-1
            this += (array[k+1]*exp(-im*2*π*(j-N/2)*k/N))::ComplexF64
        end # end for
        DFT_out[j+1] = this/sq*exp(im*wₖ*A)
        DFT_w[j+1] = wₖ
    end # end for
    return DFT_out,DFT_w
end

function make_shift_freq(n,step)
    #f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
    if (n % 2 ==1)
    out = collect(Iterators.flatten([collect(0:(n-1/2)),collect(-(n-1)/2:-1)])) / (step * n)
    elseif (n % 2 ==0)
    #f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
    out = collect(Iterators.flatten([collect(0:n/2), collect(-n/2:-1)]))/(step * n)
    end
    #shift freq time
    return fftshift(out)
end

function check_DFT(N)
    println("check_DFT time:")
    @time begin
    subplots()
    title("check DFT")
    arr = collect(range(-A,A,length=N))
    in_array = func.(arr)
    out_array,x_array = DFT(in_array)
    plot(x_array,out_array,"x",label ="DFT")
    plot(x_array,func_ft.(x_array),"+",label ="original")
    legend()
    grid()
    end
end

function check_FFT(N)
    println("check_FFT time:")
    @time begin
    subplots()
    title("check FFT")
    Δ = 2A/(N-1)
    # time 
    arr = collect(range(-A,A,length=N))
    # function before FFT
    in_array = func.(arr)
    # get x_array
    x_array = make_shift_freq(N,Δ).*2pi
    # absolute & fftshift & fft
    out = abs.(fftshift(fft(in_array)))
    # normalize
    fft_array = out*Δ/sqrt(2pi)
    # plot
    print("fft size:",size(fft_array),"x size:",size(x_array))
    plot(x_array[1:end-1],fft_array,"x",label ="FFT")
    plot(x_array,func_ft.(x_array),"+",label ="original")
    legend()
    grid()
    end
end

function show_func(N)
    arr = collect(range(-A,A,length=N))
    out_array,x_array = DFT(arr)
    subplot()
    title("-$A~$A N=$N")
    subplot(121)
    title("before ft")
    plot(arr,func.(arr),"-")
    grid()

    subplot(122)
    title("after ft")
    plot(x_array,func_ft.(x_array),"-")
    grid()
end
function main()
    N = 512
    #show_func(N)
    # check_DFT(N) # testing DFT
    # check_FFT(N) # testing FFT
    # time testing
    N_array = collect(5:15)
    DFT_t = Array{Float64,1}(undef,length(N_array))
    FFT_t = Array{Float64,1}(undef,length(N_array))
    for (i,N) in enumerate(N_array)
        this_array = func.(collect(range(-A,A,length=2^N)))
        DT = @elapsed DFT(this_array)
        FT =  @elapsed fftshift(fft(this_array))
        DFT_t[i] = DT
        FFT_t[i] = FT
    end # end for
    subplots()
    xlabel("log(N)")
    ylabel("time(s)")
    plot(N_array,DFT_t,label="DFT")
    plot(N_array,FFT_t,label="FFT")
    legend()
    plot(N_array,DFT_t,"o")
    plot(N_array,FFT_t,"o")
    println(" N     FFT      DFT")
    for (N,x,y) in zip(N_array,FFT_t,DFT_t)
        println("$N : $x   $y")
    end # end for
    grid()
end


# main()
# check_DFT(2048)
check_FFT(204800)