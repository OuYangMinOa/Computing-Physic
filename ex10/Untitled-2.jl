using PyPlot
pygui(true)

y(t) = sin(t)
y1(t) = cos(t)
y2(t) = -sin(t)
y3(t) = -cos(t)

down = 0
N = 1000
up = 5

x_asix = collect(range(down,up,length=N))
plot(x_asix,y.(x_asix),label = "y")
plot(x_asix,y1.(x_asix),label = "y'")
plot(x_asix,y2.(x_asix),label = "y''")
plot(x_asix,y3.(x_asix),label = "y'''")
legend()