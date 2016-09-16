function doTimeSimp(x, window)
    N = 9
    M = floor(Int,length(x)/window)
    x = convert(Vector{Float64}, x)
    out = zeros(N,M)
    Threads.@threads for i=1:M
        loc      = x[(1:window)+window*(i-1)]
        out[1,i] = mean(loc)
        out[2,i] = maximum(loc)
        out[3,i] = minimum(loc)
        out[4,i] = std(loc)
        out[5,i] = detrend_fluctuation_analysis(loc)
        out[6,i] = fractal_dimension(loc)
        out[7,i] = shannon_entropy(loc)
        out[8,i] = approximate_entropy(loc)
        out[9,i] = sample_entropy(loc)
    end
    out
end

#This definition differs from the meta paper, but it should be consistent with
#the original source work
function fractal_dimension(x)
    M    = 100
    vals = [1:M...]
    out  = zeros(M)
    N    = length(x)
    for t=vals
        for i=1:t:(N-t)
            out[t] += sqrt((t/100)^2 + (x[i]-x[t+i])^2)
        end
    end
    #Add in the linear regression stage
    get_slope(log(out))
end

function get_slope(y::Vector{Float64})
    N = length(y)
    x = ones(N,2)
    o = zeros(N)
    for i=1:N
        x[i,2] = i
    end
    regression = x\y
    regression[2]
end

function linear_est(y::Vector{Float64})
    linear_est!(y, zeros(length(y)))
end

function linear_est!(y::Vector{Float64}, o::Vector{Float64})
    N::Int             = length(y)
    x::Matrix{Float64} = ones(N,2)
    @inbounds for i=1:N
        x[i,2] = i
    end
    regression::Vector{Float64} = x\y
    @inbounds for i=1:N
        o[i] = regression[1] + i*regression[2]
    end
end

function dfa_sub(y::Vector{Float64}, L::Int)
    N::Int                = length(y)
    yy::Vector{Float64}   = zeros(N)
    ytmp::Vector{Float64} = zeros(L)
    tmp::Vector{Float64}  = zeros(L)
    for i=1:L:N
        for j=1:L
            if(i+j-1 < N)
                ytmp[j] = y[i+j-1]
            end
        end
        linear_est!(ytmp, tmp)
        for j=1:L
            if(i+j-1 < N)
                yy[i+j-1] = tmp[j]
            end
        end
    end
    yy
end

function detrend_fluctuation_analysis(x::Vector{Float64})
    y::Vector{Float64} = cumsum(x.-mean(x))
    N::Int             = length(y)
    Lrange             = 10:5:100
    M::Int             = length(Lrange)
    pts                = zeros(M)
    for L=1:M
        yL     = dfa_sub(y, Lrange[L])
        pts[L] = sqrt(sum((y.-yL).^2)/N)
    end
    get_slope(pts)
end

function shannon_entropy(x)
    -sum((x.^2).*log(1e-8+x.^2))
end

function appx_d(x::Vector{Float64},m::Int,i::Int,j::Int)
    d::Float64 = 0
    for k=1:m
        d = max(d,abs(x[i+k-1]-x[j+k-1]))
    end
    d
end

function appx_c(x::Vector{Float64},m::Int,r::Float64,i::Int)
    c::Float64 = 0
    N::Int = length(x)-m
    for j=1:N
        c += appx_d(x,m,i,j)>r
    end
    c/N
end

function appx_phi(x::Vector{Float64},m::Int,r::Float64)
    #calculate c at each position
    o::Float64 = 0
    N::Int = length(x)-m 
    for i=1:N
        o += log(appx_c(x,m,r,i))
    end
    o/N
end

function samp_c(x::Vector{Float64},m::Int,r::Float64,i::Int)
    c::Float64 = 0
    N::Int = length(x)-m
    for j=1:N
        c += appx_d(x,m,i,j)<r
    end
    c
end

function samp_phi(x::Vector{Float64},m::Int,r::Float64)
    #calculate c at each position
    o::Float64 = 0
    N::Int = length(x)-m 
    for i=1:N
        o += samp_c(x,m,r,i)
    end
    o
end

#how in the world are they calculating the m/r/N parameters??
#Let's just arbitrarily say m = 4, r = 0.2
function approximate_entropy(x::Vector{Float64})
    r = 0.2
    appx_phi(x,3,r) - appx_phi(x,4,r)
end

#this depends upon the appx entropy func
function sample_entropy(x::Vector{Float64})
    r = 0.2*std(x)
    -log(samp_phi(x,3,r)/samp_phi(x,2,r))
end

function doTimeFeats()
    srange = [1:27...]
    for i=11:27
        print("processing subject $i @ ")
        println(now())
        h   = h5open("general-sleep/patient$i.h5")
        dat = read(h["2"])
        oo  = doTimeSimp(dat, 30*200)
        writecsv("feat-time-patient$i.csv", oo')
    end
end

function cheap_upsample(x)
    N = length(x)*2
    y = zeros(Float64, N)
    for i=1:N-1
        if(i%2 == 1)
            y[i] = x[Int((i+1)/2)]
        else
            ii = Int(i/2)
            y[i] = 0.5*(x[ii]+x[ii+1])
        end
    end
    y
end

froot = "/home/mark/current/general-sleep/physionet/non-edf/"
files = readdir(froot)
re    = r"\.h5$"
for ind=1:length(files)
    i = files[ind]
    if(match(re,i) != nothing)
        h = h5open(string(froot,i))
        dat = cheap_upsample(read(h["1"]))

        println("Processing ", i)
        println("   Starting Process @ ", now())
        println("   Hours To Process = ", length(dat)/200/60/60, " hours")
        if(isfile("data/feat-time-physionet-$i.csv"))
            println("       Skipping Already Processed File...")
        else
            println("       Working on file ", ind, "/", length(files))
            oo = doTimeSimp(dat, 30*200)
            writecsv("data/feat-time-physionet-$i.csv", oo)
        end
        close(h)
    end
end

#Keep in mind that physionet is 100Hz data unlike DREAMS's 200Hz
#just cheaply resample?

#obtain channel 1 for fpz-cz
#process
#store
#later grab processed data & labels via the convert-label functionality
# (more file globs needed)
#use forests to train/test
#find cross performance with other two data sets
#A-A
#B-B
#C-C
#A-B
#B-A
#A-C
#C-A
#B-C
#C-B
#TODO find some typical numbers for DREAMS <-> Physionet
