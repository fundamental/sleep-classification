function getSpectra(SubjectID, workDir, ext="")
    Sp = raw_spec(SubjectID)
    Cl = readcsv("$workDir/Dejunked$SubjectID$ext-cols.csv")
    Sp[:,find(Cl)]
    #imread("$workDir/DejunkedSpectra$SubjectID.png")
end

function obtainGradient(samples::Vector{Float64})
    #Optimize for the form of y=a*x+b
    N = length(samples)
    if(N == 1)
        return samples
    end
    (a,b) = [linspace(0,1,N) ones(N)]\samples
    out = zeros(N)
    for i=1:N
        out[i] = a*(i-1)/(N-1)+b
    end
    out
end

function zoneifyGradient(img, X, Y)
    I = float(copy(img))
    xedge = find(X)
    #println("xedge = ", xedge)
    yedge = find(Y)
    prepend!(xedge, [1])
    prepend!(yedge, [1])
    push!(yedge, size(img)[1])
    push!(xedge, size(img)[2])
    #println(size(img))
    #println(xedge)
    #println(yedge)
    #println(maximum(yedge))
    for yy=1:length(yedge)-1, xx=1:length(xedge)-1
        I[yedge[yy]:yedge[yy+1],xedge[xx]:xedge[xx+1]] = obtainGradient(img[yedge[yy]:yedge[yy+1],xedge[xx]:xedge[xx+1]][:])
    end
    I
end

function viewStuff(SubjectID::String, workDir::String, amount::Float64,
    misc=false,alt_work=nothing)
    ep  = readcsv("$workDir/edgeParameter$SubjectID.csv")
    dat = nothing
    if(alt_work != nothing)
        sub = alt_work[2]
        ext = alt_work[1]
        DD = readcsv("physionet-noise/$sub-noise$ext-dB.csv")[:]
        figure(999)
        (Spectra,_) = specgram(DD, 4096, 100, noverlap=0)

        PyPlot.close()
        Spectra = log(abs(Spectra.^2))

        Sp = Spectra
        Cl = readcsv("$workDir/Dejunked$SubjectID-cols.csv")
        Sp[:,find(Cl)]
        dat = Sp[:,find(Cl)]
    else
        dat = getSpectra(SubjectID, workDir)[1:1000,:]
    end
    ll  = readcsv("$workDir/Dejunkedlabels$SubjectID.csv")

    X = find(maximum(ep,2)[:])
    Y = find(ones(size(dat)[1]))
    prepend!(X,[1])
    push!(X,size(dat,2))

    if(misc)
        mn = zoneifyGradient(dat, X, Y)
        md = zoneify(dat, X, Y, operator=median)
    else
        mn = zoneify(dat, X, Y, operator=mean)
        md = zoneify(dat, X, Y, operator=median)
    end
    out = (1.0-amount)*dat+(amount)*(mn+md)/2.0
    (ll, out)
end
