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
    (c,a,b) = [linspace(0,1,N).^2 linspace(0,1,N) ones(N)]\samples
    out = zeros(N)
    for i=1:N
        out[i] = c*((i-1)/(N-1))^2 + a*(i-1)/(N-1)+b
    end
    out
end

function obtainGradient2D(samples::Matrix{Float64})
    svec = zeros(length(samples))
    #row quad term
    #row linear term
    #column quad term
    #row linear term
    #constant bias
    lp = zeros(length(svec),5)
    M = size(samples,1)
    N = size(samples,2)
    for i=1:M,j=1:N
        ind = i+(j-1)*M
        svec[ind] = samples[i,j]
        lp[ind,1] = ((j-1)/N)^2
        lp[ind,2] = ((j-1)/N)
        lp[ind,3] = ((i-1)/N)^2
        lp[ind,4] = ((i-1)/N)
        lp[ind,5] = 1
    end
    (cr,ar,cc,ac,b) = lp\svec

    out = zeros(M,N)
    for i=1:M, j=1:N
        out[i,j] = cr*((j-1)/N)^2 + ar*((j-1)/N) + cc*((i-1)/N)^2 + ac*((i-1)/N) + b
    end
    out
end

function zoneifyGradient(img, X, Y; operator=mean)
    I = float(copy(img))
    I[:,:] = 0
    xedge = X;find(X)
    #println("xedge = ", xedge)
    yedge = Y;find(Y)
    #prepend!(xedge, [1])
    #prepend!(yedge, [1])
    #push!(yedge, size(img)[1])
    #push!(xedge, size(img)[2])
    #println(size(img))
    #println(xedge)
    #println(yedge)
    #println(maximum(yedge))

    #XXX this old code only worked for when Y was densely sampled
    #for yy=1:length(yedge)-1, xx=1:length(xedge)-1
    #    I[yedge[yy]:yedge[yy+1],xedge[xx]:xedge[xx+1]] = obtainGradient(img[yedge[yy]:yedge[yy+1],xedge[xx]:xedge[xx+1]][:])
    #end
    for yy=1:length(yedge)-1, xx=1:length(xedge)-1
        y1 = yedge[yy]
        y2 = yedge[yy+1]
        x1 = xedge[xx]
        x2 = xedge[xx+1]
        #println([y1, y2, x1, x2])
        #println((y2-y1)*(x2-x1))
        row_mean = operator(img[y1:y2,x1:x2],2)
        row_grad = obtainGradient(row_mean[:])
        col_mean = operator(img[y1:y2,x1:x2],1)
        col_grad = obtainGradient(col_mean[:])
        eh = obtainGradient2D(img[y1:y2,x1:x2])
        for y=y1:y2, x=x1:x2
            #I[y,x] = (row_grad[1+y-y1]+col_grad[1+x-x1])/2
            I[y,x] = eh[1+y-y1,1+x-x1]
        end
    end
    I
end

function viewStuff(SubjectID::String, workDir::String, amount::Float64,
    misc=false,alt_work=nothing;x_thresh=0,y_thresh=-1)
    ep  = readcsv("$workDir/edgeParameter$SubjectID.csv")
    dat = nothing
    ll  = nothing
    if(alt_work != nothing && length(alt_work) == 2)
        sub = alt_work[2]
        ext = alt_work[1]
        DD = readcsv("physionet-thresh/$sub-noise$ext-dB.csv")[:]
        figure(999)
        (Spectra,_) = specgram(DD, 4096, 100, noverlap=0)

        PyPlot.close()
        Spectra = log(abs(Spectra.^2))

        Sp = Spectra
        Cl = readcsv("$workDir/Dejunked$SubjectID-cols.csv")
        Sp[:,find(Cl)]
        dat = Sp[:,find(Cl)]
        ll  = readcsv("$workDir/Dejunkedlabels$SubjectID.csv")
    elseif(alt_work != nothing && length(alt_work) == 1)
        dat = getSpectra("ST7192", "physionet")[1:1000,:]
        ll  = readcsv("physionet/DejunkedlabelsST7192.csv")
    else
        dat = getSpectra(SubjectID, workDir)[1:1000,:]
        ll  = readcsv("$workDir/Dejunkedlabels$SubjectID.csv")
    end
    #ll  = readcsv("$workDir/Dejunkedlabels$SubjectID.csv")

    X = nothing
    Y = nothing
    X = find(maximum(ep,2)[:].>x_thresh)
    if(y_thresh == -1)
        Y = find(ones(size(dat)[1]))
    else
        coverhigh = PyPlot.imread("$workDir/coverhigh$SubjectID.png")
        bandshigh = identifyInterestBands(map(Float64, coverhigh[1:1000,:]))
        covermed = PyPlot.imread("$workDir/covermed$SubjectID.png")
        bandsmed = identifyInterestBands(map(Float64, covermed[1:1000,:]))
        coverlow = PyPlot.imread("$workDir/coverlow$SubjectID.png")
        bandslow = identifyInterestBands(map(Float64, coverlow[1:1000,:]))
        Y     = find((bandshigh[1:end-1] .!= bandshigh[2:end]) |
                     (bandsmed[1:end-1]  .!= bandsmed[2:end])  |
                     (bandslow[1:end-1]  .!= bandslow[2:end]))
        prepend!(Y, [1])
        push!(Y, 1000)
    end
    println("X.length = ", length(X)*100/size(dat,2),"%")
    println("Y.length = ", length(Y)*100/size(dat,1),"%")
    prepend!(X,[1])
    push!(X,size(dat,2))

    if(misc)
        mn = zoneifyGradient(dat, X, Y, operator=mean)
        md = zoneifyGradient(dat, X, Y, operator=median)
    else
        mn = zoneify(dat, X, Y, operator=mean)
        md = zoneify(dat, X, Y, operator=median)
    end
    out = (1.0-amount)*dat+(amount)*(mn+md)/2.0
    (ll, out)
end
