using PyPlot
using ProgressMeter
using Clustering


"""
Identify Local Maxima
"""
function findPeakSeq(est::Vector{Float64})
    tmp = (est[1:end-2].<est[2:end-1]) & (est[2:end-1].>est[3:end])
    prepend!(tmp, [false])
    append!(tmp, [false])
    tmp
end

"""
Identify Local Maxima

In areas where a flat peak is found, this method returns the middle of the flat
region
"""
function findPeakSeq2(X::Vector{Float64})
    N::Int            = length(X)
    edge::Vector{Int} = zeros(Int,N)
    prev::Float64     = X[1]
    for i=2:N-1
        cur = X[i]
        if(cur > prev)
            edge[i] = +1
        elseif(cur < prev)
            edge[i] = -1
        end
        prev = cur
    end

    #Prune fall after fall
    state::Int = edge[1]
    for i=2:N-1
        if(state == -1 && edge[i] == -1)
            edge[i] = 0
        end
        if(edge[i] != 0)
            state = edge[i]
        end
    end

    edge = flipdim(edge,1)
    #Prune rise after rise
    state = edge[1]
    for i=2:N-1
        if(state == +1 && edge[i] == +1)
            edge[i] = 0
        end
        if(edge[i] != 0)
            state = edge[i]
        end
    end
    edge = flipud(edge)

    peaks::Bool    = zeros(Bool,N)
    state          = 0
    rise_time::Int = 0

    for i=2:N-1
        if((state == 0 || state == -1) && edge[i] == +1)
            state = +1
            rise_time = i
        elseif(state == +1 && edge[i] == -1)
            peaks[int(floor(Int,(rise_time+i)/2))] = true
            state = -1
        end
    end
    peaks
end

"""
Get Initial estimate of where various states begin and end based upon region
transitions
"""
function getSpectralLines(seq::Vector{Float64}, scale::Real, N::Int; baseFig::Int=-1)
    scale_ = 1
    est = zeros(scale_*N)

    N = length(seq)
    p = Progress(length(est)*N)
    for i=1:N, j=1:length(est)
        next!(p)
        evalPt = (j-1)/scale_
        est[j] += (scale[i])*exp(-(evalPt-seq[i])^2/8)
    end

    pe = findPeakSeq(est)

    (est,pe)
end

"""
Use K-means to Generate a Rough Estimate of low rank band structure
"""
function identifyInterestBands(img::Matrix{Float64})
    c  = cov(img')
    c -= mean(c)
    mc = mapslices(x->linearMedian(x,9),c,2)
    k  = kmeans(mc, 8, maxiter=20)
    diffsig = sign(k.centers[1:end-1,:]).!=sign(k.centers[2:end,:])
    roi = sum(diffsig,2)[:]
    roi = conv([1,1,1], roi)
    eliminateTransitionNoise(max(0,roi-1))
end

"""
Eliminate modality transitions which seem to occur too frequently
(this is a heuristic method to reduce the total variation with a prior that
state variations are expected to occur more frequently at lower frequencies)
"""
function eliminateTransitionNoise(bands::Vector{Float64}, Hz::Real=100, thresh::Real=0.3)
    N = length(bands)
    classOut = zeros(N)
    logHz = linspace(0, Hz, N)
    switchPoint = 0
    for i=4:N
        trySwap = bands[i] != 0
        atSwap  = switchPoint < logHz[i]
        if(atSwap && trySwap)
            classOut[i] = classOut[i-1]+1
        else
            classOut[i] = classOut[i-1]
        end

        if(trySwap)
            switchPoint = logHz[i] + thresh
        end
    end
    classOut
end


"""
Get estimate of the likelyhood of a transition between two temporal frames
"""
function getFreqEst(img::Matrix{Float64}, cover::Matrix{Float64},
    x::Matrix{Int}, gate::Union{Void,Vector{Float64}}=nothing)
    flush(STDOUT)
    println("getFreqEst...")
    doPlot = false
    sequence_x      = vcat(x[:,1][:],x[:,2][:])
    sequence_scalex = vcat((x[:,4]-x[:,3])[:], (x[:,4]-x[:,3])[:])
    sequence_y      = vcat(x[:,3][:],x[:,4][:])
    sequence_scaley = vcat((x[:,2]-x[:,1])[:], (x[:,2]-x[:,1])[:])

    estx::Vector{Float64}
    esty::Vector{Float64}
    px::Vector{Float64}
    py::Vector{Float64}
    (estx,px) = getSpectralLines(sequence_x, (sequence_scalex.*sequence_scaley).^0.5, size(img,2), baseFig=0)
    (esty,py) = getSpectralLines(sequence_y, (sequence_scaley.*sequence_scalex).^0.5, size(img,1), baseFig=10)
    if(gate != nothing)
        println("gate.length = ", length(gate));
        println("gate.sum = ",    sum(gate));
        println("px.length = ", length(px));
        px[2:end] = (1.0*gate)
        estx+=1e-6
    end

    pex = estx.*px
    pey = esty.*py
    py[pey.<sort(pey[py.>0])[round(Int,0.7end)]] = 0

    #Use the prior that variance over frequencies is coarse
    intervalY = find(py)[2:end]-find(py)[1:end-1]
    Fy::Vector{Int} = find(py)
    for i=1:length(intervalY)
        if(intervalY[i] < 20)
            py[Fy[i+1]] = 0
        end
    end


    if(doPlot)
        figure(13);plot(px.*estx)
        if(gate != nothing)
            plot(gate)
        end
    end
    X::Vector{Float64} = px.*estx
    Y::Vector{Float64} = py.*esty

    Y = identifyInterestBands(cover)
    Y = Y[1:end-1] .!= Y[2:end]
    Z = X'.*Y
    if(doPlot)
        figure(21);
        plt.clf();
        imshow(img,aspect="auto",interpolation="none");
        imshow(log(Z+0.1),aspect="auto",interpolation="none",alpha=0.5)
    end

    println("Extracting Features")
    feat = Array(Any,10)
    p = Progress(length(feat))

    #Create the list of edges
    xedge = find(X)
    yedge = find(Y)
    prepend!(xedge, [1])
    prepend!(yedge, [1])
    push!(yedge, size(img)[1])
    push!(xedge, size(img)[2])

    feat[1] = zoneify(img,xedge,yedge,x->sort(x[:])[round(Int,0.75end)])
    next!(p)
    feat[2] = zoneify(img,xedge,yedge,x->median(x[:]))
    next!(p)
    feat[3] = zoneify(img,xedge,yedge,x->sort(x[:])[round(Int,0.25end)])
    next!(p)
    feat[4] = zoneify(img,xedge,yedge,mean)
    next!(p)
    feat[5] = zoneify(img-mean(img),xedge,yedge,x->mean(sign(x[:])))
    feat[6] = float(mapslices(sortperm,  feat[1], 1))
    feat[7] = float(mapslices(sortperm,  feat[2], 1))
    feat[8] = float(mapslices(sortperm,  feat[3], 1))
    feat[9] = float(mapslices(sortperm,  feat[4], 1))
    feat[10] = float(mapslices(sortperm, feat[5], 1))


    diff_sig::Vector{Float64} = zeros(size(feat[1],2)-1)

    for j=1:20
        for i=1:length(feat)
            r::Vector{Int} = kClassify_nway(feat[i], 8, plotID=256, maxiter=20)
            diff_sig_   = r[1:end-1].!=r[2:end]
            diff_sig += diff_sig_#.*estx[1:end-1]#1.0/sum(diff_sig)
        end
    end

    diff_freq::Vector{Float64} = diff_sig.*estx[1:end-1]./sum(estx[1:end-1])

    (diff_freq, diff_sig, estx[1:end-1])
end

"""
Apply an operator to sequence of data using a fixed size odd window

Arguments:
- data   - The 1D Sequence of samples
- op     - The operator on the window R^N->R
- window - The length (N) of the window to sample from
"""
function windowedOperator{T}(data::Vector{T}, op, window::Int)
    @assert isodd(window)
    M = int((window-1)/2)
    N = length(data)
    access(x) = x<1?1:x>N?N:x
    result = zeros(N)
    for i=1:N
        input = zeros(window)
        for j=1:window
            input[j] = data[access(i+j-M)]
        end
        result[i] = op(input)
    end
    result

end

function likelyHoodSummaryPlot(title_string::ASCIIString, d1r, d2r, d3r,
    ll::Vector{Int}, figureId::Int)
    figure(figureId)
    plt.clf();
    title(title_string)
    println("red mean   = ", mean(d1r[d1r.!=0]))
    println("blue mean  = ", mean(d2r[d2r.!=0]))
    println("green mean = ", mean(d3r[d3r.!=0]))
    scatter(1:N, 10(d1r./maximum(d1r)), marker=".", color="r")
    scatter(1:N, 10(d2r./maximum(d1r)), marker=".", color="g")
    scatter(1:N, 10(d3r./maximum(d1r)), marker=".", color="b")
    scatter(linspace(1,N,length(ll)), ll, marker="|", color="k")
end

function runPeakIter(cover::Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64}},
                interest::Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64}},
                S::Matrix{Float64},
                iterId::Int,
                plotId::Int)

    println("Initial Classification.(+$iterId)..")
    peaks = findPeakSeq2(mean(S,2).*median(S,2))[:]

    (d1,d1r,er1) = getFreqEst(I, cover[1], interest[1],  peaks)
    (d2,d2r,er2) = getFreqEst(I, cover[2], interest[2],  peaks)
    (d3,d3r,er3) = getFreqEst(I, cover[3], interest[3], peaks)

    N = length(d1)
    if(doPlot)
        likelyHoodSummaryPlot("Transition Likelyhood(+$iterId)",
                              d1r, d2r, d3r, ll, plotId)
    end

	S=hcat(NN(d1r),NN(d2r),NN(d3r));
    if(doPlot)
        figure(plotId+100);
        plt.clf();
        plot(mean(S,2).*median(S,2))
        scatter(linspace(1,N,length(ll)), ll, marker="|", color="k")
    end
    S
end

"""
Estimate Which Transitions are likely to be real changepoints
When in doubt oversegment
"""
function doLikelyHoodEst(SubjectID::Int, workingDir::ASCIIString, doPlot::Bool)
    ll::Vector{Int}    = readcsv("$workingDir/Dejunkedlabels$SubjectID.csv")
    I::Matrix{Float64} = PyPlot.imread("$workingDir/DejunkedSpectra$SubjectID.png")

    coverLow::Matrix{Float64} = PyPlot.imread("$workingDir/coverlow$SubjectID.png")
    interestLow::Matrix{Int}  = readcsv("$workingDir/interest$SubjectID-low.csv")
    coverMed::Matrix{Float64} = PyPlot.imread("$workingDir/covermed$SubjectID.png")
    interestMed::Matrix{Int}  = readcsv("$workingDir/interest$SubjectID-med.csv")
    coverHigh::Matrix{Float64}= PyPlot.imread("$workingDir/coverhigh$SubjectID.png")
    interestHigh::Matrix{Int} = readcsv("$workingDir/interest$SubjectID-high.csv")

    #Define Types
    d1::Vector{Float64}
    d1r::Vector{Int}
    er1::Vector{Float64}
    d2::Vector{Float64}
    d2r::Vector{Int}
    er2::Vector{Float64}
    d3::Vector{Float64}
    d3r::Vector{Int}
    er3::Vector{Float64}

    (d1,d1r,er1) = getFreqEst(I, coverLow,  interestLow)
    (d2,d2r,er2) = getFreqEst(I, coverMed,  interestMed)
    (d3,d3r,er3) = getFreqEst(I, coverHigh, interestHigh)
    N = length(d1)
    println("Initial Classification...")
    if(doPlot)
        likelyHoodSummaryPlot("Transition Likelyhood(Orig)",
                              d1r, d2r, d3r, ll, 405)
    end

    #Normalizing Operators
    NN(x) = x./maximum(x)
    F(x)  = NN(windowedOperator(x./maximum(x),maximum,5))


	s=hcat(F(d1r),F(d2r),F(d3r));
    if(doPlot)
        figure(505);plt.clf();plot(mean(s,2).*median(s,2))
        scatter(linspace(1,N,length(ll)), ll, marker="|", color="k")
    end

    println("Initial Classification.(+1)..")
    ##plt.close("all")
    (d1,d1r,er1) = getFreqEst(I, coverLow,  interestLow,  (d1r+d2r+d3r).>1mean(d1r))
    (d2,d2r,er2) = getFreqEst(I, coverMed,  interestMed,  (d1r+d2r+d3r).>1mean(d2r))
    (d3,d3r,er3) = getFreqEst(I, coverHigh, interestHigh, (d1r+d2r+d3r).>1mean(d3r))
    N = length(d1)
    if(doPlot)
        likelyHoodSummaryPlot("Transition Likelyhood(+1)",
                              d1r, d2r, d3r, ll, 406)
    end
	s=hcat(F(d1r),F(d2r),F(d3r));
    if(doPlot)
        figure(506);plt.clf();plot(mean(s,2).*median(s,2))
        scatter(linspace(1,N,length(ll)), ll, marker="|", color="k")
    end
    println("Initial Classification.(+2)..")
    (d1,d1r,er1) = getFreqEst(I, coverLow,  interestLow,  (d1r+d2r+d3r).>1mean(d1r))
    (d2,d2r,er2) = getFreqEst(I, coverMed,  interestMed,  (d1r+d2r+d3r).>1mean(d2r))
    (d3,d3r,er3) = getFreqEst(I, coverHigh, interestHigh, (d1r+d2r+d3r).>1mean(d3r))
    
    N = length(d1)
    if(doPlot)
        likelyHoodSummaryPlot("Transition Likelyhood(+2)",
                              d1r, d2r, d3r, ll, 407)
    end
	S::Matrix{Float64}=hcat(F(d1r),F(d2r),F(d3r));
    if(doPlot)
        figure(507);plt.clf();plot(mean(S,2).*median(S,2))
        scatter(linspace(1,N,length(ll)), ll, marker="|", color="k")
    end

    #Refine peak sequence
    Co = (coverLow,    coverMed,    coverHigh)
    In = (interestLow, interestMed, interestHigh)

    S = runPeakIter(Co, In, S, 3, 408)
    S = runPeakIter(Co, In, S, 4, 409)
    S = runPeakIter(Co, In, S, 5, 410)
    
    #Last iteration inline
    peaks = findPeakSeq2(mean(S,2).*median(S,2))[:]
    println("Initial Classification.(+6)..")
    peaks[find((mean(S,2).*median(S,2)) .< 0.02)] = false
    (d1,d1r,er1) = getFreqEst(I, coverLow,  interestLow,  peaks)
    (d2,d2r,er2) = getFreqEst(I, coverMed,  interestMed,  peaks)
    (d3,d3r,er3) = getFreqEst(I, coverHigh, interestHigh, peaks)

    if(doPlot)
        likelyHoodSummaryPlot("Transition Likelyhood(+6)",
                              d1r, d2r, d3r, ll, 411)
    end

	S=hcat(NN(d1r),NN(d2r),NN(d3r));

    if(doPlot)
        figure(511);plt.clf();plot(mean(S,2).*median(S,2))
        title("Transition Likelyhood(+6)")
        scatter(linspace(1,N,length(ll)), ll, marker="|", color="k")
    end

    writecsv("$workingDir/likelyhood$SubjectID-low.csv", d1)
    writecsv("$workingDir/likelyhood$SubjectID-med.csv", d2)
    writecsv("$workingDir/likelyhood$SubjectID-high.csv",d3)
    writecsv("$workingDir/rawlikelyhood$SubjectID-low.csv", d1r)
    writecsv("$workingDir/rawlikelyhood$SubjectID-med.csv", d2r)
    writecsv("$workingDir/rawlikelyhood$SubjectID-high.csv",d3r)
    writecsv("$workingDir/edgeRate$SubjectID-low.csv", er1)
    writecsv("$workingDir/edgeRate$SubjectID-med.csv", er2)
    writecsv("$workingDir/edgeRate$SubjectID-high.csv",er3)
    writecsv("$workingDir/edgeParameter$SubjectID.csv", S)

    D1 = windowedOperator(d1r./maximum(d1r),maximum,5);
    D2 = windowedOperator(d2r./maximum(d1r),maximum,5);
    D3 = windowedOperator(d3r./maximum(d1r),maximum,5);
end
