#Input LSpectra#ID.png
using Images
using PyPlot
using HDF5
using FixedPointNumbers

#Load The Data
function import_data(SubjectID)
    f = HDF5.h5open("/home/mark/current/general-sleep/subject$SubjectID.h5")
    DataSet = read(f["2"])
    close(f)
    println("Dataset is ", length(DataSet), " elements")
    DataSet
end

function raw_spec(SubjectID)
    DD = import_data(SubjectID)
    figure(999)
    (Spectra,_) = specgram(DD, 4096, 100, noverlap=0)
    plt.close()
    Spectra = log(abs(Spectra.^2))
end

#Generate The Image
function generate_spectra(DataSet, SubjectID, doPlot)
    figure(101)
    (Spectra,_) = specgram(DataSet, 4096, 100, noverlap=0)
    Spectra = log(abs(Spectra))
    title("Unnormalized Spectra")
    imshow(Spectra, aspect="auto", interpolation="none")
    if(!doPlot)
        plt.close()
        Base.yield()
    end
    Spectra = mapslices(x->x-sort(x)[floor(Int,end/2)], Spectra, 1)
    Spectra = mapslices(x->x-sort(x)[floor(Int,end/2)], Spectra, 2)
    #Spectra = mapslices(x->x./norm(sort(x)[floor(end/2)]), Spectra, 2)
    Spectra -= minimum(Spectra)
    Spectra /= mean(Spectra)
    #Spectra = log(abs(Spectra))
    #Spectra = downsample(Spectra, 8)
    #spectra = mapslices(x->sort(x), spectra, 2)
    LSpectra = log(Spectra)
    if(sum(LSpectra.<-10) > 0)
        LSpectra[LSpectra.<-10] = minimum(LSpectra[LSpectra.>-10][:])
    end
    Lmean = mean(LSpectra)
    Lstd  = std(LSpectra)
    Lmin  = mean(LSpectra)-3Lstd
    Lmax  = mean(LSpectra)+3Lstd
    Smean = mean(Spectra)
    Sstd  = std(Spectra)
    Smin  = mean(Spectra)-3Sstd
    Smax  = mean(Spectra)+3Sstd

    if(doPlot)
        figure(102)
        plt.clf()
        imshow(Spectra[100:end,:], aspect="auto", interpolation="none",
        vmin=Smin, vmax=Smax, cmap=gray())
        title("Normalized Spectra")
    end
    #figure(102)
    #plt.clf()
    #imshow(LSpectra[100:end,:], aspect="auto", interpolation="none",
    #vmin=Lmin, vmax=Lmax, cmap=gray())
    #title("Normalized Log Spectra")
    #figure(103)
    #plt.clf()
    #plt.hist(Spectra[:],256)
    #figure(104)
    #plt.clf()
    #plt.hist(LSpectra[:],256)
    Images.imwrite(convert(Array{Ufixed16},Spectra./maximum(Spectra)), "Spectra$SubjectID.png")
    #Images.imwrite(log(Spectra), "LSpectra$SubjectID.png")
    (Spectra, LSpectra)
end

#Identify Spikes
MaxSpikeLen = 3
function refold(data)
    t1 = sort(data)
    tout = Float64[]
    mid = ceil(Int,length(t1)/2)
    for i=1:mid
        if(i+mid < length(t1))
            push!(tout, t1[mid+i])
        end
        if(i-mid > 0)
            push!(tout, t1[mid-i])
        end
    end
    tout
end
function robustSpikeElimination(Spectra)
    excitation = sum((Spectra[:,1:end-1].-Spectra[:,2:end]).^2,1)[:]
    ex = refold(excitation)
    #figure(105);
    #plt.clf()
    #plot(excitation)

    minmean = Inf
    minstd  = Inf
    minelms = -1
    tmp = Float64[]
    #Run a min variance alg to find a suitable gausian description of this
    #spectral distance
    N = length(excitation)
    for i=int(0.8N):N
        samples = excitation[1:i]

        m = mean(samples)
        s = std(samples)
        if(s<minstd)
            minstd  = s
            minmean = m
            minelms = i
        end
        push!(tmp, s)
        #println(i, " m=", m, " s=", s, " ?=",sum(samples.^2)/2-(sum(samples)/2).^2)
    end
    #figure(105);
    #plot([1:N],ones(N)*(minmean+1.5minstd), color="red")
    #figure(106);
    #plt.clf();
    #plot(tmp)
    #figure(107);
    #plt.clf()
    #plot(excitation)
    #plt.hist(excitation,128)
    println("Excluding ", sum(excitation.>(minmean+1.5minstd)), " Elements...")
    println("Of Total Possible ", size(Spectra)[2], "...")
    mixed = zeros(N)
    good = Int[]
    for i=1:N
        if(excitation[i] .< (minmean+1.5minstd) && 
            (excitation[max(1,i-1)] .< (minmean+1.5minstd) ||
             excitation[min(N,i+1)] .< (minmean+1.5minstd)))
            push!(good, i)
        end
    end
    mixed[good] = 1
    (good,mixed)
end


function runAquisition(SubjectID, workingDir, doPlot)
    Data = import_data(SubjectID)
    (Spectra, LSpectra) = generate_spectra(Data, SubjectID, doPlot)
    (goodElms,seq) = robustSpikeElimination(Spectra)

    if(doPlot)
        figure(103)
        title("Dejunked Spectra")
        imshow(Spectra[1:end,goodElms], aspect="auto", interpolation="none",
        cmap=gray())
    end

    #figure(109)
    #plt.clf()
    #imshow(Spectra, aspect="auto", interpolation="none",
    #cmap=gray())

    Images.imwrite(convert(Array{Ufixed16},Spectra[:,goodElms]./maximum(Spectra[:,goodElms])),
    "$workingDir/DejunkedSpectra$SubjectID.png")
    writecsv("$workingDir/Dejunked$SubjectID-cols.csv", seq)

    raw_labels = readcsv("general-sleep/HypnogramAASM_subject$SubjectID.txt")[2:end]
    println("length(seq)       =", length(seq))
    println("length(raw_labels)=", length(raw_labels))
    dejunked_labels = raw_labels[int(find(seq)*4.096)]
    if(doPlot)
        figure(104)
        title("Raw Labels")
        plot(raw_labels)
        figure(105)
        title("Dejunked Labels")
        plot(dejunked_labels)
    end
    writecsv("$workingDir/Dejunkedlabels$SubjectID.csv", dejunked_labels)
end

#Run The Results Through the sobal operator

#Save

#Manually Tweak With Gimp
# equalize
# median
# normalization
# something else?

