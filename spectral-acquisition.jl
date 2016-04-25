#Input LSpectra#ID.png
using Images
using PyPlot
using HDF5
using FixedPointNumbers

"""
Import Data from a HDF5 Serialization

These files are derived the original .edf files supplied from the DREAMS
database
"""
function import_data(SubjectID::Int)
    f = HDF5.h5open("/home/mark/current/general-sleep/subject$SubjectID.h5")
    DataSet::Vector{Float64} = read(f["2"])
    close(f)
    println("Dataset is ", length(DataSet), " elements")
    DataSet
end

function import_data(SubjectID::Int)
    f = HDF5.h5open("/home/mark/current/general-sleep/patient$SubjectID.h5")
    DataSet::Vector{Float64} = read(f["2"])
    close(f)
    println("Dataset is ", length(DataSet), " elements")
    DataSet
end

"""
Obtain the raw spectra from SubjectID
"""
function raw_spec(SubjectID::Int)
    DD = import_data(SubjectID)
    figure(999)
    (Spectra,_) = specgram(DD, 4096, 100, noverlap=0)
    PyPlot.close()
    Spectra = log(abs(Spectra.^2))
end

"""
Generate the spectral image
"""
function generate_spectra(DataSet::Vector{Float64}, SubjectID::Int, doPlot::Bool)
    figure(101)
    (Spectra,_) = specgram(DataSet, 4096, 100, noverlap=0)
    Spectra = log(abs(Spectra))
    title("Unnormalized Spectra")
    imshow(Spectra, aspect="auto", interpolation="none")
    if(!doPlot)
        PyPlot.close()
        Base.yield()
    end
    Spectra = mapslices(x->x-sort(x)[floor(Int,end/2)], Spectra, 1)
    Spectra = mapslices(x->x-sort(x)[floor(Int,end/2)], Spectra, 2)
    Spectra -= minimum(Spectra)
    Spectra /= mean(Spectra)
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
        PyPlot.clf()
        imshow(Spectra[100:end,:], aspect="auto", interpolation="none",
        vmin=Smin, vmax=Smax, cmap=gray())
        title("Normalized Spectra")
    end
    Images.imwrite(convert(Array{Ufixed16},Spectra./maximum(Spectra)), "Spectra$SubjectID.png")
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


"""
Eliminate series of outliers in the data
"""
function robustSpikeElimination(Spectra::Matrix{Float64})
    excitation = sum((Spectra[:,1:end-1].-Spectra[:,2:end]).^2,1)[:]
    ex = refold(excitation)

    minmean = Inf
    minstd  = Inf
    minelms = -1
    tmp = Float64[]
    #Run a min variance alg to find a suitable gausian description of this
    #spectral distance
    N = length(excitation)
    for i=round(Int,0.8N):N
        samples = excitation[1:i]

        m = mean(samples)
        s = std(samples)
        if(s<minstd)
            minstd  = s
            minmean = m
            minelms = i
        end
        push!(tmp, s)
    end
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


"""
Create Image for the next stage of processing

Arguments:

- SubjectID  - the identifying number for the subject
- workingDir - the directory to save intermediate data to
- doPlot     - plot intermediate results
"""
function runAquisition(SubjectID::Int, workingDir::ASCIIString, doPlot::Bool)
    Data = import_data(SubjectID)
    (Spectra, LSpectra) = generate_spectra(Data, SubjectID, doPlot)
    (goodElms,seq) = robustSpikeElimination(Spectra)

    if(doPlot)
        figure(103)
        title("Dejunked Spectra")
        imshow(Spectra[1:end,goodElms], aspect="auto", interpolation="none",
        cmap=gray())
    end

    Images.imwrite(convert(Array{Ufixed16},Spectra[:,goodElms]./maximum(Spectra[:,goodElms])),
    "$workingDir/DejunkedSpectra$SubjectID.png")
    writecsv("$workingDir/Dejunked$SubjectID-cols.csv", seq)

    #raw_labels = readcsv("general-sleep/HypnogramAASM_subject$SubjectID.txt")[2:end]
    raw_labels = readcsv("/home/mark/current/general-sleep/patients/HypnogramAASM_patient$SubjectID.txt")[2:end]
    println("length(seq)       =", length(seq))
    println("length(raw_labels)=", length(raw_labels))
    dejunked_labels = raw_labels[round(Int, find(seq)*4.096)]
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
