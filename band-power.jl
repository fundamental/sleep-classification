using HDF5
using PyPlot
using DecisionTree

function import_data(SubjectID)
    f = HDF5.h5open("/home/mark/current/general-sleep/subject$SubjectID.h5")
    DataSet = read(f["2"])
    close(f)
    println("Dataset is ", length(DataSet), " elements")
    convert(Vector{Float64},DataSet)
end

function import_labels(SubjectID)
    convert(Vector{Int},readcsv("/home/mark/current/general-sleep/HypnogramAASM_subject$SubjectID.txt")[2:end])
end

function extract_band(stream::Vector{Float64}, chunk::Int, Fs::Int, Fl::Float64, Fh::Float64)
    #Find samples per chunk
    smp_per_chunk = chunk*Fs;
    #Find total number of chunks per recording
    num_chunks    = floor(Int,length(stream)/smp_per_chunk)

    #Find bin ids which correspond to the frequency range
    bin_low  = round(Int,smp_per_chunk*Fl/(2Fs))
    bin_high = round(Int,smp_per_chunk*Fh/(2Fs))

    out = zeros(num_chunks)
    for i=1:num_chunks
        c = stream[smp_per_chunk*(i-1) + (1:smp_per_chunk)]
        C = fft(c)
        for j=bin_low:bin_high
            out[i] += abs(C[j+1].^2)
        end
    end
    out
end

function make_features(stream, window=30)
    Delta     = extract_band(stream, 30, 200, 3., 4.);
    Theta     = extract_band(stream, 30, 200, 4., 8.);
    Alpha     = extract_band(stream, 30, 200, 8., 13.);
    Low_Beta  = extract_band(stream, 30, 200, 13.,16.);
    High_Beta = extract_band(stream, 30, 200, 16.,30.);
    Gamma     = extract_band(stream, 30, 200, 30.,58.)
    hcat(Delta, Theta, Alpha, Low_Beta, High_Beta, Gamma)'
end

function validate(train, train_label, test, test_label)
    model = build_forest(train_label, train', 4,20)#, 200, 20, 0.5)
    out = apply_forest(model, test')
    println("classification accuracy = ", mean(out.==test_label))
    out
end

function load_all()
    C = Any[];for i=1:20; push!(C, make_features(import_data(i)));end
    L = Any[];for i=1:20; push!(L, import_labels(i)[1:6:end]);end
    #fixup L
    for i=1:20;
        if(length(L[i]) != size(C[i],2))
            L[i] = L[i][1:size(C[i],2)];
        end;
    end

    (C,L)
end


function make_conf(A,B,values)
    N    = length(values)
    conf = zeros(N,N)
    rm   = Dict{Any, Int}() #Remapper
    for i=1:N
        rm[values[i]] = i
    end

    M = length(A)
    for j=1:M
        a = A[j]
        b = B[j]
        if(a in values && b in values)
            conf[rm[a], rm[b]] += 1
        end
    end
    conf
end
