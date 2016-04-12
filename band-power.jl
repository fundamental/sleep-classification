using HDF5
using PyPlot
using DecisionTree

"""
Import a EEG datastream

Example:
import_data(5) #Import 5th subject
"""
function import_data(SubjectID)
    f = HDF5.h5open("/home/mark/current/general-sleep/subject$SubjectID.h5")
    DataSet = read(f["2"])
    close(f)
    println("Dataset is ", length(DataSet), " elements")
    convert(Vector{Float64},DataSet)
end

"""
Import a EEG labeling

Example:
import_labels(5) #Import 5th subject
"""
function import_labels(SubjectID)
    convert(Vector{Int},readcsv("/home/mark/current/general-sleep/HypnogramAASM_subject$SubjectID.txt")[2:end])
end

"""
Extract Energy from frequency band

Arguments:
stream - EEG datastream
chunk  - length of window (seconds)
Fs     - sampling rate    (Hz)
Fl     - low  frequency bin
Fh     - high frequency bin
"""
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

"""
Generate Bandpower features for an EEG datastream

    Arguments:

    stream  - EEG data
    window  - fixed time window (seconds)
    in_freq - sample rate of EEG data (Hz)
"""
function make_features(stream::Vector{Float64}; window::Int=30,
    in_freq::Float64=200.0)
    Delta     = extract_band(stream, window, in_freq, 3., 4.);
    Theta     = extract_band(stream, window, in_freq, 4., 8.);
    Alpha     = extract_band(stream, window, in_freq, 8., 13.);
    Low_Beta  = extract_band(stream, window, in_freq, 13.,16.);
    High_Beta = extract_band(stream, window, in_freq, 16.,30.);
    Gamma     = extract_band(stream, window, in_freq, 30.,58.)
    hcat(Delta, Theta, Alpha, Low_Beta, High_Beta, Gamma)'
end

"""
Load all data and extract features and labels
"""
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
