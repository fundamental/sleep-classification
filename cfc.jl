using DSP
    
"""
Calculate Kulback-Leiber Divergence from uniform distribution
"""
function kl_uniform(x::Vector{Float64})
    N=length(x)
    expect = 1.0/N
    x += 1e-6
    xx = x/sum(x)
    sum(xx.*log(xx/expect))
end


"""
Extract the level of cross frequency coupling between two frequency bands
#d = data, fs=sampling rate,l=low band, h=high band
"""
function extract_cfc(x::Vector{Float64}, fs::Float64, l::Vector{Float64}, h::Vector{Float64})
    SN = length(x)#4096*32
    window_sec = 30
    period = 200*window_sec

    fl1 = digitalfilter(Bandpass((h*2.0./fs)...), Butterworth(16))
    fl2 = digitalfilter(Bandpass((l*2.0./fs)...), Butterworth(16))

    y1::Vector{Float64} = filt(fl1, x)
    y2::Vector{Float64} = filt(fl2, x)

    Y1::Vector{Float64} = abs(hilbert(y1))
    Y2::Vector{Float64} = angle(hilbert(y2))

    N::Int     = 256
    M          = floor(Int, SN/period)
    bins       = zeros(N,M)
    bins_count = zeros(N,M)
    vals       = N*(Y2+pi)/(2pi)
    for i=1:M, j= 1:period
        idx = (i-1)*period + j
        bin = floor(Int, vals[idx])+1


        bins[bin ,i]      += Y1[idx]
        bins_count[bin,i] += 1
    end

    bins_count[bins_count .== 0] = 1

    b = bins./bins_count

    out::Vector{Float64} = mapslices(kl_uniform, b, 1)[:]
end

"""
Generate CFC Vectors
"""
function make_features(stream::Vector{Float64}; window::Int=30,
    in_freq::Float64=200.0)
    freqs = [
    3.0  4.0
    4.0  8.0
    8.0  13.0
    13.0 16.0
    16.0 30.0
    30.0 58.0
    ]

    cfc_feats = Any[]
    Nfreq = size(freqs,1)
    for i=1:Nfreq, j=(i+1):Nfreq
    #for i=1:Nfreq, j=1:Nfreq
        println("i = $i, j=$j")
        push!(cfc_feats, extract_cfc(stream, 200.0, freqs[i,:], freqs[j,:]))
    end
    hcat(cfc_feats...)'
end


#cfc_out = cfc(d, 200.0, [0.1,10.0], [30.0,60.0])



#function extract_band(stream::Vector{Float64}, chunk::Int, Fs::Int, Fl::Float64, Fh::Float64)
#    #Find samples per chunk
#    smp_per_chunk = chunk*Fs;
#    #Find total number of chunks per recording
#    num_chunks    = floor(Int,length(stream)/smp_per_chunk)
#
#    #Find bin ids which correspond to the frequency range
#    bin_low  = round(Int,smp_per_chunk*Fl/(2Fs))
#    bin_high = round(Int,smp_per_chunk*Fh/(2Fs))
#
#    out = zeros(num_chunks)
#    for i=1:num_chunks
#        c = stream[smp_per_chunk*(i-1) + (1:smp_per_chunk)]
#        C = fft(c)
#        for j=bin_low:bin_high
#            out[i] += abs(C[j+1].^2)
#        end
#    end
#    out
#end

function demo_cfc()
    data     = randn(1000000)
    cfc_feat = make_features(data)
    figure(1)
    imshow(cfc_feat, interpolation="none", aspect="auto")
    cfc_feat
end
