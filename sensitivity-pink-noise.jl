#Setup experiment parameters
dB = 1
max_snr =  60dB
min_snr =  20dB
int_snr =  2.5dB
drange  =  data_range()[38:end]
sub_id  =  16 #ST7192
subject =  drange[sub_id]

function modified_spec_aquire(SubjectID, Data, workingDir, doPlot)
    (Spectra, LSpectra) = generate_spectra(Data, SubjectID, doPlot)
    (goodElms,seq) = robustSpikeElimination(Spectra)

    if(doPlot)
        figure(103)
        title("Dejunked Spectra")
        imshow(Spectra[1:end,goodElms], aspect="auto", interpolation="none",
        cmap=gray())
    end

    ss  = Spectra[:,goodElms]
    ss -= minimum(ss)
    ss  = ss./maximum(ss)
    Images.imwrite(convert(Array{Ufixed16},ss),
    "$workingDir/DejunkedSpectra$SubjectID.png")
    writecsv("$workingDir/Dejunked$SubjectID-cols.csv", seq)

    #raw_labels = readcsv("general-sleep/HypnogramAASM_subject$SubjectID.txt")[2:end]
    #raw_labels = readcsv("/home/mark/current/general-sleep/patients/HypnogramAASM_patient$SubjectID.txt")[2:end]
    raw_labels = get_raw_labels(subject)
    println("length(seq)       =", length(seq))
    println("length(raw_labels)=", length(raw_labels))
    #dejunked_labels = raw_labels[round(Int, find(seq)*4.096)]
    #inds = linspace(1,length(raw_labels),size(Spectra,2))[find(seq)]

    ll = raw_labels
    gs = goodElms
    sec_dat = size(Spectra,2)*4096/200
    sec_lab = length(ll)*10
    println("   [INFO] mapping ", sec_lab, " sec to ", sec_dat, " sec")
    println("   [INFO] ", size(gs), " vs ", size(ss))
    nll = zeros(Int, size(Spectra,2))
    for i=1:length(nll)
        sec    = (i-0.5)*4096/200
        llpos  = round(Int, sec/10)
        if(llpos == 0)
            nll[i] = 5
        elseif(llpos <= length(ll))
            nll[i] = ll[llpos]
        else
            nll[i] = 5
        end
    end
    dejunked_labels = nll[gs];
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


#For each value in the snr range:
# - generate the noisy input
# - save the noisy input
# - process the noisy input
# - save the processed noisy output
data = import_data(subject)
base = std(data)
if(false)
    for i=min_snr:int_snr:max_snr
        if(i==20 || i==30 || i==40 || i==50 || i==60)
            continue
        end
        noise_mult = base/(10^(i/20.0))
        println("cur_snr    = ", i, "dB")
        println("base       = ", base)
        println("noise_mult = ", noise_mult)
        dd = data + noise_mult*pink_noise(randn(length(data))/3.06)
        PyPlot.close("all")
        writecsv("physionet-noise/$subject-noise$i-dB.csv", dd)
        workingDir = "physionet-denoise/"
        Data      = dd
        doPlot    = false
        SubjectID = string(subject, i)
        modified_spec_aquire(SubjectID, Data, "physionet-denoise/", true)
        #runAquisition(i, "physionet-denoise/", true)
        doMakeThresholds(SubjectID,            "physionet-denoise/", true)
        doRectSegment(SubjectID, "low",  100, "physionet-denoise/", true)
        doRectSegment(SubjectID, "med",  200, "physionet-denoise/", true)
        doRectSegment(SubjectID, "high", 300, "physionet-denoise/", true)
        doLikelyHoodEst(SubjectID,"physionet-denoise/", true)
        #viewStuff(SubjectID, "physionet-denoise/")
    end
end
    
return
include("misc.jl")

PyPlot.close("all")

#Supervised stage:
# - Load all non-experimental samples
# - Load all experimental samples
# - Generate the random forest model from training set
# - Perform Classification
# - save the classification result

FF = Matrix{Float64}[]
LL = Vector{Int}[]
if(true)
    for j=drange
        if(j != "ST7221" && j != subject)
            tmp = viewStuff(j, "physionet/", 0.7, false)
            push!(LL, map(Int, tmp[1][:]))
            push!(FF, tmp[2])
        end
    end
end

EX = Matrix{Float64}[]
EL = Vector{Int}[]
for i=min_snr:int_snr:max_snr
    println("i = $i")
    if(i == 20 || i == 30 || i == 40 || i == 50 || i== 60)
        i = round(Int, i)
    end
    println(subject)
    tmp = viewStuff(string(subject,i), "physionet-denoise/", 0.7, false,
    [string(i),subject])
    push!(EL, map(Int, tmp[1][:]))
    push!(EX, tmp[2])
end
#for i=1:length(EX)
#    figure(2000+i)
#    imshow(EX[i],aspect="auto",interpolation="none")
#end

figure(1010101)
model = build_forest(vcat(LL...), hcat(FF...)', 20, 40)
PyPlot.close("all")

class_result = Vector{Int}[]
for i=1:length(EX)
    tmp = apply_forest(model, EX[i]')
    println([(min_snr:int_snr:max_snr)...][i])
    println("classification accuracy = ", mean(EL[i].==tmp))
    push!(class_result, tmp)
end
