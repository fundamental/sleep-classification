#feat_ds_dds = Matrix{Float64}[]
#labl_ds_dds = Vector{Int}[]
#println("Loading Subject Features...")
#for i=1:20
#    stuff = viewStuff(string(i), "/home/mark/current/work/", 0.75,
#                      dataDir="/home/mark/current/general-sleep/subject")
#    push!(labl_ds_dds, stuff[1])
#    push!(feat_ds_dds, stuff[2])
#end

println("Running Trials...")
if(false)
#Iterate over training amount
for i=1:19
    println("Evaluating performance with $i training examples")
    tmp = []
    for j=1:4
        print("Getting sample...")
        ordering   = sortperm(rand(20))
        s1         = ordering[1:i]
        s2         = ordering[(i+1):end]
        if(i == 1)
            labl_train = labl_ds_dds[s1[1]]
            feat_train = feat_ds_dds[s1[1]]'
            model = build_forest(labl_train, feat_train, 20, 40)
        else
            labl_train = vcat(labl_ds_dds[s1]...)
            feat_train = hcat(feat_ds_dds[s1]...)'
            model = build_forest(labl_train, feat_train, 20, 40)
        end
         
        res   = []
        for k=1:length(s2)
            push!(res,mean(labl_ds_dds[s2[k]].==apply_forest(model,feat_ds_dds[s2[k]]')))
        end
        println("(",mean(res),")")
        append!(tmp,res)
    end
    println("performance at $i = ", median(tmp),"/", mean(tmp))
end
end

#Running Trials...
#Evaluating performance with 1 training examples
#Getting sample...(0.6198845076530749)
#Getting sample...(0.5630149056435356)
#Getting sample...(0.6179823400818494)
#Getting sample...(0.6532679597784167)
#performance at 1 = 0.6316185968324387/0.6135374282892192
#Evaluating performance with 2 training examples
#Getting sample...(0.6962847443225395)
#Getting sample...(0.6568313289151546)
#Getting sample...(0.6565250133945512)
#Getting sample...(0.6932294236190232)
#performance at 2 = 0.6911747037670839/0.6757176275628172
#Evaluating performance with 3 training examples
#Getting sample...(0.7058486295274878)
#Getting sample...(0.6898935429209869)
#Getting sample...(0.6932906570580606)
#Getting sample...(0.6581785720064025)
#performance at 3 = 0.6989956751899126/0.6868028503782345
#Evaluating performance with 4 training examples
#Getting sample...(0.7263191208847743)
#Getting sample...(0.6914153127301494)
#Getting sample...(0.7092682709081103)
#Getting sample...(0.6859194467409561)
#performance at 4 = 0.6999956217875568/0.7032305378159978
#Evaluating performance with 5 training examples
#Getting sample...(0.7033777770437852)
#Getting sample...(0.7208085673357788)
#Getting sample...(0.7269341580089781)
#Getting sample...(0.6995895829272406)
#performance at 5 = 0.7306831509716882/0.7126775213289457
#Evaluating performance with 6 training examples
#Getting sample...(0.710086021172373)
#Getting sample...(0.754580959876837)
#Getting sample...(0.7311220725679252)
#Getting sample...(0.7327628394629672)
#performance at 6 = 0.7498222680364516/0.7321379732700255
#Evaluating performance with 7 training examples
#Getting sample...(0.7400567641149376)
#Getting sample...(0.7366673431711797)
#Getting sample...(0.7435306657212665)
#Getting sample...(0.7207846681955268)
#performance at 7 = 0.735810658335766/0.7352598603007275
#Evaluating performance with 8 training examples
#Getting sample...(0.717401143235679)
#Getting sample...(0.7388606632915856)
#Getting sample...(0.7294460623406679)
#Getting sample...(0.6836644367907344)
#performance at 8 = 0.7417193342063808/0.7173430764146668
#Evaluating performance with 9 training examples
#Getting sample...(0.7624440508225612)
#Getting sample...(0.7756688223524638)
#Getting sample...(0.7325284011681183)
#Getting sample...(0.7381807375449491)
#performance at 9 = 0.7558936609793807/0.752205502972023
#Evaluating performance with 10 training examples
#Getting sample...(0.6921451715655043)
#Getting sample...(0.7287775938280594)
#Getting sample...(0.7436374071065994)
#Getting sample...(0.7292588452973318)
#performance at 10 = 0.739384248451606/0.7234547544493739
#Evaluating performance with 11 training examples
#Getting sample...(0.729269747836881)
#Getting sample...(0.7625101493770047)
#Getting sample...(0.7365278118265408)
#Getting sample...(0.7509314178584559)
#performance at 11 = 0.7439758756453589/0.7448097817247205
#Evaluating performance with 12 training examples
#Getting sample...(0.7709035552821426)
#Getting sample...(0.7259049259509058)
#Getting sample...(0.7293319263936117)
#Getting sample...(0.7316010998799899)
#performance at 12 = 0.7322670341657251/0.7394353768766623
#Evaluating performance with 13 training examples
#Getting sample...(0.7167260365871245)
#Getting sample...(0.7764839368048084)
#Getting sample...(0.7707636399514775)
#Getting sample...(0.7576067708454283)
#performance at 13 = 0.7719315707586969/0.7553950960472096
#Evaluating performance with 14 training examples
#Getting sample...(0.7158438016626115)
#Getting sample...(0.7512791554239809)
#Getting sample...(0.7280121701900861)
#Getting sample...(0.7279673075111388)
#performance at 14 = 0.7375158472748433/0.7307756086969542
#Evaluating performance with 15 training examples
#Getting sample...(0.7771126951031037)
#Getting sample...(0.7540229781880521)
#Getting sample...(0.7488574924845237)
#Getting sample...(0.7726572282249251)
#performance at 15 = 0.779299826422519/0.7631625985001512
#Evaluating performance with 16 training examples
#Getting sample...(0.747779988909929)
#Getting sample...(0.7287128836923135)
#Getting sample...(0.7571820911731901)
#Getting sample...(0.7552597451635672)
#performance at 16 = 0.7434713503943344/0.74723367723475
#Evaluating performance with 17 training examples
#Getting sample...(0.6854588475755755)
#Getting sample...(0.745118377255244)
#Getting sample...(0.806822296515031)
#Getting sample...(0.7602695358190431)
#performance at 17 = 0.7729342875731945/0.7494172642912235
#Evaluating performance with 18 training examples
#Getting sample...(0.7578443962765562)
#Getting sample...(0.8084690280648419)
#Getting sample...(0.7086672489289129)
#Getting sample...(0.7508423673268325)
#performance at 18 = 0.7508423673268325/0.7564557601492858
#Evaluating performance with 19 training examples
#Getting sample...(0.7270624518118736)
#Getting sample...(0.7059585492227979)
#Getting sample...(0.8065902578796562)
#Getting sample...(0.8324761204996326)
#performance at 19 = 0.7668263548457649/0.7680218448534901

acc_samples = [
1  0.6198845076530749
1  0.5630149056435356
1  0.6179823400818494
1  0.6532679597784167
2  0.6962847443225395
2  0.6568313289151546
2  0.6565250133945512
2  0.6932294236190232
3  0.7058486295274878
3  0.6898935429209869
3  0.6932906570580606
3  0.6581785720064025
4 0.7263191208847743
4 0.6914153127301494
4 0.7092682709081103
4 0.6859194467409561
5 0.7033777770437852
5 0.7208085673357788
5 0.7269341580089781
5 0.6995895829272406
6 0.710086021172373
6 0.754580959876837
6 0.731122072567925
6 0.732762839462967
7 0.7400567641149376
7 0.7366673431711797
7 0.7435306657212665
7 0.7207846681955268
8 0.717401143235679
8 0.738860663291585
8 0.729446062340667
8 0.683664436790734
9 0.7624440508225612
9 0.7756688223524638
9 0.7325284011681183
9 0.7381807375449491
10 0.6921451715655043
10 0.7287775938280594
10 0.7436374071065994
10 0.7292588452973318
11 0.729269747836881
11 0.762510149377004
11 0.736527811826540
11 0.750931417858455
12 0.7709035552821426
12 0.7259049259509058
12 0.7293319263936117
12 0.7316010998799899
13 0.7167260365871245
13 0.7764839368048084
13 0.7707636399514775
13 0.7576067708454283
14 0.7158438016626115
14 0.7512791554239809
14 0.7280121701900861
14 0.7279673075111388
15 0.7771126951031037
15 0.7540229781880521
15 0.7488574924845237
15 0.7726572282249251
16 0.747779988909929
16 0.728712883692313
16 0.757182091173190
16 0.755259745163567
17 0.685458847575575
17 0.745118377255244
17 0.806822296515031
17 0.760269535819043
18 0.7578443962765562
18 0.8084690280648419
18 0.7086672489289129
18 0.7508423673268325
19 0.7270624518118736
19 0.7059585492227979
19 0.8065902578796562
19 0.8324761204996326]
