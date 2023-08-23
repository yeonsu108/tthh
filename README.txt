[Work Flow] 
___Cutflow and Histogram___
1. Object Selection > "ana_tthh.C" -> Output : ./sample1/
2. Event Selection > "drawHisto.py" -> Output : ./plots/OS/
                     "drawHisto_Same.py" -> Output : ./plots/OS/Same/ 
3. Cutflow > Cutflow.py

___DNN_____________________
1. Object Selection > "ana_tthh.C" -> Output : ./sample1/
2. Add Branches > "DNN_variable.py" -> Output : ./sample2/
3. DNN > train.py -> Output : ./DNN_result/OS/


[Updates]
- Add files : classes/,  external/, drawHisto_Same.py, Cutflow.py, DNN_variable.py
- Change file names : drawHisto.py -> drawHisto_Same.py (tthh, ttbbbb, ttbbcc, ttbb all in one canvas.)
                      and NOW, drawHisto.py makes single histograms.

- Newly defined functions :
               utility.h -> ::MakeLastTag, ::FromMotherExact, ::FromWhere, ::dR, ::ConcatFloat, ::RecoHiggs, ::Order
                            ::ConcatLep -> ::ConcatVector 

- Newly defined branches (main) : Higgs1,2 and there bJets. 

- Sample directory : ./sample1/ & ./sample2/

- Tidy up some codes :

               Unify Input/Output naming : You may only change "OS" tag. 

               ana_tthh_df.C -> df0 : Cconstants, Basic Branches
                                df1 : Gen_bi
                                df2 : nGen
                                df3 : Gen 4-vector
                                df4 : Gen dR
                                df5 : Reco
                                Sorting : Common -> From -> Add && Quark -> Jet

