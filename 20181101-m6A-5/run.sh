## Kidney Age
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_KI_age.txt _data/NK_m6Agenes.txt plots/ Kidney_age.NK_m6A
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_KI_age.txt _data/AK_m6Agenes.txt plots/ Kidney_age.AK_m6A
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_KI_age.txt _data/NK_m6Agenes.txt,_data/AK_m6Agenes.txt plots/ Kidney_age.AK_NK_m6A
python2 scripts/make.Revigo.Input.py plots/GO_result_table_Kidney_age.AK_NK_m6A.txt > Revigo/RevGO.Kidney_age.txt
sort -k2,2gr Revigo/RevGO.Kidney_age.txt | head -n 300 > Revigo2/RevGO.Kidney_age.txt
# put to Revigo server
Rscript scripts/REVIGO_plot_read_data.r Revigo/Kidney_Age.csv Revigo/Kidney_Age2.pdf
Rscript scripts/REVIGO_plot_read_data_2.r Revigo/Kidney_Age.csv Revigo2/Kidney_Age2.pdf


## Heart Age
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_V_age.txt _data/NH_m6Agenes.txt plots/ Heart_age.NH_m6A
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_V_age.txt _data/AH_m6Agenes.txt plots/ Heart_age.AH_m6A
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_V_age.txt _data/AH_m6Agenes.txt,_data/NH_m6Agenes.txt plots/ Heart_age.AH_NH_m6A
python2 scripts/make.Revigo.Input.py plots/GO_result_table_Heart_age.AH_NH_m6A.txt > Revigo/RevGO.Heart_age.txt
# put to Revigo server
Rscript scripts/REVIGO_plot_read_data.r Revigo/Heart_Age.csv Revigo/Heart_Age2.pdf
Rscript scripts/REVIGO_plot_read_data_2.r Revigo/Heart_Age.csv Revigo2/Heart_Age2.pdf



## Adult Kidney-Heart
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_adult_KI_V.txt _data/AK_m6Agenes.txt plots/ AHK.AK_m6A
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_adult_KI_V.txt _data/AH_m6Agenes.txt plots/ AHK.AH_m6A
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_adult_KI_V.txt _data/AH_m6Agenes.txt,_data/AK_m6Agenes.txt plots/ AHK.AH_AK_m6A
python2 scripts/make.Revigo.Input.py plots/GO_result_table_AHK.AH_AK_m6A.txt > Revigo/RevGO.AHK.txt
# put to Revigo server
Rscript scripts/REVIGO_plot_read_data.r Revigo/AHK.csv Revigo/AHK2.pdf
Rscript scripts/REVIGO_plot_read_data_2.r Revigo/AHK.csv Revigo2/AHK_2.pdf


## Neonatal Kidney-Heart
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_neonatal_KI_V.txt _data/NK_m6Agenes.txt plots/ NHK.NK_m6A
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_neonatal_KI_V.txt _data/NH_m6Agenes.txt plots/ NHK.NH_m6A
#Rscript scripts/GO_plot.R _data/DiffSiteGenes_neonatal_KI_V.txt _data/NH_m6Agenes.txt,_data/NK_m6Agenes.txt plots/ NHK.NH_NK_m6A
python2 scripts/make.Revigo.Input.py plots/GO_result_table_NHK.NH_NK_m6A.txt > Revigo/RevGO.NHK.txt
# put to Revigo server
Rscript scripts/REVIGO_plot_read_data.r Revigo/NHK.csv Revigo/NHK2.pdf
Rscript scripts/REVIGO_plot_read_data_2.r Revigo/NHK.csv Revigo2/NHK_2.pdf


# Adult Heart
#Rscript scripts/GO_plot.R _data/LV.top500.txt,_data/RV.top500.txt _data/AH_allExpressed.txt plots/ AHtop500.AH_allExpressed
#Rscript scripts/GO_plot.R _data/LV.top1000.txt,_data/RV.top1000.txt _data/AH_allExpressed.txt plots/ AHtop1000.AH_allExpressed
#Rscript scripts/GO_plot.R _data/LV.top2000.txt,_data/RV.top2000.txt _data/AH_allExpressed.txt plots/ AHtop2000.AH_allExpressed
python2 scripts/make.Revigo.Input.py plots/GO_result_table_AHtop1000.AH_allExpressed.txt > Revigo/RevGO.AH_top1000.txt
python2 scripts/make.Revigo.Input.py plots/GO_result_table_AHtop2000.AH_allExpressed.txt > Revigo/RevGO.AH_top2000.txt

Rscript scripts/REVIGO_plot_read_data_2.r Revigo/AH_top1000.csv Revigo2/AH_2.pdf


# Adult Kidney
#Rscript scripts/GO_plot.R _data/KI.top500.txt _data/AK_allExpressed.txt plots/ AKtop500.AK_allExpressed
#Rscript scripts/GO_plot.R _data/KI.top1000.txt _data/AK_allExpressed.txt plots/ AKtop1000.AK_allExpressed
#Rscript scripts/GO_plot.R _data/KI.top2000.txt _data/AK_allExpressed.txt plots/ AKtop2000.AK_allExpressed
python2 scripts/make.Revigo.Input.py plots/GO_result_table_AKtop1000.AK_allExpressed.txt > Revigo/RevGO.AK_top1000.txt

python2 scripts/make.Revigo.Input.py plots/GO_result_table_AKtop2000.AK_allExpressed.txt > Revigo/RevGO.AK_top2000.txt

Rscript scripts/REVIGO_plot_read_data_2.r Revigo/AK_top1000.csv Revigo2/AK_2.pdf


# Neonatal Heart
#Rscript scripts/GO_plot.R _data/NH.top500.txt _data/NH_allExpressed.txt plots/ NHtop500.NH_allExpressed
#Rscript scripts/GO_plot.R _data/NH.top1000.txt _data/NH_allExpressed.txt plots/ NHtop1000.NH_allExpressed
#Rscript scripts/GO_plot.R _data/NH.top2000.txt _data/NH_allExpressed.txt plots/ NHtop2000.NH_allExpressed
python2 scripts/make.Revigo.Input.py plots/GO_result_table_NHtop1000.NH_allExpressed.txt > Revigo/RevGO.NH_top1000.txt

python2 scripts/make.Revigo.Input.py plots/GO_result_table_NHtop2000.NH_allExpressed.txt > Revigo/RevGO.NH_top2000.txt

Rscript scripts/REVIGO_plot_read_data_2.r Revigo/NH_top1000.csv Revigo2/NH_2.pdf


# Neonatal Kidney
#Rscript scripts/GO_plot.R _data/P1_NKI.top500.txt _data/NK_allExpressed.txt plots/ NKtop500.NK_allExpressed
#Rscript scripts/GO_plot.R _data/P1_NKI.top1000.txt _data/NK_allExpressed.txt plots/ NKtop1000.NK_allExpressed
#Rscript scripts/GO_plot.R _data/P1_NKI.top2000.txt _data/NK_allExpressed.txt plots/ NKtop2000.NK_allExpressed
python2 scripts/make.Revigo.Input.py plots/GO_result_table_NKtop1000.NK_allExpressed.txt > Revigo/RevGO.NK_top1000.txt

python2 scripts/make.Revigo.Input.py plots/GO_result_table_NKtop2000.NK_allExpressed.txt > Revigo/RevGO.NK_top2000.txt

Rscript scripts/REVIGO_plot_read_data_2.r Revigo/NK_top1000.csv Revigo2/NK_2.pdf
