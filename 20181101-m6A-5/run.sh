## Kidney Age

# query enrichR web and export txt
# MDM
python2 scripts/make.Revigo.Input.py \
enrichr/GO_BP.Kidney_Age.MDM.txt \
> Revigo/RevGO.Kidney_Age.MDM.txt

# put to Revigo server

Rscript scripts/REVIGO_plot_read_data_2.r Revigo/Kidney_Age.MDM.csv Revigo/Kidney_Age.MDM.pdf

# EDM
python2 scripts/make.Revigo.Input.py \
enrichr/GO_BP.Kidney_Age.EDM.txt \
> Revigo/RevGO.Kidney_Age.EDM.txt

# put to Revigo server

Rscript scripts/REVIGO_plot_read_data_2.r Revigo/Kidney_Age.EDM.csv Revigo/Kidney_Age.EDM.pdf


## Heart Age

# MDM
python2 scripts/make.Revigo.Input.py \
enrichr/GO_BP.Heart_Age.MDM.txt \
> Revigo/RevGO.Heart_Age.MDM.txt

# put to Revigo server

Rscript scripts/REVIGO_plot_read_data_2.r Revigo/Heart_Age.MDM.csv Revigo/Heart_Age.MDM.pdf



# EDM
python2 scripts/make.Revigo.Input.py \
enrichr/GO_BP.Heart_Age.EDM.txt \
> Revigo/RevGO.Heart_Age.EDM.txt

# put to Revigo server

Rscript scripts/REVIGO_plot_read_data_2.r Revigo/Heart_Age.EDM.csv Revigo/Heart_Age.EDM.pdf
