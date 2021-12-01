import os

#Luca CHN_DICT_SR = {
#Luca     "Wen": [["vhcc_Wen_1_13TeV2017","Signal Region","#splitline{Resolved-jet}{1L (e)}"]],
#Luca     "Wmn": [["vhcc_Wmn_1_13TeV2017","Signal Region","#splitline{Resolved-jet}{1L (#mu)}"]],
#Luca     "Zee": [["vhcc_Zee_1_13TeV2017","Signal Region","#splitline{Resolved-jet}{2L (ee), High V-p_{T}}"],["vhcc_Zee_2_13TeV2017","Signal Region","#splitline{Resolved-jet}{2L (ee), Low V-p_{T}}"]],
#Luca     "Zmm": [["vhcc_Zmm_1_13TeV2017","Signal Region","#splitline{Resolved-jet}{2L (#mu#mu), High V-p_{T}}"],["vhcc_Zmm_2_13TeV2017","Signal Region","#splitline{Resolved-jet}{2L (#mu#mu), Low V-p_{T}}"]],
#Luca     "Znn": [["vhcc_Znn_1_13TeV2017","Signal Region","#splitline{Resolved-jet}{0L}"]]
#Luca }


#PRELIMINARY - FOR PAS
CHN_DICT_SR = {
    #"Wen": [["vhcc_Wen_1_13TeV2017","Signal Region","#splitline{#scale[1.2]{Preliminary}}{#splitline{Resolved-jet}{1L (e)}}"]],
    #"Wmn": [["vhcc_Wmn_1_13TeV2017","Signal Region","#splitline{#scale[1.2]{Preliminary}}{#splitline{Resolved-jet}{1L (#mu)}}"]],
    #"Zee": [["vhcc_Zee_1_13TeV2017","Signal Region","#splitline{#scale[1.2]{Preliminary}}{#splitline{Resolved-jet}{2L (ee), High V-p_{T}}}"],["vhcc_Zee_2_13TeV2017","Signal Region","#splitline{#scale[1.2]{Preliminary}}{#splitline{Resolved-jet}{2L (ee), Low V-p_{T}}}"]],
    #"Zmm": [["vhcc_Zmm_1_13TeV2017","Signal Region","#splitline{#scale[1.2]{Preliminary}}{#splitline{Resolved-jet}{2L (#mu#mu), High V-p_{T}}}"],["vhcc_Zmm_2_13TeV2017","Signal Region","#splitline{#scale[1.2]{Preliminary}}{#splitline{Resolved-jet}{2L (#mu#mu), Low V-p_{T}}}"]],
    #"Znn": [["vhcc_Znn_1_13TeV2017","Signal Region","#splitline{#scale[1.2]{Preliminary}}{#splitline{Resolved-jet}{0L}}"]]
	"Zetau": [["vhcc_Zetau_1_13TeV2017","Signal Region","#splitline{#scale[1.2]{Preliminary}}{Z(e#tau_{h})H(c#bar{c})}"]],
	"Zmutau": [["vhcc_Zmutau_1_13TeV2017","Signal Region","#splitline{#scale[1.2]{Preliminary}}{Z(#mu#tau_{h})H(c#bar{c})}"]],
}



#for MODE in ['prefit']:
#    for CHN in ['Zmutau','Zetau']:
#    for CHN in ['Zmutau']:
#    for CHN in ['Wen','Wmn','Zee','Zmm','Znn']:
#        for i in range(0,len(CHN_DICT_SR[CHN])):
#            LABEL = "%s" % CHN_DICT_SR[CHN][i][1]
#            OUTNAME = "%s" % CHN_DICT_SR[CHN][i][0]
#            EXTRALABEL = CHN_DICT_SR[CHN][i][2]
#            os.system(('./scripts/postFitPlot_vhcc_alt.py' \
#                       ' --file=shapes_postfit.root --ratio --extra_pad=0.53 --file_dir=%(OUTNAME)s --outdir="test_Zmutau1/" ' \
                       #' --ratio_range 0.0,200.0 --empty_bin_error --channel=%(CHN)s --blind --x_blind_min 0.0 --x_blind_max 200.0 --x_title "m_{jj}" --y_title "Events" \
                      # ' --ratio_range 0.0,2.0 --empty_bin_error --channel=%(CHN)s --x_title "BDT output" --y_title "Events" --mu ", #mu=1" --doVZ=True' \
#                       ' --outname %(OUTNAME)s --mode %(MODE)s --log_y --custom_y_range --y_axis_min "5E-3" --keepPreFitSignal True --lumi="41.5 fb^{-1} (13 TeV)"'\
#                       ' --outname %(OUTNAME)s --mode %(MODE)s --custom_y_range --y_axis_min "5E-3" --keepPreFitSignal True --lumi="41.5 fb^{-1} (13 TeV)"'\
#                       ' --channel_label "%(LABEL)s" --extralabel "%(EXTRALABEL)s"' % vars()))


for MODE in ['postfit']:
#    for CHN in ['Zmutau','Zetau']:
    for CHN in ['Zetau']:
#    for CHN in ['Wen','Wmn','Zee','Zmm','Znn']:
        for i in range(0,len(CHN_DICT_SR[CHN])):
              LABEL = "%s" % CHN_DICT_SR[CHN][i][1]
              OUTNAME = "%s" % CHN_DICT_SR[CHN][i][0]
              EXTRALABEL = CHN_DICT_SR[CHN][i][2]
              os.system(('./scripts/postFitPlot_vhcc_alt.py' \
                  ' --file=shapes_postfit.root --ratio --extra_pad=0.53 --file_dir=%(OUTNAME)s --outdir="test_post_sigx50_Zetau/"' \
                  ' --ratio_range 0.5,1.5 ' \
#                  ' --ratio_range 0.0,2.0 --empty_bin_error --channel=%(CHN)s --x_title "BDT output" --y_title "Events" --mu ", #mu=-9.5" ' \
#                  ' --ratio_range 0.0,2.0 --empty_bin_error --channel=%(CHN)s --x_title "BDT output" --y_title "Events" --mu ", #mu=0.96" --doVZ=True' \
#                  ' --outname %(OUTNAME)s --mode %(MODE)s --log_y --custom_y_range --y_axis_min "5E-3" --keepPreFitSignal True --lumi="41.5 fb^{-1} (13 TeV)"'\
                  ' --channel_label "%(LABEL)s" --extralabel "%(EXTRALABEL)s"' % vars()))


#option --ratio_justb

# shapes_VZincl_15May.root
#shapes_VZpt300_15May.root


