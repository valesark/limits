import CombineHarvester.CombineTools.ch as ch

def AddCommonSystematics(cb):
  
  signal = cb.cp().signals().process_set()
  # rateParams


########################################################################################################################################
### Uncertainties common to all Run-2 years
########################################################################################################################################
  
  # Theory uncertainties: signal
  cb.cp().AddSyst(cb,
                  'pdf_Higgs_qqbar', 'lnN', ch.SystMap('process')
                  (['ZH_hbb','ZH_hcc'],1.016)
                  (['WH_hbb','WH_hcc'],1.019))
  
  cb.cp().process(['ggZH_hbb','ggZH_hcc']).AddSyst(cb,'pdf_Higgs_gg', 'lnN', ch.SystMap()(1.024))

  cb.cp().process(['ggZH_hbb','ggZH_hcc']).AddSyst(cb,'QCDscale_ggZH', 'lnN',ch.SystMap()((1.251,0.811)))

  cb.cp().AddSyst(cb,'QCDscale_VH', 'lnN', ch.SystMap('process') 
                  (['ZH_hbb','ZH_hcc'], (1.038,0.969)) 
                  (['WH_hbb','WH_hcc'], (1.005,0.993)))

  
  cb.cp().process(['ZH_hcc','WH_hcc','ggZH_hcc']).AddSyst(cb,'BR_hcc', 'lnN', ch.SystMap()((1.05,0.97)))
  cb.cp().process(['ZH_hbb','WH_hbb','ggZH_hbb']).AddSyst(cb,'BR_hbb', 'lnN', ch.SystMap()(1.005))
  #Hbb branching ratio as measured by CMS 
  #cb.cp().process(['ZH_hbb','WH_hbb','ggZH_hbb']).AddSyst(cb,'BR_hbb', 'lnN', ch.SystMap()(1.20))
  
  #NLO EWK pt(V) correction to the VH and ggZH processes
  cb.cp().AddSyst(cb,
      'CMS_vhcc_boost_EWK_13TeV', 'lnN', ch.SystMap('channel','process') 
                  (['Zee','Zmm'],['ZH_hcc','ggZH_hcc','ZH_hbb','ggZH_hbb'], 1.02)
                  (['Znn'],['ZH_hcc','WH_hcc','ggZH_hcc','ZH_hbb','WH_hbb','ggZH_hbb'],1.02)
                  (['Wen','Wmn'],['WH_hcc','ZH_hcc','WH_hbb','ZH_hbb'],1.02)) 
  

  # Measured cross section uncertainties because we don't have SF
  cb.cp().process(['VVother','VWother','VZcc','VWcq']).AddSyst(cb,'CMS_vhcc_VV', 'lnN', ch.SystMap()(1.05)) 
  cb.cp().process(['s_Top']).AddSyst(cb,'CMS_vhcc_ST', 'lnN', ch.SystMap()(1.15)) 


  # Uncertainty on di-boson NNLO reweighting 
  cb.cp().process(['VVother','VWother','VZcc','VWcq']).AddSyst(cb, 'CMS_VV_NNLOWeights_13TeV', 'shape', ch.SystMap()(1.0))


  # Theoretical PDF uncertainties
  cb.cp().process(['ggZH_hbb','ggZH_hcc']).AddSyst(cb,'CMS_LHE_pdf_ggZH', 'lnN', ch.SystMap()(1.023)) 
  cb.cp().process(['ZH_hbb','ZH_hcc']).AddSyst(cb,'CMS_LHE_pdf_ZH', 'lnN', ch.SystMap()(1.018)) 
  cb.cp().process(['WH_hbb','WH_hcc']).AddSyst(cb,'CMS_LHE_pdf_WH', 'lnN', ch.SystMap()(1.018))
  cb.cp().process(['TT']).AddSyst(cb,'CMS_LHE_pdf_TT', 'lnN', ch.SystMap()(1.0265))
  cb.cp().process(['s_Top']).AddSyst(cb,'CMS_LHE_pdf_ST', 'lnN', ch.SystMap()(1.0288))
  cb.cp().process(['Zj_ll']).AddSyst(cb,'CMS_LHE_pdf_Zj_ll', 'lnN', ch.SystMap()(1.027))
  cb.cp().process(['Zj_bj']).AddSyst(cb,'CMS_LHE_pdf_Zj_bj', 'lnN', ch.SystMap()(1.027))
  cb.cp().process(['Zj_cj']).AddSyst(cb,'CMS_LHE_pdf_Zj_cj', 'lnN', ch.SystMap()(1.027))
  cb.cp().process(['Wj_ll']).AddSyst(cb,'CMS_LHE_pdf_Wj_ll', 'lnN', ch.SystMap()(1.027))
  cb.cp().process(['Wj_bj']).AddSyst(cb,'CMS_LHE_pdf_Wj_bj', 'lnN', ch.SystMap()(1.027))
  cb.cp().process(['Wj_cj']).AddSyst(cb,'CMS_LHE_pdf_Wj_cj', 'lnN', ch.SystMap()(1.027))
  cb.cp().process(['VZcc']).AddSyst(cb,'CMS_LHE_pdf_VZcc', 'lnN', ch.SystMap()(1.015))
  cb.cp().process(['VVother']).AddSyst(cb,'CMS_LHE_pdf_VVother', 'lnN', ch.SystMap()(1.0135)) 
  cb.cp().process(['VWother','VWcq']).AddSyst(cb,'CMS_LHE_pdf_VW', 'lnN', ch.SystMap()(1.0135)) 


  # Theoretical Renormalization and Factorization scale uncertainties
  cb.cp().process(['ZH_hbb','ZH_hcc']).AddSyst(cb,'CMS_LHE_weights_scale_muR_ZH','shape',ch.SystMap()(1.0))
  cb.cp().process(['WH_hbb','WH_hcc']).AddSyst(cb,'CMS_LHE_weights_scale_muR_WH','shape',ch.SystMap()(1.0))
  cb.cp().process(['ggZH_hbb','ggZH_hcc']).AddSyst(cb,'CMS_LHE_weights_scale_muR_ggZH','shape',ch.SystMap()(1.0))
  cb.cp().process(['ZH_hbb','ZH_hcc']).AddSyst(cb,'CMS_LHE_weights_scale_muF_ZH','shape',ch.SystMap()(1.0))
  cb.cp().process(['WH_hbb','WH_hcc']).AddSyst(cb,'CMS_LHE_weights_scale_muF_WH','shape',ch.SystMap()(1.0))
  cb.cp().process(['ggZH_hbb','ggZH_hcc']).AddSyst(cb,'CMS_LHE_weights_scale_muF_ggZH','shape',ch.SystMap()(1.0))
  cb.cp().process(['Zj_ll']).AddSyst(cb,'CMS_LHE_weights_scale_muR_Zj_ll','shape',ch.SystMap()(1.0))
  cb.cp().process(['Zj_ll']).AddSyst(cb,'CMS_LHE_weights_scale_muF_Zj_ll','shape',ch.SystMap()(1.0))
  cb.cp().process(['Zj_bj']).AddSyst(cb,'CMS_LHE_weights_scale_muR_Zj_bj','shape',ch.SystMap()(1.0))
  cb.cp().process(['Zj_bj']).AddSyst(cb,'CMS_LHE_weights_scale_muF_Zj_bj','shape',ch.SystMap()(1.0))
  cb.cp().process(['Zj_cj']).AddSyst(cb,'CMS_LHE_weights_scale_muR_Zj_cj','shape',ch.SystMap()(1.0))
  cb.cp().process(['Zj_cj']).AddSyst(cb,'CMS_LHE_weights_scale_muF_Zj_cj','shape',ch.SystMap()(1.0))
  cb.cp().process(['Wj_ll']).AddSyst(cb,'CMS_LHE_weights_scale_muR_Wj_ll','shape',ch.SystMap()(1.0))
  cb.cp().process(['Wj_ll']).AddSyst(cb,'CMS_LHE_weights_scale_muF_Wj_ll','shape',ch.SystMap()(1.0))
  cb.cp().process(['Wj_bj']).AddSyst(cb,'CMS_LHE_weights_scale_muR_Wj_bj','shape',ch.SystMap()(1.0))
  cb.cp().process(['Wj_bj']).AddSyst(cb,'CMS_LHE_weights_scale_muF_Wj_bj','shape',ch.SystMap()(1.0))
  cb.cp().process(['Wj_cj']).AddSyst(cb,'CMS_LHE_weights_scale_muR_Wj_cj','shape',ch.SystMap()(1.0))
  cb.cp().process(['Wj_cj']).AddSyst(cb,'CMS_LHE_weights_scale_muF_Wj_cj','shape',ch.SystMap()(1.0))
  cb.cp().process(['TT']).AddSyst(cb,'CMS_LHE_weights_scale_muR_TT','shape',ch.SystMap()(1.0))
  cb.cp().process(['TT']).AddSyst(cb,'CMS_LHE_weights_scale_muF_TT','shape',ch.SystMap()(1.0))
  cb.cp().process(['s_Top']).AddSyst(cb,'CMS_LHE_weights_scale_muR_ST','shape',ch.SystMap()(1.0))
  cb.cp().process(['s_Top']).AddSyst(cb,'CMS_LHE_weights_scale_muF_ST','shape',ch.SystMap()(1.0))
  cb.cp().process(['VVother']).AddSyst(cb,'CMS_LHE_weights_scale_muR_VVother','shape',ch.SystMap()(1.0))
  cb.cp().process(['VZcc']).AddSyst(cb,'CMS_LHE_weights_scale_muR_VZcc','shape',ch.SystMap()(1.0))
#  cb.cp().process(['VWother','VWcq']).AddSyst(cb,'CMS_LHE_weights_scale_muR_VW','shape',ch.SystMap()(1.0))
  cb.cp().process(['VVother']).AddSyst(cb,'CMS_LHE_weights_scale_muF_VVother','shape',ch.SystMap()(1.0))
  cb.cp().process(['VZcc']).AddSyst(cb,'CMS_LHE_weights_scale_muF_VZcc','shape',ch.SystMap()(1.0))
#  cb.cp().process(['VWother','VWcq']).AddSyst(cb,'CMS_LHE_weights_scale_muF_VW','shape',ch.SystMap()(1.0))


########################################################################################################################################
### Uncertainties for 2016
########################################################################################################################################
def AddSystematics2016(cb, splitJEC=False):
  #Additional lnN of 50%/100% on Wj_bj
  cb.cp().process(['Wj_bj']).AddSyst(cb,'Norm_Wj_bj_2016', 'lnN', ch.SystMap()(1.50))

####################### SCALE FACTORS RATEPARAM
  
  # TT Zll
  cb.cp().channel(['Zetau','Zmutau']).process(['TT']).AddSyst(cb,
     'SF_TT_high_ZtauX_2017', 'rateParam', ch.SystMap('bin_id')
     ([1],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['TT']).AddSyst(cb,
     'SF_TT_low_Zll_2016', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))
     
  cb.cp().channel(['Zee','Zmm']).process(['TT']).AddSyst(cb,
     'SF_TT_high_Zll_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Zj_ll Zll
  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll']).AddSyst(cb,
     'SF_Zj_ll_high_Zll_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll']).AddSyst(cb,
     'SF_Zj_ll_low_Zll_2016', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))

  # Zj_bj Zll
  cb.cp().channel(['Zee','Zmm']).process(['Zj_bj']).AddSyst(cb,
     'SF_Zj_bj_high_Zll_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['Zj_bj']).AddSyst(cb,
     'SF_Zj_bj_low_Zll_2016', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))

  # Zj_cj Zll
  cb.cp().channel(['Zee','Zmm']).process(['Zj_cj']).AddSyst(cb,
     'SF_Zj_cj_high_Zll_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['Zj_cj']).AddSyst(cb,
     'SF_Zj_cj_low_Zll_2016', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))

  # TT Znn
  cb.cp().channel(['Znn']).process(['TT']).AddSyst(cb,
     'SF_TT_Znn_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Zj_ll Znn
  cb.cp().channel(['Znn']).process(['Zj_ll']).AddSyst(cb,
     'SF_Zj_ll_Znn_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Zj_bj Znn
  cb.cp().channel(['Znn']).process(['Zj_bj']).AddSyst(cb,
     'SF_Zj_bj_Znn_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Zj_cj Znn
  cb.cp().channel(['Znn']).process(['Zj_cj']).AddSyst(cb,
     'SF_Zj_cj_Znn_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # TT Wln
  cb.cp().channel(['Wen','Wmn']).process(['TT']).AddSyst(cb,
     'SF_TT_Wln_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Wj_ll Wln
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_ll']).AddSyst(cb,
     'SF_Wj_ll_Wln_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Wj_cj Wln
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_cj','Wj_bj']).AddSyst(cb,
     'SF_Wj_cj_Wln_2016', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))


  #Set a sensible range for the rate params
  for syst in cb.cp().syst_type(["rateParam"]).syst_name_set():
    if not (syst=='SF_Wj_ll_Znn_2016' or syst=='SF_Wj_cj_Znn_2016'):
      cb.GetParameter(syst).set_range(0.0,5.0)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%% EXPERIMENTAL UNCERTAINTIES

  cb.cp().AddSyst(cb,'lumi_13TeV_2016','lnN', ch.SystMap()(1.010))
  cb.cp().AddSyst(cb,'lumi_13TeV_correlated','lnN', ch.SystMap()(1.006))
  cb.cp().AddSyst(cb,'CMS_vhcc_puWeight_2016','shape',ch.SystMap()(1.0))
#  cb.cp().AddSyst(cb,'CMS_vhcc_puBackupWeight_2016','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_PUIDWeight_13TeV_2016','shape',ch.SystMap()(1.0))
 
#============= Prefire efficiencies
  cb.cp().channel(['Zee','Zmm','Wmn','Wen']).AddSyst(cb,'CMS_PrefireWeight','shape',ch.SystMap()(1.0))



#============= lepton efficiencies

  cb.cp().channel(['Wmn']).AddSyst(cb,'CMS_vhcc_eff_m_Wln_13TeV_2016','lnN',ch.SystMap()(1.02))
  cb.cp().channel(['Wen']).AddSyst(cb,'CMS_vhcc_eff_e_Wln_13TeV_2016','lnN',ch.SystMap()(1.02))
  cb.cp().channel(['Zmm']).AddSyst(cb,'CMS_vhcc_eff_m_Zll_13TeV_2016','lnN',ch.SystMap()(1.04))
  cb.cp().channel(['Zee']).AddSyst(cb,'CMS_vhcc_eff_e_Zll_13TeV_2016','lnN',ch.SystMap()(1.04))


#=============  met efficiencies
  cb.cp().channel(['Znn']).AddSyst(cb,'CMS_vhcc_trigger_MET_13TeV_2016','lnN',ch.SystMap()(1.02))

#=============  VpT reweightings
  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_ptzweights_13TeV_2016','shape',ch.SystMap()(1.0))


#============= Jet energy scale and resolution

#Luca  cb.cp().AddSyst(cb,'CMS_res_j_13TeV_2016','shape',ch.SystMap()(1.0))
#Luca  cb.cp().channel(['Wen','Wmn','Znn']).AddSyst(cb,'CMS_METUnclustEn','shape',ch.SystMap()(1.0)) 

#Luca  cb.cp().process(['ZH_hcc','WH_hcc','ggZH_hcc','ZH_hbb','WH_hbb','ggZH_hbb','VVother','VZcc','VWother','VWcq']).AddSyst(cb,'CMS_j_PtCReg_Scale','shape',ch.SystMap()(1.0))
#Luca  cb.cp().process(['ZH_hcc','WH_hcc','ggZH_hcc','ZH_hbb','WH_hbb','ggZH_hbb','VVother','VZcc','VWother','VWcq']).AddSyst(cb,'CMS_j_PtCReg_Smear','shape',ch.SystMap()(1.0))

#Luca  if splitJEC:
#Luca    #split as JET/MET recommends - NEW 11-splitting scheme
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_Absolute','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_BBEC1','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_EC2','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_FlavorQCD','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_HF','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_RelativeBal','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_Absolute_2016','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_BBEC1_2016','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_EC2_2016','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_HF_2016','shape',ch.SystMap()(1.0))
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_RelativeSample_2016','shape',ch.SystMap()(1.0))
#Luca     
#Luca  else:
#Luca    cb.cp().AddSyst(cb,'CMS_scale_j_13TeV_2016','shape',ch.SystMap()(1.0))

#============= c-tagger uncertainties - inclusive in pt/eta

  cb.cp().AddSyst(cb,'CMS_cTagWeight_JER','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_JES','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_PU','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_muR','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_muF','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_Stat_2016','shape',ch.SystMap()(1.0))    
  cb.cp().AddSyst(cb,'CMS_cTagWeight_EleId','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_MuId','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecDYJets','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecST','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecWJets','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecTTbar','shape',ch.SystMap()(1.0))

#============= post-cTagWeight uncertainties

  cb.cp().AddSyst(cb,'CMS_PostCTagWeight_13TeV_2016','shape',ch.SystMap()(1.0))

#fit uncertainties on DY(ll)+jets process (taken from 2L, but applied to DY+jets in all the channels)
  cb.cp().channel(['Zee','Zmm','Wen','Wmn','Znn']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Zll_2016','shape',ch.SystMap()(1.0))

#fit uncertainties on W+jets process (taken from 1L, but applied to W+jets in all the channels: 1L and 0L)
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_ll','Wj_bj','Wj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Wln_2016','shape',ch.SystMap()(1.0))

#fit uncertainties on Z(nn)+jets process (taken from 0L, but applied to Z(nn)+jets in all the channels: 0L only in this case)
  cb.cp().channel(['Znn']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Znn_2016','shape',ch.SystMap()(1.0))
  
  cb.cp().channel(['Zee','Zmm']).process(['Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_cj_2016','shape',ch.SystMap()(1.0))
  cb.cp().channel(['Zee','Zmm']).process(['Zj_bj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_bj_2016','shape',ch.SystMap()(1.0))  
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_cj_2016','shape',ch.SystMap()(1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Zj_bj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_bj_2016','shape',ch.SystMap()(1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Wj_cj_2016','shape',ch.SystMap()(1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_bj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Wj_bj_2016','shape',ch.SystMap()(1.0))




########################################################################################################################################
### Uncertainties for 2017
########################################################################################################################################
def AddSystematics2017(cb, splitJEC=False):

  cb.cp().process(['Wj_bj']).AddSyst(cb,'Norm_Wj_bj_2017', 'lnN', ch.SystMap()(1.50))

####################### SCALE FACTORS RATEPARAM

#=================================================================
#=================================================================
  
   # TT Zll
  cb.cp().channel(['Zee','Zmm']).process(['TT']).AddSyst(cb,
     'SF_TT_high_Zll_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
  cb.cp().channel(['Zee','Zmm']).process(['TT']).AddSyst(cb,
     'SF_TT_low_Zll_2017', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))
 
  # Zj_ll Zll
  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll']).AddSyst(cb,
     'SF_Zj_ll_high_Zll_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll']).AddSyst(cb,
     'SF_Zj_ll_low_Zll_2017', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))
 
  # Zj_bj Zll
  cb.cp().channel(['Zee','Zmm']).process(['Zj_bj']).AddSyst(cb,
     'SF_Zj_bj_high_Zll_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
  cb.cp().channel(['Zee','Zmm']).process(['Zj_bj']).AddSyst(cb,
     'SF_Zj_bj_low_Zll_2017', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))
 
  # Zj_cj Zll
  cb.cp().channel(['Zee','Zmm']).process(['Zj_cj']).AddSyst(cb,
     'SF_Zj_cj_high_Zll_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
  cb.cp().channel(['Zee','Zmm']).process(['Zj_cj']).AddSyst(cb,
     'SF_Zj_cj_low_Zll_2017', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))
 
  # TT Znn
  cb.cp().channel(['Znn']).process(['TT']).AddSyst(cb,
     'SF_TT_Znn_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
  # Zj_ll Znn
  cb.cp().channel(['Znn']).process(['Zj_ll']).AddSyst(cb,
     'SF_Zj_ll_Znn_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
  # Zj_bj Znn
  cb.cp().channel(['Znn']).process(['Zj_bj']).AddSyst(cb,
     'SF_Zj_bj_Znn_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
  # Zj_cj Znn
  cb.cp().channel(['Znn']).process(['Zj_cj']).AddSyst(cb,
     'SF_Zj_cj_Znn_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
 
  # TT Wln
  cb.cp().channel(['Wen','Wmn']).process(['TT']).AddSyst(cb,
     'SF_TT_Wln_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
  # Wj_ll Wln
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_ll']).AddSyst(cb,
     'SF_Wj_ll_Wln_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))
 
 #BASELINE-RESTORE  # Wj_bj Wln
 #BASELINE-RESTORE  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_bj']).AddSyst(cb,
 #BASELINE-RESTORE     'SF_Wj_bj_Wln_2017', 'rateParam', ch.SystMap('bin_id')
 #BASELINE-RESTORE     ([1,3,5,7,9],1.0))
 
  # Wj_cj Wln
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_cj','Wj_bj']).AddSyst(cb,
     'SF_Wj_cj_Wln_2017', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))


  #Set a sensible range for the rate params
  for syst in cb.cp().syst_type(["rateParam"]).syst_name_set():
    if not (syst=='SF_Wj_ll_Znn_2017' or syst=='SF_Wj_cj_Znn_2017'):
      cb.GetParameter(syst).set_range(0.0,5.0)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%% EXPERIMENTAL UNCERTAINTIES
  cb.cp().AddSyst(cb,'lumi_13TeV_2017','lnN', ch.SystMap()(1.020))
  cb.cp().AddSyst(cb,'lumi_13TeV_1718','lnN', ch.SystMap()(1.006))
  cb.cp().AddSyst(cb,'lumi_13TeV_correlated','lnN', ch.SystMap()(1.009))  
  cb.cp().AddSyst(cb,'CMS_vhcc_puWeight_2017','shape',ch.SystMap()(1.0))
##  cb.cp().AddSyst(cb,'CMS_vhcc_puBackupWeight_2017','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_PUIDWeight_13TeV_2017','shape',ch.SystMap()(1.0))

#============= Prefire efficiencies
  cb.cp().channel(['Zee','Zmm']).AddSyst(cb,'CMS_PrefireWeight','shape',ch.SystMap()(1.0))

#============= lepton efficiencies

  cb.cp().channel(['Wmn']).AddSyst(cb,'CMS_vhcc_eff_m_Wln_13TeV_2017','lnN',ch.SystMap()(1.02))
  cb.cp().channel(['Wen']).AddSyst(cb,'CMS_vhcc_eff_e_Wln_13TeV_2017','lnN',ch.SystMap()(1.02))
  cb.cp().channel(['Zmm']).AddSyst(cb,'CMS_vhcc_eff_m_Zll_13TeV_2017','lnN',ch.SystMap()(1.04))
  cb.cp().channel(['Zee']).AddSyst(cb,'CMS_vhcc_eff_e_Zll_13TeV_2017','lnN',ch.SystMap()(1.04))

#=============  met efficiencies
  cb.cp().channel(['Znn']).AddSyst(cb,'CMS_vhcc_trigger_MET_13TeV_2017','lnN',ch.SystMap()(1.02))

#=============  VpT reweightings - to comment for NLO
  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_ptzweights_13TeV_2017','shape',ch.SystMap()(1.0))

  
#============= Jet energy scale and resolution
  #cb.cp().AddSyst(cb,'CMS_res_j_13TeV_2017','shape',ch.SystMap()(1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).AddSyst(cb,'CMS_METUnclustEn','shape',ch.SystMap()(1.0)) 

  #cb.cp().process(['ZH_hcc','WH_hcc','ggZH_hcc','ZH_hbb','WH_hbb','ggZH_hbb','VVother','VZcc','VWother','VWcq']).AddSyst(cb,'CMS_j_PtCReg_Scale','shape',ch.SystMap()(1.0))
  #cb.cp().process(['ZH_hcc','WH_hcc','ggZH_hcc','ZH_hbb','WH_hbb','ggZH_hbb','VVother','VZcc','VWother','VWcq']).AddSyst(cb,'CMS_j_PtCReg_Smear','shape',ch.SystMap()(1.0))


  if splitJEC:
    #split as JET/MET recommends - NEW 11-splitting scheme
    cb.cp().AddSyst(cb,'CMS_scale_j_Absolute','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_BBEC1','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_EC2','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_FlavorQCD','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_HF','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_RelativeBal','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_Absolute_2017','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_BBEC1_2017','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_EC2_2017','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_HF_2017','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_RelativeSample_2017','shape',ch.SystMap()(1.0))

    
  #else:
    #cb.cp().AddSyst(cb,'CMS_scale_j_13TeV_2017','shape',ch.SystMap()(1.0))

#============= tagger uncertainties

  #cb.cp().AddSyst(cb,'CMS_cTagWeight_JER','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_JES','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_PU','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_EleId','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_MuId','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_muR','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_muF','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecDYJets','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecST','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecWJets','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecTTbar','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_Stat_2017','shape',ch.SystMap()(1.0))

#============= post-cTagWeight uncertainties

  #cb.cp().AddSyst(cb,'CMS_PostCTagWeight_13TeV_2017','shape',ch.SystMap()(1.0))


#============= dRjj-Reweighting uncertainties

#fit uncertainties on DY(ll)+jets process (taken from 2L, but applied to DY+jets in all the channels)
#  cb.cp().channel(['Zee','Zmm','Wen','Wmn','Znn']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Zll_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6],1.0))
  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Zll_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))

#fit uncertainties on W+jets process (taken from 1L, but applied to W+jets in all the channels: 1L and 0L)
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_ll','Wj_bj','Wj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Wln_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))

#fit uncertainties on Z(nn)+jets process (taken from 0L, but applied to Z(nn)+jets in all the channels: 0L only in this case)
  cb.cp().channel(['Znn']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Znn_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_cj_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))
  cb.cp().channel(['Zee','Zmm']).process(['Zj_bj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_bj_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0)) 
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_cj_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Zj_bj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_bj_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Wj_cj_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_bj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Wj_bj_2017','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))



########################################################################################################################################
### Uncertainties for 2018
########################################################################################################################################
def AddSystematics2018(cb, splitJEC=False):

  cb.cp().process(['Wj_bj']).AddSyst(cb,'Norm_Wj_bj_2018', 'lnN', ch.SystMap()(1.50))

####################### SCALE FACTORS RATEPARAM
  
  # TT Zll
  cb.cp().channel(['Zee','Zmm']).process(['TT']).AddSyst(cb,
     'SF_TT_high_Zll_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['TT']).AddSyst(cb,
     'SF_TT_low_Zll_2018', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))

  # Zj_ll Zll
  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll']).AddSyst(cb,
     'SF_Zj_ll_high_Zll_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll']).AddSyst(cb,
     'SF_Zj_ll_low_Zll_2018', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))

  # Zj_bj Zll
  cb.cp().channel(['Zee','Zmm']).process(['Zj_bj']).AddSyst(cb,
     'SF_Zj_bj_high_Zll_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['Zj_bj']).AddSyst(cb,
     'SF_Zj_bj_low_Zll_2018', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))

  # Zj_cj Zll
  cb.cp().channel(['Zee','Zmm']).process(['Zj_cj']).AddSyst(cb,
     'SF_Zj_cj_high_Zll_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['Zj_cj']).AddSyst(cb,
     'SF_Zj_cj_low_Zll_2018', 'rateParam', ch.SystMap('bin_id')
     ([2,4,6,8,10],1.0))

  # TT Znn
  cb.cp().channel(['Znn']).process(['TT']).AddSyst(cb,
     'SF_TT_Znn_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Zj_ll Znn
  cb.cp().channel(['Znn']).process(['Zj_ll']).AddSyst(cb,
     'SF_Zj_ll_Znn_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Zj_bj Znn
  cb.cp().channel(['Znn']).process(['Zj_bj']).AddSyst(cb,
     'SF_Zj_bj_Znn_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Zj_cj Znn
  cb.cp().channel(['Znn']).process(['Zj_cj']).AddSyst(cb,
     'SF_Zj_cj_Znn_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # TT Wln
  cb.cp().channel(['Wen','Wmn']).process(['TT']).AddSyst(cb,
     'SF_TT_Wln_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

  # Wj_ll Wln
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_ll']).AddSyst(cb,
     'SF_Wj_ll_Wln_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))

#BASELINE-RESTORE  # Wj_bj Wln
#BASELINE-RESTORE  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_bj']).AddSyst(cb,
#BASELINE-RESTORE     'SF_Wj_bj_Wln_2018', 'rateParam', ch.SystMap('bin_id')
#BASELINE-RESTORE     ([1,3,5,7,9],1.0))

  # Wj_cj Wln
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_cj','Wj_bj']).AddSyst(cb,
     'SF_Wj_cj_Wln_2018', 'rateParam', ch.SystMap('bin_id')
     ([1,3,5,7,9],1.0))


  #Set a sensible range for the rate params
  for syst in cb.cp().syst_type(["rateParam"]).syst_name_set():
    if not (syst=='SF_Wj_ll_Znn_2018' or syst=='SF_Wj_cj_Znn_2018'):
      cb.GetParameter(syst).set_range(0.0,5.0)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%% EXPERIMENTAL UNCERTAINTIES
  cb.cp().AddSyst(cb,'lumi_13TeV_2018','lnN', ch.SystMap()(1.015))
  cb.cp().AddSyst(cb,'lumi_13TeV_1718','lnN', ch.SystMap()(1.002))
  cb.cp().AddSyst(cb,'lumi_13TeV_correlated','lnN', ch.SystMap()(1.020))
  cb.cp().AddSyst(cb,'CMS_vhcc_puWeight_2018','shape',ch.SystMap()(1.0))
#  cb.cp().AddSyst(cb,'CMS_vhcc_puBackupWeight_2018','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_PUIDWeight_13TeV_2018','shape',ch.SystMap()(1.0))

#============= HEM 15/16 uncertainties

  cb.cp().AddSyst(cb,'CMS_scale_j_HEMIssue_13TeV_2018','shape',ch.SystMap()(1.0))

#============= lepton efficiencies

  cb.cp().channel(['Wmn']).AddSyst(cb,'CMS_vhcc_eff_m_Wln_13TeV_2018','lnN',ch.SystMap()(1.02))
  cb.cp().channel(['Wen']).AddSyst(cb,'CMS_vhcc_eff_e_Wln_13TeV_2018','lnN',ch.SystMap()(1.02))
  cb.cp().channel(['Zmm']).AddSyst(cb,'CMS_vhcc_eff_m_Zll_13TeV_2018','lnN',ch.SystMap()(1.04))
  cb.cp().channel(['Zee']).AddSyst(cb,'CMS_vhcc_eff_e_Zll_13TeV_2018','lnN',ch.SystMap()(1.04))

#=============  met efficiencies
  cb.cp().channel(['Znn']).AddSyst(cb,'CMS_vhcc_trigger_MET_13TeV_2018','lnN',ch.SystMap()(1.02))

#=============  VpT reweightings - to comment for NLO
  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_ptzweights_13TeV_2018','shape',ch.SystMap()(1.0))


#============= Jet energy scale and resolution
  cb.cp().AddSyst(cb,'CMS_res_j_13TeV_2018','shape',ch.SystMap()(1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).AddSyst(cb,'CMS_METUnclustEn','shape',ch.SystMap()(1.0)) 

  cb.cp().process(['ZH_hcc','WH_hcc','ggZH_hcc','ZH_hbb','WH_hbb','ggZH_hbb','VVother','VZcc','VWother','VWcq']).AddSyst(cb,'CMS_j_PtCReg_Scale','shape',ch.SystMap()(1.0))
  cb.cp().process(['ZH_hcc','WH_hcc','ggZH_hcc','ZH_hbb','WH_hbb','ggZH_hbb','VVother','VZcc','VWother','VWcq']).AddSyst(cb,'CMS_j_PtCReg_Smear','shape',ch.SystMap()(1.0))

  if splitJEC:
    #split as JET/MET recommends - NEW 11-splitting scheme
    cb.cp().AddSyst(cb,'CMS_scale_j_Absolute','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_BBEC1','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_EC2','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_FlavorQCD','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_HF','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_RelativeBal','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_Absolute_2018','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_BBEC1_2018','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_EC2_2018','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_HF_2018','shape',ch.SystMap()(1.0))
    cb.cp().AddSyst(cb,'CMS_scale_j_RelativeSample_2018','shape',ch.SystMap()(1.0))
    
  else:
    cb.cp().AddSyst(cb,'CMS_scale_j_13TeV_2018','shape',ch.SystMap()(1.0))

#============= tagger uncertainties
  cb.cp().AddSyst(cb,'CMS_cTagWeight_JER','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_JES','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_PU','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_EleId','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_MuId','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_muR','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_muF','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecDYJets','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecST','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecWJets','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_XSecTTbar','shape',ch.SystMap()(1.0))
  cb.cp().AddSyst(cb,'CMS_cTagWeight_Stat_2018','shape',ch.SystMap()(1.0))  


#============= post-cTagWeight uncertainties

  cb.cp().AddSyst(cb,'CMS_PostCTagWeight_13TeV_2018','shape',ch.SystMap()(1.0))


#============= dRjj-Reweighting uncertainties

#fit uncertainties on DY(ll)+jets process (taken from 2L, but applied to DY+jets in all the channels)
#  cb.cp().channel(['Zee','Zmm','Wen','Wmn','Znn']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Zll_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6],1.0))
  cb.cp().channel(['Zee','Zmm']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Zll_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))

#fit uncertainties on W+jets process (taken from 1L, but applied to W+jets in all the channels: 1L and 0L)
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_ll','Wj_bj','Wj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Wln_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))

#fit uncertainties on Z(nn)+jets process (taken from 0L, but applied to Z(nn)+jets in all the channels: 0L only in this case)
  cb.cp().channel(['Znn']).process(['Zj_ll','Zj_bj','Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_fit_Znn_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))

  cb.cp().channel(['Zee','Zmm']).process(['Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_cj_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))
  cb.cp().channel(['Zee','Zmm']).process(['Zj_bj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_bj_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0)) 
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Zj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_cj_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Zj_bj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Zj_bj_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_cj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Wj_cj_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))
  cb.cp().channel(['Wen','Wmn','Znn']).process(['Wj_bj']).AddSyst(cb,'CMS_vhcc_dRjjReweight_Flavour_Wj_bj_2018','shape',ch.SystMap('bin_id')([1,2,3,4,5,6,9,10],1.0))

