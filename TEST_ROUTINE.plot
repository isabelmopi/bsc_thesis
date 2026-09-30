#global settings for all plots
BEGIN PLOT /TEST_ROUTINE/.*
LogY=0
NormalizeToSum=1
YLabel=Normalized
RatioPlotYMin=0.73
RatioPlotYMax=1.27
LegendXPos=0.71
LegendYPos=0.85
LegendSize=0.035
LineColor=black
LineWidth=0.9
PlotTopMargin=0.155
PlotResolution=600
RatioPlot=0
END PLOT

#local settings for separate plots
BEGIN PLOT /TEST_ROUTINE/smallRjets_n
Title=Small-R jet multiplicity (at least 3 per event)
XLabel=Number of jets
END PLOT

BEGIN PLOT /TEST_ROUTINE/smallRjets_pT_all
Title=All small-R jets transverse momenta
XLabel=$p_{T}^{\mathrm{smallRjets}}$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/smallRjets_pT_lead
Title=Leading small-R jets transverse momenta
XLabel=$p_{T}$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/bT_smallRjets_n_precut
Title=B-tagged small-R jet multiplicity (precut)
XLabel=Number of jets
END PLOT

BEGIN PLOT /TEST_ROUTINE/bT_smallRjets_n
Title=B-tagged small-R jet multiplicity (at least 2 per event)
XLabel=Number of jets
END PLOT

BEGIN PLOT /TEST_ROUTINE/bT_smallRjets_pT_all
Title=B-tagged small-R jet transverse momenta
XLabel=$p_{T}$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/largeRjets_n
Title=Large-R jet multiplicity (after first selection cut)
XLabel=Number of jets [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/largeRjets_mass_all
Title=All large-R jet masses
XLabel=$m$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/lep_pT
Title=Prompt lepton transverse momenta
XLabel=$p_{T}$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/elec_pT
Title=Prompt electron transverse momenta
XLabel=$p_{T}$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/muon_pT
Title=Prompt muon transverse momenta
XLabel=$p_{T}$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/hadrTop_smallRjets_n
Title=Small-R jet constituent multiplicity (inside large-R, 2 or 3 per event)
XLabel=Number of jets
END PLOT

BEGIN PLOT /TEST_ROUTINE/hadrTopJet_mass
Title=Reconstructed hadronic top mass
XLabel=$m_t^{\text{had}}$ [GeV]
RatioPlotYMin=0.9
RatioPlotYMax=1.15
LegendXPos=0.05
END PLOT

BEGIN PLOT /TEST_ROUTINE/hadrTopJet_pT
Title=Reconstructed large-R (hadr.) Top-jet transverse momenta
XLabel=$p_{T}$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/sumLep_bTjet_mass
Title=Reconstructed leptonic Top mass (without neutrino)
XLabel=$m$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/sumLep_bTjet_pT
Title=Reconstructed leptonic Top pT (without neutrino)
XLabel=$p_{T}$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/lepTop_mass
Title=Reconstructed leptonic top mass
RatioPlotYMin=0.92
RatioPlotYMax=1.04
LegendXPos=0.52
XLabel=$m_t^\ell$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/lepTop_pT
Title=Reconstructed leptonic Top pT
XLabel=$p_{T}$ [GeV]
END PLOT

BEGIN PLOT /TEST_ROUTINE/ttbar_mass
Title=Reconstructed mass of $t\bar{t}$ system
RatioPlotYMin=0.95
RatioPlotYMax=1.08
LegendXPos=0.5
XLabel=$m_{t\bar{t}}$ [GeV]
END PLOT


#Polarisation plots.
BEGIN PLOT /TEST_ROUTINE/cosTheta_k
Title=Top polarization along $k$-axis 
XLabel=$\cos\theta_k$
LegendXPos=0.5
LegendYPos=0.3
END PLOT

BEGIN PLOT /TEST_ROUTINE/cosTheta_r
Title=Top polarization along $r$-axis
XLabel=$\cos\theta_r$
LegendXPos=0.5
LegendYPos=0.3
END PLOT

BEGIN PLOT /TEST_ROUTINE/cosTheta_n
Title=Top polarization along $n$-axis 
XLabel=$\cos\theta_n$
LegendXPos=0.5
LegendYPos=0.3
END PLOT

BEGIN PLOT /TEST_ROUTINE/cosThetaBar_k
Title=Antitop polarization along $\bar{k}$-axis 
XLabel=$\cos\bar{\theta}_k$
LegendXPos=0.5
LegendYPos=0.3
END PLOT

BEGIN PLOT /TEST_ROUTINE/cosThetaBar_r
Title=Antitop polarization along $\bar{r}$-axis
XLabel=$\cos\bar{\theta}_r$
LegendXPos=0.5
LegendYPos=0.3
END PLOT

BEGIN PLOT /TEST_ROUTINE/cosThetaBar_n
Title=Antitop polarization along $\bar{n}$-axis 
XLabel=$\cos\bar{\theta}_n$
LegendXPos=0.5
LegendYPos=0.3
END PLOT

#Correlation plots.
BEGIN PLOT /TEST_ROUTINE/C_kk
Title=Angular distribution for $C(k,k)$
XLabel=$\cos\theta_k \cdot \cos\bar{\theta}_k$
LegendXPos=0.06
LegendYPos=0.5
RatioPlotYMin=0.98
RatioPlotYMax=1.03
END PLOT

BEGIN PLOT /TEST_ROUTINE/C_rr
Title=Angular distribution for $C(r,r)$
XLabel=$\cos\theta_r \cdot \cos\bar{\theta}_r$
RatioPlotYMin=0.98
RatioPlotYMax=1.02
LegendXPos=0.051
LegendYPos=0.6
END PLOT

BEGIN PLOT /TEST_ROUTINE/C_nn
Title=Angular distribution for $C(n,n)$
XLabel=$\cos\theta_n \cdot \cos\bar{\theta}_n$
RatioPlotYMin=0.99
RatioPlotYMax=1.015
LegendYPos=0.25
LegendXPos=0.05
END PLOT

#Non-diagonals.
BEGIN PLOT /TEST_ROUTINE/C_nk_plus
Title=$C(n,k)+C(k,n)$
XLabel=$\cos\theta_n\cos\bar{\theta}_k+\cos\theta_k\cos\bar{\theta}_n$
END PLOT

BEGIN PLOT /TEST_ROUTINE/C_nk_minus
Title=$C(n,k)-C(k,n)$
XLabel=$\cos\theta_n\cos\bar{\theta}_k-\cos\theta_k\cos\bar{\theta}_n$
END PLOT

BEGIN PLOT /TEST_ROUTINE/C_nr_plus
Title=$C(n,r)+C(r,n)$
XLabel=$\cos\theta_n\cos\bar{\theta}_r+\cos\theta_r\cos\bar{\theta}_n$
END PLOT

BEGIN PLOT /TEST_ROUTINE/C_nr_minus
Title=$C(n,r)-C(r,n)$
XLabel=$\cos\theta_n\cos\bar{\theta}_r-\cos\theta_r\cos\bar{\theta}_n$
END PLOT

BEGIN PLOT /TEST_ROUTINE/C_rk_plus
Title=$C(r,k)+C(k,r)$
XLabel=$\cos\theta_r\cos\bar{\theta}_k+\cos\theta_k\cos\bar{\theta}_r$
END PLOT

BEGIN PLOT /TEST_ROUTINE/C_rk_minus
Title=$C(r,k)+C(k,r)$
XLabel=$\cos\theta_r\cos\bar{\theta}_k-\cos\theta_k\cos\bar{\theta}_r$
END PLOT

