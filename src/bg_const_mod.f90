module bg_const_mod
!----------------------------------------------------------------------------
!  Purpose: biological and biogeochemical constants
!           Notice: normally [DOC] is measured on the basis of carbon, e.g. 
!           mg C L-1 or umol C L-1, not the total mass of DOM if not clearly 
!           indicated!!
!           All photochemical parameters are from "Seasonality of 
!           photochemical dissolved organic carbon mineralization and its 
!           relative contribution to pelagic CO2 production in northern 
!           lakes" (Vachon et al., 2016).
!----------------------------------------------------------------------------
   use shr_kind_mod,    only : r8
   use shr_ctrl_mod,    only : NLAKTYPE  ! NPOC
   public
   ! lake type identifier
   !integer, parameter :: temperate_lake = 001   
   !--------------------------------------------------------------
   !
   ! Vegetaion type has proved to be one of the most important factor for 
   ! thermokarst lakes carbon budget on QTP. Here Lin has changed the 
   ! lake type identifier from climate identification to vegetation.
   !
   !--------------------------------------------------------------
   !integer, parameter :: small_ppk = 001, large_ppk = 002  ! Lin has canceled the distinguish between different types of phytoplankton
   ! integer, parameter :: pasC = 001, actC = 002  ! closed by Lin, so as follows
   real(r8), parameter :: fwlnd(NLAKTYPE) = (/1.0,0.75,0.375,0.125/)     ! added by Lin, 20250218
   ! freshwater lake pH
   real(r8), parameter :: SMpH(NLAKTYPE) = (/8.1,8.3,8.0,7.8/)  ! added by Lin for July&August; according to Mu et al., 2022
   real(r8), parameter :: LKpH(NLAKTYPE) = (/7.4,7.3,7.2,7.1/)  ! changed by Lin, according to Jia,2021;Mu et al., 2022
   real(r8), parameter :: SedpH(NLAKTYPE) = (/7.95,8.31,8.44,8.53/)! changed by Lin, according to Chi, 2023
   ! CH4 oxidation Reference temperature (K)
   real(r8), parameter :: Tor(2) = (/278.15, 283.15/)  ! for anaerobic and aerobic, according to Feng,2018;Gu,2021;Xu et al.,2024
   ! CH4 production reference temperature for passive matter (K)
   real(r8), parameter :: Tpr = 278.15
   ! CH4 production reference temperature for active matter
   ! real(r8), parameter :: Tpr_act = 273.15
   ! O2, CH4 and CO2 minimum dissolved concentrations (umol/m3)
   real(r8), parameter :: minDo2 = 1.0d+2
   real(r8), parameter :: minDch4 = 1.0d+1
   real(r8), parameter :: minDco2 = 1.0d+1
   ! the initial density of 14C-depleted organic matter (kg/m3)
   ! real(r8), parameter :: oldCarb0 = 29.3
   ! the suppression coefficient of O2 to methanogenesis (m3 water mol-1)
   real(r8), parameter :: etaO2 = 400.0
   ! minimum, maximum and optimum pH for CH4 production
   real(r8), parameter :: PpHmin = 6.0, PpHmax = 10.0, PpHopt = 7.0
   ! minimum, maximum and optimum pH for CH4 oxidation
   real(r8), parameter :: OpHmin = 5.0, OpHmax = 8.0, OpHopt = 6.0
   ! The carbon-specific absorption coefficient by NAP at 440 nm (m2 gC-1)
   real(r8), parameter :: aNAP_440 = 0.1
   ! The chlorophill specific coefficient of phytoplankton (m2 mg-1)
   !real(r8), parameter :: apico_676 = 0.023
   real(r8), parameter :: amicro_674 = 0.0086
   ! The diameter (m) and density (kg/m3) of detritus
   real(r8), parameter :: daDetrs = 8.0d-5
   real(r8), parameter :: dsDetrs = 1.04d+3
   ! The diameter (m) of small and large phytoplankton
   !real(r8), parameter :: daPico = 3.0d-6
   real(r8), parameter :: daMicro = 1.0d-5
   ! The background backscattering coefficient (m-1)
   real(r8), parameter :: Bbbg = 1.7d-4
   ! The backscattering ratio of sea water
   real(r8), parameter :: Bbsw = 0.5_r8
   ! Proportion of biosynthate allocated to synthesis of biomass of the
   ! biosynthetic machinery (Geider et al., 1996)
   real(r8), parameter :: kE = 0.6_r8
   ! The minimum and maximum carbon to chlorophill ratio (mmol C mg chl-1)
   ! Wang et al. (2009)
   real(r8), parameter :: C2Chlmin = 12.0 !C2Chlmin(NPOC) = (/12.0, 12.0/)
   real(r8), parameter :: C2Chlmax = 30.0 !C2Chlmax(NPOC) = (/40.0, 30.0/)
   ! The slope of C:Chl ratio vs. growth rate (mg C mg chl-1 d)
   real(r8), parameter :: Kpc2chl = 70 !Kpc2chl(NPOC) = (/95, 70/)
   ! Apparant quantum yield of photochemical degradation
   ! units: AQYem1 (mol C mol photons-1), AQYm2 (nm-1)
   real(r8), parameter :: AQYem1(NLAKTYPE) = (/2.6d-4, 2.6d-4, 2.6d-4, 2.6d-4/)! Lin didn't distinguish the value among different vegetation
   real(r8), parameter :: AQYm2(NLAKTYPE) = (/0.018, 0.018, 0.018, 0.018/)
   ! The ratio of partial photo-oxidation to photo-mineralization
   !real(r8), parameter :: PPO2PM(NLAKTYPE) = (/2.547/)
   ! The spectral slope of CDOM absorption (nm-1)
   real(r8), parameter :: SCDOM(NLAKTYPE) = (/0.018, 0.018, 0.018, 0.018/)
   ! The reference CDOM-specific UV absorption at 305 nm (m3 g C-1 m-1)
   real(r8), parameter :: SUVA305(NLAKTYPE) = (/2.0, 2.0, 2.0, 2.0/)
   ! the maximum growth rate of phytoplankton at 0 celcius (d-1) 
   real(r8), parameter :: mu0 = 1.0 !mu0(NPOC) = (/0.4, 1.0/) ! changed from "Vm0" to "mu0" according to ALBM-update
   ! temperature multiplier for phytoplankton growth
   real(r8), parameter :: ThetaG = 1.08
   ! temperature function parameters for phytoplankton growth
   real(r8), parameter :: kt_ppk = 12.76830  !kt_ppk(NPOC) = (/1.93841, 12.76830/)
   real(r8), parameter :: at_ppk = 21.67022  !at_ppk(NPOC) = (/29.27777, 21.67022/)
   real(r8), parameter :: bt_ppk = 0.21632  !bt_ppk(NPOC) = (/0.28991, 0.21632/)
   ! default sinking velocity (m/s)
   real(r8), parameter :: Vs0 = 8.68d-6  !Vs0(NPOC) = (/9.838d-8, 8.68d-6/)
   ! temperature multiplier for metabolic loss
   real(r8), parameter :: ThetaML = 1.073
   ! fraction of respiratioin relative to total metabolic loss
   real(r8), parameter :: Fres = 0.5  !Fres(NPOC) = (/0.8, 0.5/)
   ! fraction of excretion relative to non-respiration metabolic loss
   real(r8), parameter :: Fdom = 0.5  !Fdom(NPOC) = (/0.1, 0.5/)
   ! temperature multiplier for DOC microbial mineralization
   real(r8), parameter :: ThetaCM = 1.073
   ! The O2 half-saturation constant for DOC microbial degradation (umol/m3)
   real(r8), parameter :: Ko2CM = 4.6875d+4
   ! The O2 half-saturation constant for POM degradation (umol/m3)
   !real(r8), parameter :: Ko2PM = 1.5625d+4
   ! temperature multiplier for POC decomposition
   real(r8), parameter :: ThetaPM = 1.073
   ! algae mortality rate due to hypoxia (day-1)
   real(r8), parameter :: Kdhyp = 0.8
   ! The CO2 half-saturation constant for photosynthesis (umol/m3)
   real(r8), parameter :: Kco2 = 6.163d+4
   ! stoichemistry of C:P in POM and in allochthonous DOM
   real(r8), parameter :: YC2P_POM = 106.0
   real(r8), parameter :: YC2P_DOM = 199.0
   ! Redfield C : DOM mass ratio (C106H175O42N16P) (g/g)
   real(r8), parameter :: YC2DOM = 0.5358_r8 
   ! Mole mass of DOM (g/mol)
   real(r8), parameter :: MasDOM = 2374.0_r8 
   ! The rate coefficient of density change to irradiance (kg/m3/s)
   real(r8), parameter :: dsc1 = 7.816d-3  !dsc1(NPOC) = (/7.155d-3, 7.816d-3/)
   ! The minimum rate of density increase (kg/m3/s)
   real(r8), parameter :: dsc3 = 3.83d-4
   ! The concentration of DOC in precipitation (gC/m3)
   real(r8), parameter :: DOCrf = 1.0 ! Lin has changed from 2.0 according to Li et al., 2017
   ! aerial OC loading factor (gC/m/d)
   real(r8), parameter :: DOCae = 2.0
   ! ratio of DIC to DOC loading
   real(r8), parameter :: rDIC2DOC = 0.25
   ! The ratio of Chla:POC in streams (mg chla mmol C-1)
   real(r8), parameter :: rChla2POC = 0.036
   ! POC resuspension rate (umol/m2/s), if for chla this rate ranges
   ! at 0.1 to 1 mg/m2/d [Saloranta and Anderson, 2007]
   real(r8), parameter :: Scsed = 5.787d-3 
   ! flocculation rate (s-1)
   real(r8), parameter :: floc = 2.2d-8
   ! SOC(kg C/m2) added by Lin according to Wei et al., 2022
   real(r8), parameter :: soc(NLAKTYPE) = (/25.27,15.91,9.07,7.25/)
   real(r8), parameter :: salinity(NLAKTYPE) = (/0.05, 0.20, 0.85, 0.75/)
   ! according to Gui, 2023, master thesis
   real(r8), parameter :: AlphaChl = 7.5d-6   ! gC m2 (umol photonss g Chla)-1
   real(r8), parameter :: Ea = 7.0d4  ! J mol-1
   real(r8), parameter :: Ks = 0.1    ! uM
   real(r8), parameter :: PCmT = 2.6d-5    ! s-1 
   real(r8), parameter :: PQ10(NLAKTYPE) = (/4.91,5.33,3.82,1.86/) ! according to Chi, 2023 ! & Xu et al., 2024 (not used anymore)
 !  real(r8), parameter :: PQ10(NLAKTYPE) = (/2.54,6.60,2.37,1.90/)
   ! added by Lin, according to Gu, 2021
   real(r8), parameter :: OQ10 = 1.7 ! Q10 for CH4 oxidation
   real(r8), parameter :: AOMQ10(NLAKTYPE) = (/2.51,4.05,2.68,1.50/)!(/11.18,5.18,2.64,2.03/) ! according to Xu et al., 2024  //changed: considering Chi, 2023
   !real(r8), parameter :: AOMQ10(NLAKTYPE) = (/2.50,4.21,1.69,1.03/)
   real(r8), parameter :: Pco2Q10 = 2.7 ! Q10 for CO2 production due to CH4 oxidation
   real(r8), parameter :: Pn2oQ10 = 1.1 ! Q10 for N2O production due to CH4 oxidation
   ! according to Mu et al., 2023
   real(r8), parameter :: LKSO4(NLAKTYPE) = (/9.44, 25.28, 172.67, 183.95/)  ! SO42- in water, mg/L
   real(r8), parameter :: SedSO4(NLAKTYPE) = (/0.09, 0.08, 0.53, 0.23/) ! SO42- in sediment, mg/g, Chi, 2023
   real(r8), parameter :: SedNO3(NLAKTYPE) = (/0.019, 0.018, 0.007, 0.01/) ! NO3- in sediment, mg/g, Chi, 2023
   real(r8), parameter :: SedTN(NLAKTYPE) = (/0.067, 0.038, 0.036, 0.031/) ! TN(%)
   real(r8), parameter :: LKCl(NLAKTYPE) = (/54.45,119.97,478.46,364.45/) ! Cl- (mg/L)
   real(r8), parameter :: C2N(NLAKTYPE) = (/8.06,6.13,5.47,4.47/)
   real(r8), parameter :: LKNO3(NLAKTYPE) = (/1.25, 1.61, 0.98, 0.42/) !mg/L
   real(r8), parameter :: PLFA(NLAKTYPE) = (/2102.25, 804.94, 627.56, 524.82/) ! Total phospholipid fatty acids, ug/g 
   real(r8), parameter :: mcrA(NLAKTYPE) = (/18,12,13,10/) ! according to Mu et al., 2023
   real(r8), parameter :: pmoA(NLAKTYPE) = (/20,18,17,15/)
end module bg_const_mod
