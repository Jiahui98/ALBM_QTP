module shr_param_mod
!----------------------------------------------------------------------------
! Sensitive parameters list below (Param Id should be identical with 
! optpar.dat)
!
!1 Param_Ks  : thermal conductivity of soil constituent (W/(m*K))
!2 Param_Cps : heat capacity of soil constituent (J/(kg*K))
!3 Param_Por : porosity of lake sediment (0~1)
!4 Param_Rous: density of soil solid particle (kg/m3)
!5 BetaCH4   : CH4 concentration exponent coefficient for methanotrophy (0.5~1.5)
!6 LamO2     : O2 inhibition coefficient for methanotrophy (2~20)
!7 Kan       : Anaerobic oxidation to production rate ratio (0.15~0.27)
!8 Param_Re  : ebullition rate (unit: s-1)
!9 Param_DMP : recalcitrant organic matter dampening rate (m-1)
!10 Param_Rca : the fraction of aerobic decomposed carbon (s-1)
!11 Param_Vch : Chla-specific light saturated growth rate (mg C mg Chl-1 d-1)
!12 Param_Klr : metabolic loss rate coefficient (day-1)
!13 Param_RDOMaq: aquatic DOM microbial degradation rate (d-1)
!14 Param_RDOCtr: terrestrail DOM microbial degradation rate (d-1)
!15 Param_DOCwt: groundwater DOC concentration (mol/m3) 
!16 Param_Roun: snow density (kg/m3)
!17 Param_Feta: light attenuation correction factor for chla and CDOM
!18 Param_Wstr: wind shielding factor of mixing 
!19 Param_Ktscale: turbulence diffusivity scaling factor 
!20 Param_Hscale: sensible and latent heat transfer coefficent scaling factor
!----------------------------------------------------------------------------
   use shr_kind_mod,    only : r8, cx => SHR_KIND_CX 
   
   implicit none

   integer, parameter :: NPARAM = 20
   ! thermal related parameters
   integer, parameter :: Param_Ks = 1, Param_Cps = 2, Param_Por = 3, Param_Rous = 4
   ! methane oxidation related parameters
  ! integer, parameter :: !Param_OQ10 = 5, Param_Qch4 = 6,  &
                         !Param_Kch4 = 7, Param_Ko2 = 8
   integer, parameter :: Param_BetaCH4 = 5, Param_LamO2 = 6 ! changed according to Tan et al., 2024
   integer, parameter :: Param_Kan = 7 ! added by Lin 20250903
                                ! 0.15~0.27 according to Xu et al., 2024; Feng, 2018.
   ! 14C-enriched pool related parameters
   integer, parameter :: Param_Re = 8, &!, Param_PQ10n = 10, Param_Rcn = 11, & 
                         Param_DMP = 9
   ! aerobic decomposition parameters
   integer, parameter :: Param_Rca = 10
   ! phytoplankton related parameters    ! deleted by Lin
   integer, parameter :: Param_Vch = 11 !l = 15, Param_Vchs = 14 
   integer, parameter :: Param_Klr = 12 !l = 17, Param_Klrs = 16
   ! allochthonous carbon parameters
   integer, parameter :: Param_RDOMaq = 13, Param_RDOMtr = 14
   integer, parameter :: Param_DOCwt = 15
   ! photo-response parameters     ! deleted by Lin
 !  integer, parameter :: Param_phAlphas = 21, Param_phAlphal = 22
 !  integer, parameter :: Param_phBetal = 24, Param_phBetas = 23 
   ! phosphorus limitation     ! deleted by Lin
 !  integer, parameter :: Param_Ksrps = 25, Param_Ksrpl = 26
   ! ice phenology
   integer, parameter :: Param_Roun = 16
   ! thermodynamics
   integer, parameter :: Param_Feta = 17
   integer, parameter :: Param_Wstr = 18, Param_Ktscale = 19
   integer, parameter :: Param_Hscale = 20
   ! sensitive parameter collection
   real(r8) :: sa_params(NPARAM)
   ! maximum sample size
   integer :: NMAXSAMPLE
   ! first and last sample id to run
   integer :: sample_range(2)
   ! observation file and standard deviation
   character(cx) :: obs_dir, obs_var
   character(cx) :: obs_weight
   ! Monte Carlo sample file and result file
   character(cx) :: mc_file, sa_file

contains
   subroutine LoadSensitiveParameters(params)
      implicit none
      real(r8), intent(in) :: params(NPARAM)

      sa_params = params
   end subroutine

end module shr_param_mod
