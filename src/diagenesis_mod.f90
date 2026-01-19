module diagenesis_mod
!---------------------------------------------------------------------------------
! Purpose: this module governs the diagenetic reactions of organic carbon in the 
!          sediments, which converts organic carbon to either CO2 and CH4.
!
!          [Stepanenko et al., 2011; Izvestiya, Atmospheric and 
!           Oceanic Physics], [Kessler et al., 2012; JGR] and [Hanson et al.,
!           2011; PLus One].
!
!---------------------------------------------------------------------------------
   use shr_ctrl_mod, inft => INFINITESIMAL_E8
   use shr_typedef_mod,       only : RungeKuttaCache2D
   use shr_param_mod
   use phy_utilities_mod
   use bg_utilities_mod
   use data_buffer_mod
   use bg_const_mod

   implicit none
   private
   public :: InitializeDiagenesisModule, DestructDiagenesisModule
   public :: DiagenesisModuleSetup, DiagenesisModuleCallback 
   public :: DiagenesisEquation, GetSedCLossRate 
   !public :: ConstructActCarbonPool  ! Lin has closed this subroutine, as 14C-depleted C & 
   ! participates little in QTP thrmkst lakes.
   public :: mem_ch4
   ! memory cache for Runge-Kutta
   type(RungeKuttaCache2D) :: mem_ch4
   ! CH4 anaerobic production rates (umol C m-3 s-1)
   real(r8), allocatable :: rPCH4(:)  ! changed by Lin, so as follows
   ! CH4 aerobic and anaerobic oxidation rates (umol C m-3 s-1)
   real(r8), allocatable :: rOCH4(:), rOACH4(:), rONCH4(:)
   ! CO2 anaerobic production rates (umol C m-3 s-1)
   real(r8), allocatable :: rnPCO2(:)
   ! CO2 aerobic production rates (umol C m-3 s-1)
   real(r8), allocatable :: rPCO2(:)
   ! water substance dynamcis (umol/m3/s)
   real(r8), allocatable :: rDYN(:,:)
   ! gas ebullition rates in the sediments
   real(r8), allocatable :: Ebch4(:,:)
   ! lagged ice percentage
   real(r8), allocatable :: lagged_ice(:)
   ! diffusivity in the sediments (m2/s)
   real(r8), allocatable :: Dch4(:)
   ! passive pool mobilization rate (umol/m3/s)
   real(r8) :: pas2actC

contains
   subroutine InitializeDiagenesisModule()
      implicit none

      allocate(rDYN(NSSUB,NSLAYER+1))
      allocate(rPCH4(NSLAYER+1))   ! deleted NPOOL by Lin, so as follows
      allocate(rOCH4(NSLAYER+1))
      allocate(rOACH4(NSLAYER+1))
      allocate(rONCH4(NSLAYER+1))
      allocate(rPCO2(NSLAYER+1))
      allocate(rnPCO2(NSLAYER+1))
      allocate(Ebch4(NSSUB,NSLAYER+1))
      allocate(lagged_ice(NSLAYER+1))
      allocate(Dch4(NSLAYER+1))
      allocate(mem_ch4%K1(NSSUB,NSLAYER+1))
      allocate(mem_ch4%K2(NSSUB,NSLAYER+1))
      allocate(mem_ch4%K3(NSSUB,NSLAYER+1))
      allocate(mem_ch4%K4(NSSUB,NSLAYER+1))
      allocate(mem_ch4%K5(NSSUB,NSLAYER+1))
      allocate(mem_ch4%K6(NSSUB,NSLAYER+1))
      allocate(mem_ch4%nxt4th(NSSUB,NSLAYER+1))
      allocate(mem_ch4%nxt5th(NSSUB,NSLAYER+1))
      allocate(mem_ch4%interim(NSSUB,NSLAYER+1))
      allocate(mem_ch4%rerr(NSSUB,NSLAYER+1))

      call InitializeDiagenesisStateVariables()
      Ebch4 = 0.0_r8
      rDYN = 0.0_r8
      rPCH4 = 0.0_r8
      rOCH4 = 0.0_r8
      rOACH4 = 0.0_r8
      rONCH4 = 0.0_r8
      rPCO2 = 0.0_r8
      rnPCO2 = 0.0_r8
      Dch4 = 0.0_r8
      pas2actC = 0.0_r8
   end subroutine

   subroutine DestructDiagenesisModule()
      implicit none

      deallocate(rDYN)
      deallocate(rPCH4)
      deallocate(rOCH4)
      deallocate(rOACH4)
      deallocate(rONCH4)
      deallocate(rPCO2)
      deallocate(rnPCO2)
      deallocate(Ebch4)
      deallocate(lagged_ice)
      deallocate(Dch4)
      deallocate(mem_ch4%K1)
      deallocate(mem_ch4%K2)
      deallocate(mem_ch4%K3)
      deallocate(mem_ch4%K4)
      deallocate(mem_ch4%K5)
      deallocate(mem_ch4%K6)
      deallocate(mem_ch4%nxt4th)
      deallocate(mem_ch4%nxt5th)
      deallocate(mem_ch4%interim)
      deallocate(mem_ch4%rerr)
   end subroutine
   
   subroutine DiagenesisModuleSetup()
      implicit none

      call CalcDiagenesisRates()
      call UpdateBubbleFlux()
      call UpdateSubDiffusivity() 
   end subroutine

   subroutine DiagenesisModuleCallback(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: Porosity

      Porosity = sa_params(Param_Por)
      lagged_ice = m_sedIce / Porosity
      call AdjustNegativeConcentration()
      ! Add wind- and thermo-karst eroded C to surface C pool
      ! units: umol/m3/s
      !call UpdateActCarbonPool(dt, 0d0) ! closed by Lin
      call UpdatePasCarbonPool(dt, 0d0)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Govern the production of methane bubbles
   !
   !------------------------------------------------------------------------------
   subroutine DiagenesisEquation(con, dcon)
      implicit none
      real(r8), intent(in) :: con(NSSUB,NSLAYER+1)
      real(r8), intent(out) :: dcon(NSSUB,NSLAYER+1)
      real(r8), dimension(NSSUB) :: qb, qt, dC1, dC2
      real(r8) :: a, b
      integer :: ii, Top_Index, Bottom_Index, jj

      Top_Index = m_sedWaterTopIndex
      Bottom_Index = m_sedWaterBtmIndex
      if (Top_Index>Bottom_Index) then
         dcon = 0.0_r8
         print *,"Top>Bot"
         return
      end if
      call GetTopBoundaryFlux(con(:,1), qt)
      qb = 0.0_r8
      do ii = 1, NSLAYER+1, 1
         if (ii==1) then
            a = 0.5*(Dch4(ii) + Dch4(ii+1))
            dC1 = (con(:,ii+1) - con(:,ii)) / (m_Zs(ii+1) - m_Zs(ii))
            dcon(:,ii) = (a*dC1 - qt) / m_dZs(ii) + rDYN(:,ii) &
                     - Ebch4(:,ii)    
         else if (ii==NSLAYER+1) then
            b = 0.5*(Dch4(ii-1) + Dch4(ii))
            dC2 = (con(:,ii) - con(:,ii-1)) / (m_Zs(ii) - m_Zs(ii-1))
            dcon(:,ii) = (qb - b*dC2) / m_dZs(ii) + rDYN(:,ii) &
                     - Ebch4(:,ii)
         else
            a = 0.5*(Dch4(ii) + Dch4(ii+1))
            b = 0.5*(Dch4(ii-1) + Dch4(ii))
            dC1 = (con(:,ii+1) - con(:,ii)) / (m_Zs(ii+1) - m_Zs(ii))
            dC2 = (con(:,ii) - con(:,ii-1)) / (m_Zs(ii) - m_Zs(ii-1))
            dcon(:,ii) = (a*dC1 - b*dC2) / m_dZs(ii) + rDYN(:,ii) &
                     - Ebch4(:,ii)
         end if
      end do
      where(con<=0 .and. dcon<0) dcon = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Adjust negative sediment gas concentration (project negative values 
   !          to zero) - "Non-negative solutions of ODEs" (Shampine, L., 2005).
   !
   !------------------------------------------------------------------------------
   subroutine AdjustNegativeConcentration()
      implicit none

      where (m_sedSubCon<0) m_sedSubCon = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate substance prod/loss rates
   !
   !------------------------------------------------------------------------------
   subroutine CalcDiagenesisRates()
      implicit none
      real(r8) :: Prtot, Prnet, Hch4, Hn2
      real(r8) :: Hco2, Ho2, fn2, fo2, fco2, fch4
      real(r8) :: Porosity, Re, Ps, Psat, pH
      real(r8) :: pch4, pn2, po2, pco2, temp, Do2, depth, Sedch4 
      real(r8) :: pool   ! deleted'(NPOOL)' by Lin
      integer :: ii
      real(r8) :: gammaT, gamma_tref, Pcap_tref, Pcap

      Re = sa_params(Param_Re)
      Porosity = sa_params(Param_Por)
      Ps = m_surfData%pressure
      pH = SedpH(lake_info%itype)
      !salinity = salinity(lake_info%itype)
      do ii = 1, NSLAYER+1, 1
         if (m_sedIce(ii)>=Porosity) then
            rPCH4(ii) = 0.0_r8   ! changed by Lin, so as follows
            rPCO2(ii) = 0.0_r8
            rnPCO2(ii) = 0.0_r8
            rOCH4(ii) = 0.0_r8
            rOACH4(ii) = 0.0_r8
            rONCH4(ii) = 0.0_r8
            Ebch4(:,ii) = 0.0_r8
         else
            ! carbon mineralization rates
            temp = m_sedTemp(ii)
            Do2 = m_sedSubCon(Wo2,ii)
            Sedch4 = m_sedSubCon(Wch4, ii)
            pool = m_unfrzCarbPool(ii)  ! changed by Lin, as there is only one kind of C left
            call Methanogenesis(lake_info%itype, pool, temp, Do2, &
               rPCH4(ii), rnPCO2(ii))  ! added itype by Lin 
            call AerobicDegradation(pool, temp, Do2, rPCO2(ii))
            ! methanotrophy: aerobic/anaerobic
            if (Sedch4>minDch4) then
               call Methanotrophy(lake_info%itype, Sedch4, Do2, temp, &
                  rPCH4(ii), rOCH4(ii), rOACH4(ii), rONCH4(ii))
            end if
            
            ! ebullition rates (Eq. (A.1), Stepanenko et al.)
            depth = m_Zs(ii) + lake_info%depth
            Hch4 = CalcHenrySolubility(lake_info%itype, Wch4, temp, pH)
            Hn2 = CalcHenrySolubility(lake_info%itype, Wn2, temp, pH)
            Ho2 = CalcHenrySolubility(lake_info%itype, Wo2, temp, pH)
            Hco2 = CalcHenrySolubility(lake_info%itype, Wco2, temp, pH)
           ! Hco2 = exp(-58.0931+90.5069*(100/temp)+22.294*log(100/temp)+ &
           !    salinity(lake_info%itype)*(0.027766-0.025888*(temp/100)+ &
           !    0.0050578*(temp/100)**2)) * 1d9 /P0 
            pch4 = m_sedSubCon(Wch4,ii)/Hch4
            pn2 = m_sedSubCon(Wn2,ii)/Hn2
            po2 = m_sedSubCon(Wo2,ii)/Ho2
            pco2 = m_sedSubcon(Wco2,ii)/Hco2
            Prtot = pn2 + pch4 + po2 + pco2
            !Prnet = Prtot - Ae * Porosity * (Ps + Roul*G*depth)            
         !   Prtot = pn2 + pch4
            Psat = Ae * Porosity * (Ps+Roul*G*depth)
            
            ! introduce in Pcap - threshold pressure. Added by Lin, 20250904
            ! Incorporate capillary effects for triggering ebullition 
            ! Kellner et al., 2005; Tang et al., 2010
            gammaT = CalcSurfaceTension(temp)
            gamma_tref = CalcSurfaceTension(278.15_r8)
            Pcap_tref = Roul * G * 0.218 ! 21.8 cm for sandy loam, empirical from Clapp&Hornberger, 1978
                        ! Mu et al., 2025 & Cui et al., 2012, 2014: claim as sandy loam
                        ! Capillary intrusion pressure at base temperature
            Pcap = Pcap_tref * gammaT/max(gamma_tref, inft) ! Young-Laplace
                        ! threshold pressure
            Psat = Psat + Pcap ! Incorporating capillary.

            Prnet = max(Prtot-Psat, 0.0_r8)
            fn2 = pn2 / (Prtot + inft)
            fo2 = po2 / (Prtot + inft)
            fch4 = pch4 / (Prtot + inft)
            fco2 = pco2/ (Prtot + inft) 
            Ebch4(Wn2,ii) = Re * fn2 * Prnet * Hn2
            Ebch4(Wo2,ii) = Re * fo2 * Prnet * Ho2 
            Ebch4(Wch4,ii) = Re * fch4 * Prnet * Hch4
            Ebch4(Wco2,ii) = Re * fco2 * Prnet * Hco2
            Ebch4(Wsrp,ii) = 0.0_r8
         end if
      end do
      rDYN(Wn2,:) = 0.72*rOACH4 + 4*rONCH4  ! Gu et al., 2021 
      rDYN(Wo2,:) = -rPCO2 - 2.0*rOACH4      !-sum(rPCO2,1)    changed by Lin, so as follows
      rDYN(Wco2,:) = rPCO2 + rnPCO2+rOACH4    !sum(rPCO2,1) + sum(rnPCO2,1)
      rDYN(Wch4,:) = rPCH4-rOCH4     !sum(rPCH4,1)
      rDYN(Wsrp,:) = ( rPCH4-rOCH4 + rPCO2 + rnPCO2 ) / YC2P_POM
      !( sum(rPCH4,1) + sum(rPCO2,1) + sum(rnPCO2,1) ) / YC2P_POM
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: get top boundary flux.
   !
   !------------------------------------------------------------------------------
   subroutine GetTopBoundaryFlux(con, topsflux)
      implicit none
      real(r8), intent(in) :: con(NSSUB)        ! units: umol/m3
      real(r8), intent(out) :: topsflux(NSSUB)  ! units: umol/m2/s
      real(r8) :: Porosity ! added according to ALBM-update
      integer :: ii, bottom

      if (m_lakeWaterTopIndex>WATER_LAYER+1) then
         topsflux = 0.0_r8
      else
         ! top dissolved substance fluxes
         Porosity = sa_params(Param_Por) ! added according to ALBM-update
         bottom = WATER_LAYER + 1
         topsflux = m_Kbm / DeltaD * ( con/Porosity - &
            m_waterSubCon(Wn2:Wsrp,bottom) ) ! changed according to ALBM-update
         do ii = Wn2, Wo2, 1
            topsflux(ii) = min(0.0, topsflux(ii))
         end do
         do ii = Wco2, Wsrp, 1
            topsflux(ii) = max(0.0, topsflux(ii))
         end do
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update bubble flux rates from sediment top.
   !
   !------------------------------------------------------------------------------
   subroutine UpdateBubbleFlux()
      implicit none
      integer :: ii

      if (m_lakeWaterTopIndex>WATER_LAYER+1) then
         m_btmbflux = 0.0_r8
      else
         do ii = 1, NGAS, 1
            m_btmbflux(ii) = sum(Ebch4(ii,:)*m_dZs)
         end do
      end if
      !print *, m_btmbflux(Wn2),m_btmbflux(Wch4)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate dissolved substance diffusivity in sediments
   !
   !------------------------------------------------------------------------------
   subroutine UpdateSubDiffusivity()
      implicit none
      real(r8) :: Porosity, diffst
      integer :: ii

      Porosity = sa_params(Param_Por)
      do ii = 1, NSLAYER+1, 1
         if (m_sedIce(ii)>=Porosity) then
            diffst = CalcGasDiffusivityInIce(Wch4, m_sedTemp(ii))
            Dch4(ii) = diffst*0.66*Porosity 
         !else if (ii==1) then
         !   Dch4(ii) = 1 ! added by Lin
         else  
            ! Eq. (8) in Walter et al. (2000) gives Dch4 at the magnitude
            ! of 1e-9 but Goto et al. (2017; Marine Geophysical Research) shows
            ! this value at the magnitude of 1e-7
            Dch4(ii) = 2.0d-7*0.66*Porosity   ! [Eq. (8) in Walter 2000] ! changed from 9 to 7 according to ALBM-update
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Construct active and passive carbon pool profile in sediments.
   !          Active carbon pool: [West and plug, 2008 JGR]
   !          Passive carbon pool: [Stepanenko et al., 2011]
   !  
   !------------------------------------------------------------------------------
   !subroutine ConstructActCarbonPool()        ! closed by Lin
   !   implicit none
   !   real(r8) :: Ht, Zs, time, Ctotal
   !   real(r8) :: Ct, Rco
   !   integer  :: ii, itk
   !   
   !   if (lake_info%thrmkst/=2) then
   !      ! non-yedoma region
   !      m_frzCarbPool(actC,:) = 0.0_r8
   !      m_unfrzCarbPool(actC,:) = 0.0_r8
   !   else
   !      ! permafrost thawing rate (m/(yr^0.5))
   !      if (lake_info%margin==1) then
   !         Ct = 0.67_r8
   !      else
   !         Ct = 0.77_r8
   !      end if
   !      Rco = 2.149145d-10 
   !      do itk = NSLAYER+1, 1, -1
   !         if (m_sedIce(itk)<=0) then
   !            exit
   !         end if
   !      end do
   !      if (itk==NSLAYER+1) then
   !         m_frzCarbPool(actC,:) = 0.0_r8
   !         m_unfrzCarbPool(actC,:) = 0.0_r8
   !      else
   !         Ht = m_Zs(itk)
   !         do ii = 1, NSLAYER+1, 1
   !            Zs = m_Zs(ii)
   !            if (Zs+lake_info%depth>25) then
   !               m_frzCarbPool(actC,ii:NSLAYER+1) = 0.0_r8
   !               m_unfrzCarbPool(actC,ii:NSLAYER+1) = 0.0_r8
   !               exit
   !            end if
   !            if (ii<=itk) then
   !               ! time that talik has developed
   !               time = SECOND_OF_YEAR * (Ht**2 - Zs**2) / Ct**2
   !               Ctotal = oldcarb0 * exp(-Rco*time) ! max(0.0_r8, oldcarb0*(1.0-Rco*time)) changed according to ALBM-update
   !            else
   !               Ctotal = oldcarb0
   !            end if
   !            ! mean labile carbon density (umol/m3)
   !            Ctotal = Ctotal / 3.0 / MasC * 1.0d+9
   !            m_frzCarbPool(actC,ii) = Ctotal * lagged_ice(ii)
   !            m_unfrzCarbPool(actC,ii) = Ctotal - m_frzCarbPool(actC,ii)
   !         end do
   !      end if
   !   end if
   !end subroutine

   subroutine ConstructPasCarbonPool()
      implicit none
      real(r8), parameter :: soc_ref = 21.98_r8  ! catchment SOC of Quebec
      real(r8), parameter :: carb0 = 8.9386_r8     ! carbon constant
      real(r8) :: fsoc, fshape, fersn
      real(r8) :: labile, dampen 
      real(r8) :: carbon
      character(cx) :: strmsg
    !  integer :: ii

    !  allocate(carbon(NSLAYER+1))

      ! surface soil organic carbon density
      fshape = ( sqrt(1.0d-6*lake_info%Asurf) / lake_info%basin )**(-0.555) ! change depth into basin according to Ferland et al., 2012
      fsoc = m_soc / soc_ref
    !  if (lake_info%thrmkst==2 .and. lake_info%margin==1) then
    !     fersn = 1.0d0 / 2.0d0
    !     labile = 0.05
    !  else if (lake_info%thrmkst==1 .and. lake_info%margin==1) then
      fersn = 1.0d0 / 6.0d0
      labile = 1.0d0 !0.05 changed by Lin, 20250903 after meeting
    !  else
    !     fersn = 0.0_r8
    !     labile = 0.05
    !  end if
      !carbon = carbon*1d9/MasC      ! convert kg/m2 to umol/m2
      
      carbon = (labile*carb0*fshape*fsoc + m_soc*fersn)*1d9/MasC
      pas2actC = 4d-3*carb0*fshape*fsoc*1d9/MasC/SECOND_OF_YEAR
      dampen = sa_params(Param_DMP)
      if (dampen>0) then
         carbon = carbon * dampen / (1.0 - exp(-dampen))
         pas2actC = pas2actC * dampen / (1.0 - exp(-dampen))
         m_frzCarbPool(:) = carbon * lagged_ice * exp(-dampen*m_Zs)
         m_unfrzCarbPool(:) = carbon * (1.0 - lagged_ice) * &
                              exp(-dampen*m_Zs)
      else
         m_frzCarbPool(:) = carbon * lagged_ice
         m_unfrzCarbPool(:) = carbon * (1.0 - lagged_ice)
      end if
      m_frzCarbPool = carbon * lagged_ice
      m_unfrzCarbPool = carbon * (1.0 - lagged_ice)
     ! deallocate(carbon)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update the active and passive carbon profile in sediments.
   !          Active carbon pool: [West and plug, 2008 JGR]
   !          Passive carbon pool: [Stepanenko et al., 2011]
   !  
   !------------------------------------------------------------------------------
   !subroutine UpdateActCarbonPool(dt, rDOCer)     closed by Lin
   !   implicit none
   !   real(r8), intent(in) :: dt       ! time interval (sec)
   !   real(r8), intent(in) :: rDOCer   ! DOC gain by erosion (umol/m3/s)
   !   real(r8) :: frice_lag, frice, ratio
   !   real(r8) :: Porosity, Ctotal, decay
   !   real(r8) :: cDOC
   !   integer  :: ii
   !
   !   if (lake_info%thrmkst/=2) then
   !      return   ! non-yedoma region
   !   end if
   !   ! Update C pool due to erosion and surface exchange
   !   m_unfrzCarbPool(actC,1) = m_unfrzCarbPool(actC,1) + &
   !      rDOCer*dt
   !   ! Shrink C pool due to degradation
   !   Porosity = sa_params(Param_Por)
   !   do ii = 1, NSLAYER+1, 1
   !      Ctotal = m_unfrzCarbPool(actC,ii) + m_frzCarbPool(actC,ii)
   !      if (Ctotal>e8) then
   !         decay = rPCH4(actC,ii) + rnPCO2(actC,ii) + rPCO2(actC,ii)
   !         m_unfrzCarbPool(actC,ii) = m_unfrzCarbPool(actC,ii) - decay*dt
   !         m_unfrzCarbPool(actC,ii) = max(0.0_r8, m_unfrzCarbPool(actC,ii))
   !         ! redistribute frozen and unfrozen carbon pools
   !         frice = m_sedIce(ii) / Porosity 
   !         frice_lag = lagged_ice(ii)
   !         if (frice>frice_lag) then
   !            ! ice expand and freezing carbon
   !            ratio = (frice - frice_lag) / (1.0 - frice_lag)
   !            m_frzCarbPool(actC,ii) = m_frzCarbPool(actC,ii) + ratio * &
   !                                 m_unfrzCarbPool(actC,ii)
   !            m_unfrzCarbPool(actC,ii) = m_unfrzCarbPool(actC,ii) * &
   !                                 (1.0 - ratio)
   !         else if (frice<frice_lag) then
   !            ! ice shrink and thawing carbon
   !            ratio = (frice_lag - frice) / frice_lag
   !            m_unfrzCarbPool(actC,ii) = m_unfrzCarbPool(actC,ii) + &
   !                                    ratio * m_frzCarbPool(actC,ii)
   !            m_frzCarbPool(actC,ii) = m_frzCarbPool(actC,ii) * &
   !                                 (1.0 - ratio)
   !         end if
   !      else
   !         m_frzCarbPool(actC,ii) = 0.0_r8
   !         m_unfrzCarbPool(actC,ii) = 0.0_r8
   !      end if
   !   end do
   !end subroutine

   subroutine UpdatePasCarbonPool(dt, rDOCer)
      implicit none
      real(r8), intent(in) :: dt       ! time step (sec)
      real(r8), intent(in) :: rDOCer   ! DOC gain by erosion (umol/m3/s)
      real(r8) :: Ctotal, Porosity
      real(r8) :: ratio, decay, dampen
      real(r8) :: frice_lag, frice
      integer :: ii

      Porosity = sa_params(Param_Por)
      dampen = sa_params(Param_DMP)
      ! update C pool due to erosion
      m_unfrzCarbPool(1) = m_unfrzCarbPool(1) + rDOCer*dt
      ! aerobic and anaerobic mineralization and freeze/thaw
      do ii = 1, NSLAYER+1, 1
         Ctotal = m_unfrzCarbPool(ii) + m_frzCarbPool(ii)
         if (Ctotal>e8) then
            decay = rPCH4(ii) + rnPCO2(ii) + rPCO2(ii) 
            m_unfrzCarbPool(ii) = m_unfrzCarbPool(ii) - &
               decay*dt + pas2actC*exp(-dampen*m_Zs(ii))*dt
            m_unfrzCarbPool(ii) = max(0.0_r8, m_unfrzCarbPool(ii))
            ! redistribute frozen and unfrozen carbon pools
            frice = m_sedIce(ii) / Porosity
            frice_lag = lagged_ice(ii)
            if (frice>frice_lag) then
               ! ice expand and freezing carbon
               ratio = (frice - frice_lag) / (1.0 - frice_lag)
               m_frzCarbPool(ii) = m_frzCarbPool(ii) + ratio * m_unfrzCarbPool(ii)
               m_unfrzCarbPool(ii) = m_unfrzCarbPool(ii) * (1.0 - ratio)
            else if (frice<frice_lag) then     
               ! ice shrink and thawing carbon
               ratio = (frice_lag - frice) / frice_lag
               m_unfrzCarbPool(ii) = m_unfrzCarbPool(ii) + ratio * m_frzCarbPool(ii)
               m_frzCarbPool(ii) = m_frzCarbPool(ii) * (1.0 - ratio)
            end if
         else
            m_frzCarbPool(ii) = 0.0_r8
            m_unfrzCarbPool(ii) = 0.0_r8
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate CH4/CO2 production rate.
   !
   !------------------------------------------------------------------------------
   subroutine GetSedCLossRate(pco2o, pco2n, pch4o, pch4n)
      implicit none
      real(r8), intent(out) :: pco2o   ! units: umol/m2/s
      real(r8), intent(out) :: pco2n   ! units: umol/m2/s
      real(r8), intent(out) :: pch4o   ! units: umol/m2/s
      real(r8), intent(out) :: pch4n   ! units: umol/m2/s

      pco2n = sum( (rPCO2(:) + rnPCO2(:)) * m_dZs )
      !pco2o = sum( (rPCO2(actC,:) + rnPCO2(actC,:)) * m_dZs )  ! closed by Lin, so as follows
      pch4n = sum( (rPCH4(:)-rOCH4(:)) * m_dZs )
      !pch4o = sum( rPCH4(actC,:) * m_dZs )
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Initialize sediment gas profiles
   !
   !------------------------------------------------------------------------------
   subroutine InitializeDiagenesisStateVariables()
      implicit none
      real(r8) :: ps, temp, Porosity, pH
      integer :: ii

      if (m_radPars%month==7 .or. m_radPars%month==8) then
         pH = SMpH(lake_info%itype)
      else
         pH = LKpH(lake_info%itype) 
      end if
      Porosity = sa_params(Param_Por)
      do ii = 1, NSLAYER+1, 1
         temp = m_sedTemp(ii)
         ps = m_radPars%qCH4 * P0    ! CH4 
         m_sedSubCon(Wch4,ii) = Porosity*CalcEQConc(lake_info%itype,Wch4,temp,pH,ps)
         ps = 1.0d-6 * m_radPars%qCO2 * P0   ! CO2
         m_sedSubCon(Wco2,ii) = Porosity*CalcEQConc(lake_info%itype,Wco2,temp,pH,ps)
      !   ps = Xo2 * P0    ! o2
      !   m_sedSubCon(Wo2,ii) = Porosity*CalcEQConc(lake_info%itype,Wo2,temp,pH,ps)
      end do
      m_sedSubCon(Wsrp,:) = 1.0 / YC2P_POM * (m_sedSubCon(Wco2,:) + &
                              m_sedSubCon(Wch4,:))
      m_sedSubCon(Wn2,:) = 0.0_r8   ! initially no n2
      m_sedSubCon(Wo2,:) = 0.0_r8   ! initially no o2
      m_btmbflux = 0.0_r8

      ! initialize 14C-depleted and 14C-enrich carbon pool 
      lagged_ice = m_sedIce / Porosity
      call ConstructPasCarbonPool()
      !call ConstructActCarbonPool()
   end subroutine

end module diagenesis_mod
