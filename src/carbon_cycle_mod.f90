module carbon_cycle_mod
!---------------------------------------------------------------------------------
! Purpose: this module governs the concentration of water dissolved substances.
!          
!          References include Tang's thesis; Stepanenko et al., 2011;
!           and Hanson et al., 2004, 2011 & 2014.
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,          only : r8
   use shr_ctrl_mod,          only : WATER_LAYER, NGAS, NWSUB, lake_info, &
                                     e8 => SHR_CTRL_E8
   use shr_typedef_mod,       only : RungeKuttaCache2D, RungeKuttaCache1D ! added 1D by Lin
   use shr_param_mod
   use phy_utilities_mod 
   use bg_utilities_mod
   use data_buffer_mod

   implicit none
   private
   public :: mem_sub, mem_poc, top_fluxes
   public :: InitializeCarbonModule, DestructCarbonModule
   public :: CarbonModuleSetup, CarbonModuleCallback 
   public :: CarbonCycleEquation, ParticulateEquation 
   public :: GetO2ProductionRate, GetCO2ProductionRate
   public :: GetpCO2MixingRatio
   ! memory cache for Runge-Kutta
   type(RungeKuttaCache2D) :: mem_sub
   type(RungeKuttaCache1D) :: mem_poc
   ! gas transfer velocity
   real(r8), allocatable :: Kg(:)
   ! top surface exchange (umol/m2/s)
   real(r8), allocatable :: top_fluxes(:)
   ! particulate settling velocity (m/s)
   real(r8), allocatable :: Vsettl(:)
   ! algae density (kg/m3)
   real(r8), allocatable :: RhoA(:)
   ! DOC photo degradation (umol/m3/s)
   real(r8), allocatable :: rCDOM(:)
   ! heterotrophic respiration (umol/m3/s)
   real(r8), allocatable :: rRDOC(:,:)
   ! carbon fixation rate (umol/m3/s)
   real(r8), allocatable :: rGPP(:), rGPP_hourly(:)!added by Lin 20250903
   ! phytoplankton metabolic loss (umol/m3/s)
   real(r8), allocatable :: rLPOC(:)
   real(r8), allocatable :: rRLPOC(:)
   ! CH4 oxidation (umol/m3/s)
   real(r8), allocatable :: rOACH4(:)!, rONCH4(:), rOCH4(:)
   ! CH4 production (umol/m3/s)
   real(r8), allocatable :: rPCH4(:), rnPCO2(:)
   ! Substance exchange with side boundary (umol/m3/s)
   real(r8), allocatable :: rBdExg(:,:)
   ! water substance dynamcis (umol/m3/s) 
   real(r8), allocatable :: rDYN(:,:)
   ! POC dynamics (umol/m3/s)
   real(r8), allocatable :: rPDYN(:)
   ! PAR radiation (mol/m2/s)
   real(r8), allocatable :: Ipar(:)
   ! PAR radiation (W/m2)
   real(r8), allocatable :: pPAR(:)
   ! Chla : C ratio (mg chla mmol C-1)
   real(r8), allocatable :: rChl2C(:), rC2Chl(:)    ! added by Lin
   ! POC resuspension
   real(r8) :: RsPOCLoad
   ! DOCtr flocculation rate (umol/m3/s)
   real(r8), allocatable :: rfloc(:)
   ! DOC wetland and erosion loads
   real(r8), allocatable :: WtCLoad(:)
   ! DOC, DIC, P, DO and POC streamflow loads
   real(r8), allocatable :: SwDOCLoad(:)
   real(r8), allocatable :: SwDICLoad(:)
   real(r8), allocatable :: SwPOCLoad(:)
   real(r8), allocatable :: SwSRPLoad(:)
   real(r8), allocatable :: SwDOLoad(:)
   ! stream out flow rate (s-1)
   real(r8), allocatable :: SwQso(:)
   ! DOC Aerial and rainfall loads and P deposition
   real(r8) :: AeCLoad, RfCLoad
   real(r8) :: AePLoad
   real(r8), allocatable :: PC(:), PCm(:)   ! added by Lin

contains
   subroutine InitializeCarbonModule()
      implicit none

      allocate(Kg(NGAS))
      allocate(top_fluxes(NWSUB))
      !allocate(RsPOCLoad(NPOC))   deleted by Lin
      allocate(Vsettl(WATER_LAYER+1))
      allocate(RhoA(WATER_LAYER+1))
      allocate(rDYN(NWSUB,WATER_LAYER+1))
      allocate(rPDYN(WATER_LAYER+1))
      allocate(rBdExg(NWSUB,WATER_LAYER+1))
      allocate(rCDOM(WATER_LAYER+1))
      allocate(rGPP(WATER_LAYER+1))
      allocate(rGPP_hourly(WATER_LAYER+1)) ! added by Lin 20250903
      allocate(rRDOC(NDOC,WATER_LAYER+1))
      allocate(rLPOC(WATER_LAYER+1))
      allocate(rRLPOC(WATER_LAYER+1))   
     ! allocate(rOCH4(WATER_LAYER+1))
      allocate(rOACH4(WATER_LAYER+1))
     ! allocate(rONCH4(WATER_LAYER+1))
      allocate(rPCH4(WATER_LAYER+1))
      allocate(rnPCO2(WATER_LAYER+1))
      allocate(rfloc(WATER_LAYER+1))
      allocate(Ipar(WATER_LAYER+1))
      allocate(pPAR(WATER_LAYER+1))
      allocate(rChl2C(WATER_LAYER+1))
      allocate(rC2Chl(WATER_LAYER+1))
      allocate(WtCLoad(WATER_LAYER+1))
      allocate(SwDOCLoad(WATER_LAYER+1))
      allocate(SwDICLoad(WATER_LAYER+1))
      allocate(SwPOCLoad(WATER_LAYER+1))
      allocate(SwSRPLoad(WATER_LAYER+1))
      allocate(SwDOLoad(WATER_LAYER+1))
      allocate(SwQso(WATER_LAYER+1))
      allocate(mem_sub%K1(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K2(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K3(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K4(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K5(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K6(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%nxt4th(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%nxt5th(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%interim(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%rerr(NWSUB,WATER_LAYER+1))
      allocate(mem_poc%K1(WATER_LAYER+1))
      allocate(mem_poc%K2(WATER_LAYER+1))
      allocate(mem_poc%K3(WATER_LAYER+1))
      allocate(mem_poc%K4(WATER_LAYER+1))
      allocate(mem_poc%K5(WATER_LAYER+1))
      allocate(mem_poc%K6(WATER_LAYER+1))
      allocate(mem_poc%nxt4th(WATER_LAYER+1))
      allocate(mem_poc%nxt5th(WATER_LAYER+1))
      allocate(mem_poc%interim(WATER_LAYER+1))
      allocate(mem_poc%rerr(WATER_LAYER+1))
      allocate(PC(WATER_LAYER+1))
      allocate(PCm(WATER_LAYER+1))

      call InitializeCarbonStateVariables()
      rCDOM = 0.0_r8
      rGPP = 0.0_r8
      rGPP_hourly = 0.0_r8 ! added by Lin 20250903
      rRDOC = 0.0_r8
      rLPOC = 0.0_r8
      rRLPOC = 0.0_r8
     ! rOCH4 = 0.0_r8
      rOACH4 = 0.0_r8
     ! rONCH4 = 0.0_r8
      rPCH4 = 0.0_r8
      rnPCO2 = 0.0_r8
      rDYN = 0.0_r8
      rPDYN = 0.0_r8
      rBdExg = 0.0_r8
      rfloc = 0.0_r8
      Kg = 0.0_r8
      top_fluxes = 0.0_r8
      Vsettl = 0.0_r8
      Ipar = 0.0_r8
      pPAR = 0.0_r8
      rChl2C = 0.24_r8
      !rC2Chl = 0.0_r8
      PCm = 0.0_r8
      PC = 0.0_r8
      RsPOCLoad = 0.0_r8
      RhoA = Roul
      AeCLoad = 0.0_r8
      AePLoad = 0.0_r8
      WtCLoad = 0.0_r8
      RfCLoad = 0.0_r8
      SwDOCLoad = 0.0_r8
      SwDICLoad = 0.0_r8
      SwSRPLoad = 0.0_r8
      SwPOCLoad = 0.0_r8
      SwDOLoad = 0.0_r8
      SwQso = 0.0_r8
   end subroutine

   subroutine DestructCarbonModule()
      implicit none

      deallocate(Kg)
      deallocate(top_fluxes)
      !deallocate(RsPOCLoad)
      deallocate(Vsettl)
      deallocate(RhoA)
      deallocate(rDYN)
      deallocate(rPDYN)
      deallocate(rBdExg)
      deallocate(rCDOM)
      deallocate(rGPP)
      deallocate(rGPP_hourly)
      deallocate(rRDOC)
      deallocate(rLPOC)
      deallocate(rRLPOC)
     ! deallocate(rOCH4)
      deallocate(rOACH4)
     ! deallocate(rONCH4)
      deallocate(rPCH4)
      deallocate(rnPCO2)
      deallocate(rfloc)
      deallocate(Ipar)
      deallocate(pPAR)
      deallocate(rChl2C)
      deallocate(rC2Chl)
      deallocate(WtCLoad)
      deallocate(SwDOCLoad)
      deallocate(SwDICLOad)
      deallocate(SwPOCLoad)
      deallocate(SwSRPLoad)
      deallocate(SwDOLoad)
      deallocate(SwQso)
      deallocate(mem_sub%K1)
      deallocate(mem_sub%K2)
      deallocate(mem_sub%K3)
      deallocate(mem_sub%K4)
      deallocate(mem_sub%K5)
      deallocate(mem_sub%K6)
      deallocate(mem_sub%nxt4th)
      deallocate(mem_sub%nxt5th)
      deallocate(mem_sub%interim)
      deallocate(mem_sub%rerr)
      deallocate(mem_poc%K1)
      deallocate(mem_poc%K2)
      deallocate(mem_poc%K3)
      deallocate(mem_poc%K4)
      deallocate(mem_poc%K5)
      deallocate(mem_poc%K6)
      deallocate(mem_poc%nxt4th)
      deallocate(mem_poc%nxt5th)
      deallocate(mem_poc%interim)
      deallocate(mem_poc%rerr)
      deallocate(PC)
      deallocate(PCm)
   end subroutine

   subroutine CarbonModuleSetup(isHourNode,isDayNode)  ! added by Lin 20250903
      implicit none
      logical, intent(in) :: isHourNode,isDayNode !added by Lin 20250903
      real(r8) :: temp, w10, rho0, vv
      real(r8) :: k600SR, k600CC, k600UW
      real(r8) :: Schmidt ! Cpas, Cact  deleted by Lin
      integer  :: ii, gas, top
      real(r8) :: wind1, height      ! added by Lin

      height = 1.0_r8
      !if (lake_info%thrmkst==2 .and. lake_info%margin==1) then
      !   top = m_lakeWaterTopIndex 
      !   temp = T0
      !else
      !   top = 1
      !   temp = m_waterTemp(top) ! added according to ALBM-update
      !end if
      !top = 1 ! changed by Lin, as there is no yedoma on QTP
      top = m_lakeWaterTopIndex ! changed by Lin, 20250422
      temp = m_waterTemp(top) !
      wind1 = ConvertWindSpeed(m_surfData%wind, height) ! added by Lin to calculate wind speed at 1m height
      if (temp>=T0) then
         rho0 = m_wrho(top)
         vv = m_dVsc(top) / rho0
         k600CC = CalcPistonVelocity(m_surfData%wind)
         !k600SR = CalcPistonVelocity(m_surfData%wind, temp, rho0, vv, &
         !   m_Heff, m_Hmix, m_surfData%pressure)
         !k600UW = CalcPistonVelocity(lake_info, m_surfData%wind, temp, &
         !   rho0, m_Heff, m_Hmix)
         do gas = Wn2, Wch4, 1
            Schmidt = CalcSchmidtNumber(gas, temp)
            !Kg(gas) = k600SR / sqrt(Schmidt)
            if (wind1>3.7) then
                Kg(gas) = k600CC * (Schmidt/600.0)**(-0.5)
            else
                Kg(gas) = k600CC * (Schmidt/600.0)**(-0.67)
            end if
            ! Tan used to choose k600UW, but Lin has decided to use k600cc according to Wang et al., 2021
            ! According to Repo et al., 2007, changed the way on calculating Kg(gas) based on k600cc, which makes it dependent on wind speed at 1m.
            !Kg(gas) = k600UW / sqrt(Schmidt)
         end do
      else
         Kg = 0.0_r8
      end if
      ! update hourly chemistry fluxes 
      if (isHourNode) then
         call DeriveSubLateralFlux()
      end if
      call CalcCarbonCycleRates(isHourNode,isDayNode)
      if (isHourNode) then
         ! add burial C to sediment C pool
         !Cpas = 0.98*m_burialAlCarb   ! changed by Lin, so as follows
         !Cact = 0.02*m_burialAlCarb + m_burialAtCarb
         !m_unfrzCarbPool(:,1) = m_unfrzCarbPool(:,1) + &
         !   (/Cpas, Cact/)/m_dZs(1)
         m_unfrzCarbPool(1) = m_unfrzCarbPool(1) + &
             (m_burialAlCarb + m_burialAtCarb)/m_dZs(1)
         m_burialAlCarb = 0.0_r8
         m_burialAtCarb = 0.0_r8
      end if
      call UpdatePOMSettlingVelocity()
      if (m_Hice<e8) then
         !do ii = 1, NPOC, 1
         if (m_sinkPOCPool>e8) then
            RsPOCLoad = Scsed / lake_info%depth
         else
            RsPOCLoad = 0.0_r8
         end if
         !end do
      else
         RsPOCLoad = 0.0_r8
      end if
      call CalcSurfaceExchange(m_waterSubCon(:,top), top_fluxes)
   end subroutine

   subroutine CarbonModuleCallback(isHourNode, hindx, dt)
      implicit none
      logical, intent(in) :: isHourNode
      integer(i8), intent(in) :: hindx
      real(r8), intent(in) :: dt
      real(r8) :: avgPOC
      real(r8) :: maxPOC, cDOC, tzw
      real(r8) :: fdpoc, rDPOC, rsPOC
      integer, save :: prv_top = 1
      integer :: top, bottom, ii
 
      top = m_lakeWaterTopIndex
      bottom = WATER_LAYER + 1
      call AdjustNegativeConcentrations()
      call ConvectiveMixing()
      ! assuming no dissolved substances in ice layers
      ! coagulation precipitation
      ! m_burialAlCarb = 0.0_r8  ! deleted according to ALBM-update
      do ii = 1, top-1, 1
         cDOC = sum( m_waterSubCon(Waqdoc:Wtrdoc,ii) )
         m_burialAlCarb = m_burialAlCarb + cDOC * m_dZw(ii)
         m_waterSubCon(Waqdoc:Wtrdoc,ii) = 0.0_r8
         m_waterSubCon(Wo2,ii) = 0.0_r8
         m_waterSubCon(Wco2,ii) = 0.0_r8
         m_waterSubCon(Wch4,ii) = 0.0_r8
      end do
      if (m_Hice>e8) then
         do ii = top, m_mixTopIndex, 1
            cDOC = sum( m_waterSubCon(Waqdoc:Wtrdoc,ii) )
            m_burialAlCarb = m_burialAlCarb + cDOC*m_dZw(ii)
         end do
         m_waterSubCon(Waqdoc,top:m_mixTopIndex) = 0.0_r8
         m_waterSubCon(Wtrdoc,top:m_mixTopIndex) = 0.0_r8
        ! if (lake_info%hydroconn==1 .or. lake_info%depth>=10) then
        !    if (m_surfData%Qsi>5d3 .or. m_radPars%season==1) then
        !       m_waterSubCon(Wo2,top) = 3.5d5
        !    end if
        ! end if    ! closed by Lin in version20.6
      end if
      m_burialAlCarb = m_burialAlCarb + sum(rfloc*m_dZw)*dt
      !do ii = 1, NPOC, 1
      fdpoc = 1.0 - Fdom
      rDPOC = sum( (rLPOC(:)-rRLPOC(:))*fdpoc*m_dZw )
      m_burialAtCarb = m_burialAtCarb + rDPOC * dt
      !end do
      if (m_Hice<e8) then
         !do ii = 1, NPOC, 1
         rsPOC = m_waterPOC(bottom) * max(Vsettl(bottom),0.0) 
         m_sinkPOCPool = m_sinkPOCPool + (rsPOC-Scsed)*dt
         m_sinkPOCPool = max(0.0, m_sinkPOCPool)
         !end do
      end if
      ! redistribute phytoplankton in winter
      if (isHourNode) then
         if (prv_top/=top .and. top<WATER_LAYER) then
            if (prv_top<top) then
               avgPOC = 0.0_r8
               do ii = prv_top, top-1, 1
                  avgPOC = avgPOC + m_waterPOC(ii)*m_dZw(ii)
               end do
               avgPOC = avgPOC / m_dZw(top)
               maxPOC = 5d2 * m_chla0 / 0.24
               m_waterPOC(top) = m_waterPOC(top) + avgPOC
               !do ii = 1, NPOC, 1
               m_waterPOC(top) = min(m_waterPOC(top), maxPOC)
               !end do
            else
               tzw = m_dZw(prv_top) / sum( m_dZw(top:prv_top) )
               avgPOC = m_waterPOC(prv_top) * tzw
               do ii = top, prv_top-1, 1
                  m_waterPOC(ii) = avgPOC
               end do
            end if
            m_waterPOC(1:top-1) = 0.0_r8
         else if (top>=bottom) then
            m_waterPOC = 0.0_r8
         end if
         prv_top = top
      end if
      ! adjust chlorophyll for adaptation and inflow
      call AdjustAlgaeChla(isHourNode, hindx)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Discretize substance transport equation by Finite Difference Method.
   !          This sub-routine returns time derivatives of dissolved substance
   !          concentration.
   !
   !------------------------------------------------------------------------------
   subroutine CarbonCycleEquation(con, dcon)
      implicit none
      real(r8), intent(in) :: con(NWSUB,WATER_LAYER+1)      ! Unit: umol/m3
      real(r8), intent(out) :: dcon(NWSUB,WATER_LAYER+1)    ! Unit: umol/(m3*s)
      real(r8), dimension(NWSUB) :: qt, qb, dCa, dCb
      real(r8), dimension(NWSUB) :: Qso
      real(r8) :: aa, bb, as, bs
      integer :: ii, top, btm

      top = m_lakeWaterTopIndex
      btm = WATER_LAYER + 1
      if (top>WATER_LAYER-4) then
         dcon = 0.0_r8
         return
      end if
      call CalcBottomExchange(con(:,btm), qb)
      call CalcBoundaryExchange(qb, rBdExg)
      qt = top_fluxes
      dcon(:,1:top-1) = 0.0_r8
      do ii = top, WATER_LAYER+1, 1
         if(ii==top) then
            aa = 0.5 * (m_Kv(ii) + m_Kv(ii+1))
            bb = 1.0
            as = 0.5 * (m_Az(ii) + m_Az(ii+1))
            bs = m_Az(ii) 
            dca = (con(:,ii+1) - con(:,ii)) / (m_Zw(ii+1) - m_Zw(ii))
            dcb = qt
         else if(ii==WATER_LAYER+1) then
            aa = 1.0
            bb = 0.5 * (m_Kv(ii-1) + m_Kv(ii))
            as = m_Az(ii)
            bs = 0.5 * (m_Az(ii-1) + m_Az(ii))
            dca = qb
            dcb = (con(:,ii) - con(:,ii-1)) / (m_Zw(ii) - m_Zw(ii-1))
         else
            aa = 0.5 * (m_Kv(ii) + m_Kv(ii+1))
            bb = 0.5 * (m_Kv(ii) + m_Kv(ii-1))
            as = 0.5 * (m_Az(ii) + m_Az(ii+1))
            bs = 0.5 * (m_Az(ii) + m_Az(ii-1))
            dca = (con(:,ii+1) - con(:,ii)) / (m_Zw(ii+1) - m_Zw(ii))
            dcb = (con(:,ii) - con(:,ii-1)) / (m_Zw(ii) - m_Zw(ii-1))
         end if
         dcon(:,ii) = (aa*as*dca - bb*bs*dcb)/m_dZw(ii)/m_Az(ii) + &
            rDYN(:,ii) + rBdExg(:,ii)
         Qso = (/0d0, 0d0, SwQso(ii), 0d0, SwQso(ii), 0d0, SwQso(ii)/)
         dcon(:,ii) = dcon(:,ii) - con(:,ii)*Qso
      end do
      where(con<=0 .and. dcon<0) dcon = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Discretize substance transport equation by Finite Difference Method.
   !          This sub-routine returns time derivatives of non-dissolved 
   !          substance concentration (see Riley and Stefan, 1988; Ecologoical
   !          Modeling).
   !
   !------------------------------------------------------------------------------
   subroutine ParticulateEquation(con, dcon)
      implicit none
      real(r8), intent(in) :: con(WATER_LAYER+1)     ! Unit: umol/m3
      real(r8), intent(out) :: dcon(WATER_LAYER+1)   ! Unit: umol/(m3*s)
      real(r8) :: aa, bb
      integer :: ii, top

      top = m_lakeWaterTopIndex
      do ii = top, WATER_LAYER+1, 1
         if (ii==top) then
            aa = 0.0_r8 
            if (ii<WATER_LAYER+1) then
               bb = 0.5 * ( con(ii+1)*Vsettl(ii+1) + &
                  con(ii)*Vsettl(ii) ) !changed according to ALBM 2024 update
            else
               bb = 0.0_r8
            end if
         else if (ii==WATER_LAYER+1) then
            aa = 0.5 * ( con(ii-1)*Vsettl(ii-1) + &
               con(ii)*Vsettl(ii))
            bb = m_waterPOC(ii) * max(Vsettl(ii), 0.0) 
         else
            aa = 0.5 * con(ii-1) * Vsettl(ii-1)
            bb = 0.5 * con(ii+1) * Vsettl(ii+1)
         end if
         dcon(ii) = (aa - bb) / m_dZw(ii) + rPDYN(ii) - &
            con(ii)*SwQso(ii) + RsPOCLoad
      end do
      dcon(1:top-1) = 0.0_r8
      where(con<=0 .and. dcon<0) dcon = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Adjust negative dissolved gas concentrations (project negative 
   !          values to zero) ("Non-negative solutions of ODEs", Shampine, L., 
   !          2005).
   !
   !------------------------------------------------------------------------------
   subroutine AdjustNegativeConcentrations()
      implicit none

      where (m_waterSubCon<0) m_waterSubCon = 0.0_r8
      where (m_waterPOC<e8) m_waterPOC = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate carbon recycle reaction rates.
   !
   !          Important: for photosynthesis and CDOM photodegradation, the
   !          scalar irradiance should be used to calculate the rates.
   !
   !------------------------------------------------------------------------------
   subroutine CalcCarbonCycleRates(isHourNode, isDayNode)
      implicit none
      logical, intent(in) :: isHourNode, isDayNode
      real(r8) :: LPOC, Chla, rp
      real(r8) :: rDPOC, temp, SRP, Dco2
      real(r8) :: Do2, Dch4, aqDOC, trDOC
      integer :: ii, top, ltype

      top = m_lakeWaterTopIndex
      ltype = lake_info%itype
      rCDOM(1:top-1) = 0.0_r8
      rGPP(1:top-1) = 0.0_r8
      rLPOC(1:top-1) = 0.0_r8
      rRLPOC(1:top-1) = 0.0_r8
      rRDOC(:,1:top-1) = 0.0_r8
     ! rOCH4(1:top-1) = 0.0_r8
      rOACH4(1:top-1) = 0.0_r8
     ! rONCH4(1:top-1) = 0.0_r8
      rPCH4(1:top-1) = 0.0_r8
      rnPCO2(1:top-1) = 0.0_r8
      rfloc(1:top-1) = 0.0_r8
      do ii = top, WATER_LAYER+1, 1
         temp = max(m_waterTemp(ii), T0)
         !temp = m_waterTemp(ii) ! changed in 2025/4/16
         Dch4 = m_waterSubCon(Wch4,ii)
         Do2 = m_waterSubCon(Wo2,ii)
         Dco2 = m_waterSubCon(Wco2,ii)
         SRP = m_waterSubCon(Wsrp,ii)
         aqDOC = m_waterSubCon(Waqdoc,ii)
         trDOC = m_waterSubCon(Wtrdoc,ii)
         LPOC = m_waterPOC(ii)
         Chla = min( 1.55d2, 1d-3*rChl2C(ii)*LPOC ) ! changed according to ALBM-update
         m_chla(ii) = Chla
         rfloc(ii) = trDOC * floc
         if (isHourNode) then
            ! incident PAR radiation
            !if (m_Hsnow<e8) then
            call GetIncidentPAR(m_wvln, m_fsphot(:,ii), Ipar(ii), pPAR(ii))
            !else
            !   Ipar(ii) = 0.0_r8
            !end if

            ! update algae densities
            call UpdateAlgaeDensity(rChl2C(ii)/12.0, Ipar(ii), 3.6d3, RhoA(ii)) ! added rChl2C(:,ii)/12.0 according to ALBM-update

            ! photosynthesis, photochemical and microbial DOC degradation,
            ! algae metabolic loss, and CH4 oxidation
            call Photosynthesis(Chla, rChl2C(ii)/12.0, Dco2, temp, SRP, & ! added rChl2C(:,ii)/12.0 according to ALBM-update
                  Ipar(ii), rGPP(ii))
            rGPP_hourly(ii) = rGPP_hourly(ii) + rGPP(ii)        ! added by Lin 20250903
            if (isDayNode) then         ! added by Lin 20250903 after meeting
               rGPP(ii) = rGPP(ii) + rGPP_hourly(ii) * 0.08  ! including dark carbon sequestration
               rGPP_hourly(ii) = 0.0_r8
            end if
            !call Photosynthesis(Temp, pPAR(ii), rGPP(ii), rRLPOC(ii))
            !rp = rRLPOC(ii)
            call PhotoDegradation(ltype, trDOC, m_wvln, m_fsphot(:,ii), &
                  rCDOM(ii))
            call MicrobialRespiration((/aqDOC,trDOC/), temp, Do2, rRDOC(:,ii))
            call MetabolicLoss(LPOC, temp, rLPOC(ii))
            rRLPOC(ii)=0.02*exp(0.093*(temp-T0))*rGPP(ii)
            ! call POCDecomposition(DPOC, temp, Do2, rDPOC)
            !call Methanogenesis(lake_info%itype, aqDOC+trDOC, temp, Do2, &
            !   rPCH4(ii), rnPCO2(ii)) 
            call Methanotrophy(lake_info%itype, m_radPars%month, Dch4, Do2, temp, &
               m_Iab(ii), rOACH4(ii))
         end if

         if (Do2<minDo2) then
            rCDOM(ii) = 0.0_r8
            rRDOC(:,ii) = 0.0_r8
            rRLPOC(ii) = 0.0_r8
         end if
         
         if (Do2<minDo2 .or. Dch4<minDch4) then
            rOACH4(ii) = 0.0_r8
         end if
 
         if (Dco2<minDco2) then
            rGPP(ii) = 0.0_r8
         end if
      end do
      ! dynamics rates of dissolved substances
      ! 2.75*rRDOC(2,:) is a compensation for non-DOC oxidation
      rDYN(Wn2,:) = m_gasExchange(Wn2,:) + 0.72*rOACH4 ! Gu et al., 2021
      rDYN(Wo2,:) = rGPP - rCDOM - sum(rRDOC,1) - & ! changed according to ALBM-update
                    rRLPOC + SwDOLoad - 2.0*rOACH4 + &
                    m_gasExchange(Wo2,:)
      rDYN(Wco2,:) = -rGPP + rCDOM + sum(rRDOC,1) + & 
                     rRLPOC + rDIC2DOC * WtCLoad + &
                     SwDICLoad + rOACH4 + m_gasExchange(Wco2,:)! + rnPCO2(:)
      rDYN(Wch4,:) = -rOACH4 + m_gasExchange(Wch4,:)
      rDYN(Wsrp,:) = (rCDOM + sum(rRDOC,1))/YC2P_DOM + (rRLPOC - & ! changed according to ALBM-update
                     rGPP)/YC2P_POM + SwSRPLoad ! changed according to ALBM-update
      rDYN(Waqdoc,:) = -rRDOC(1,:) + Fdom*(rLPOC(:)-rRLPOC(:))
      rDYN(Wtrdoc,:) = -rCDOM - rRDOC(2,:) - rfloc + WtCLoad + SwDOCLoad
                        
      ! dynamics rates of particulate substances
      rPDYN = rGPP - rLPOC + SwPOCLoad
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update POM settling velocity (Stokes law)
   !          Phytoplankton motion in SBL follows MacIntyre (1998).
   !
   !------------------------------------------------------------------------------
   subroutine UpdatePOMSettlingVelocity()
      implicit none
      real(r8) :: dnw !dnp, dnm    deleted by Lin
      real(r8) :: vs2, Hcline      !vs1  deleted by Lin
      integer :: ii, top, bottom

      if (m_Hice>e8) then
         Vsettl = 0.0_r8
         return
      end if

      top = m_lakeWaterTopIndex
      bottom = WATER_LAYER + 1
      Hcline = 0.5 * (m_HbLayer(1) - m_HbLayer(2) + m_Zw(bottom))
      Vsettl(1:top-1) = 0.0_r8
      do ii = top, bottom, 1
         dnw = m_wrho(ii) 
         !dnp = RhoA(small_ppk,ii)
         !dnm = RhoA(large_ppk,ii)  !deleted by Lin
         !vs1 = CalcSettlingVelocity(dnw, dnp, daPico, m_dVsc(ii))
         !vs2 = CalcSettlingVelocity(dnw, dnm, daMicro, m_dVsc(ii))
         !vs1 = 9.838d-8   ! deleted by Lin; for small_ppk
         vs2 = 8.68d-6
         if (ii==top) then
            if (ii<m_mixTopIndex) then
               !vs1 = 0.0_r8
               vs2 = 0.0_r8
            else
               !vs1 = max(0.0_r8, vs1)
               vs2 = max(0.0_r8, vs2)
            end if
         else if (m_Zw(ii)<=Hcline) then
            !vs1 = 0.0_r8
            vs2 = 0.0_r8
         end if
         Vsettl(ii) = max( min(vs2, Vs0), -Vs0 ) 
      end do 
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate gas flux via diffusive transfer
   !      The surface gas concentration is calculated from the top 20 cm layers
   !      (Schubert et al., 2012).
   !      When the maximum snow depth is exceeded, the accumulated gas can be
   !      emitted to the atmosphere.
   !
   !------------------------------------------------------------------------------
   subroutine CalcSurfaceExchange(con, fluxes)
      implicit none
      real(r8), intent(in) :: con(NWSUB)        ! units: umol/m3
      real(r8), intent(out) :: fluxes(NWSUB)    ! units: umol/m2/s
      real(r8) :: Ceq, pH, temp, pressure
      real(r8) :: mixing_ratio(NGAS)
      real(r8) :: fehn(NGAS)
      integer :: top, gas

      if (lake_info%thrmkst==2 .and. lake_info%margin==1) then
         top = m_lakeWaterTopIndex
      else
         top = 1
      end if
     ! top = m_lakeWaterTopIndex ! Changed by Lin, 20250422
      if (top==1 .and. m_Hice>e8) then
         fluxes = 0.0_r8
      else
         if (m_radPars%month==7 .or. m_radPars%month==8) then
             pH = SMpH(lake_info%itype)
         else
             pH = LKpH(lake_info%itype)
         end if
         temp = m_waterTemp(top)
         ! maximum saturation of O2 is 50% [Holgerson et al., 2016]
         mixing_ratio = (/Xn2, Xo2, 1.0d-6*m_radPars%qCO2, 1.0d-6*m_radPars%qCH4/)
         fehn = 1.0_r8
         do gas = Wn2, Wch4, 1
            pressure = mixing_ratio(gas) * m_surfData%pressure
            Ceq = CalcEQConc(lake_info%itype, gas, temp, pH, pressure)
            fluxes(gas) = fehn(gas) * Kg(gas) * (con(gas) - Ceq)
         end do
         fluxes(Wsrp) = -AePLoad
         fluxes(Waqdoc) = 0.0_r8
         fluxes(Wtrdoc) = -RfCLoad - AeCLoad
      end if
   end subroutine
   
   subroutine CalcBottomExchange(con, fluxes)
      implicit none
      real(r8), intent(in) :: con(NWSUB)        ! units: umol/m3
      real(r8), intent(out) :: fluxes(NWSUB)    ! units: umol/m2/s
      real(r8) :: Porosity ! added according to ALBM-update
      integer :: ii

      Porosity = sa_params(Param_Por) ! added according to ALBM-update
      if (m_lakeWaterTopIndex>WATER_LAYER+1) then
         fluxes = 0.0_r8
      else
         ! bottom dissolved substance fluxes
         fluxes(Wn2:Wsrp) = m_Kbm / DeltaD * & ! changed according to ALBM-update
            (m_sedSubCon(:,1)/Porosity - con(Wn2:Wsrp)) ! changed according to ALBM-update
         fluxes(Wo2) = min(0.0, fluxes(Wo2))
         do ii = Wco2, Wsrp, 1
            fluxes(ii) = max(0.0, fluxes(ii))
         end do
         !fluxes(Waqdoc:Wtrdoc) = 0.0_r8
      end if
   end subroutine

   subroutine CalcBoundaryExchange(fluxes, exg)
      implicit none
      real(r8), intent(in) :: fluxes(:)   ! bottom fluxes (umol/m2/s)
      real(r8), intent(out) :: exg(:,:)   ! side exchange (umol/m3/s)
      real(r8) :: Vz
      integer :: ii, top

      top = m_lakeWaterTopIndex
      ! No side exchange at bottom and frozen zone
      exg(:,1:top-1) = 0.0_r8
      exg(:,WATER_LAYER+1) = 0.0_r8
      do ii = top, WATER_LAYER, 1
         Vz = m_Az(ii) * m_dZw(ii)
         exg(:,ii) = fluxes * m_dAz(ii) / Vz
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: DOC and POC inflow and outflow (Hanson et al., 2014 Limnol.
   !          Oceanogr.; Buffam et al., 2011 GCB; ...).
   !
   !------------------------------------------------------------------------------
   subroutine DeriveSubLateralFlux()
      implicit none
      real(r8) :: fcanopy, fwlnd, subflow, vol, units
      real(r8) :: DOCwt, Qsifrc, Qsofrc, pr
      integer  :: idxTop, idxBtm, ii

      pr = 2.0 * sqrt( PI*lake_info%Asurf )
      if (m_Hice<e8) then
         ! Precipitation carbon load (umol/m2/s)
         RfCLoad = DOCrf / MasC * 1.0d+6 * m_surfData%rainfall
         ! Aerial carbon load from canopy litter fall (umol/m2/s)
         fcanopy = 0.2 + 0.8 * m_ftree
         units = 1.0d+6 / MasC / SECOND_OF_DAY
         AeCLoad = fcanopy * DOCae * units * pr / lake_info%Asurf
      else
         RfCLoad = 0.0_r8
         AeCLoad = 0.0_r8
      end if

      ! Wetland carbon load from sub-surface flow (umol/m3/s)
      WtCLoad = 0.0_r8
      if (m_Hice<e8) then
         idxTop = m_mixTopIndex
         idxBtm = WATER_LAYER + 1
         vol = sum( m_Az(idxTop:idxBtm) * m_dZw(idxTop:idxBtm) )
         subflow = m_surfData%Qgw / vol
      else
         idxTop = 1
         idxBtm = WATER_LAYER + 1
         subflow = 0.0_r8
      end if
      fwlnd = 0.2 + 0.8 * m_fwlnd
      DOCwt = fwlnd * 1.0d6 * sa_params(Param_DOCwt) * subflow
      WtCLoad(idxTop:idxBtm) = DOCwt

      ! Stream-water flow carbon and phosphorus load (umol/m3/s)
      if (lake_info%hydroconn==0) then
         SwDOCLoad = 0.0_r8
         SwDICLoad = 0.0_r8
         SwSRPLoad = 0.0_r8
         SwPOCLoad = 0.0_r8
         SwDOLoad = 0.0_r8
         SwQso = 0.0_r8
         return
      end if

      call BinarySearch(m_wrho, m_surfData%dQsi, idxTop)
      idxBtm = WATER_LAYER + 1
      vol = sum( m_dZw(idxTop:idxBtm) * m_Az(idxTop:idxBtm) )
      if (m_Hice<e8) then
         Qsifrc = m_surfData%Qsi / vol
         Qsofrc = m_surfData%Qso / vol
      else
         Qsifrc = 0.0_r8
         Qsofrc = 0.0_r8
      end if
      do ii = 1, WATER_LAYER+1, 1
         if (ii>=idxTop .and. ii<=idxBtm) then
            SwDOCLoad(ii) = Qsifrc * m_surfData%DOCQsi
            SwDICLoad(ii) = Qsifrc * m_surfData%DICQsi
            SwSRPLoad(ii) = Qsifrc * m_surfData%SRPQsi
            SwQso(ii) = Qsofrc
         else
            SwDOCLoad(ii) = 0.0_r8
            SwDICLoad(ii) = 0.0_r8
            SwSRPLoad(ii) = 0.0_r8
            SwQso(ii) = 0.0_r8
         end if
      end do
      idxTop = m_lakeWaterTopIndex 
      idxBtm = m_mixTopIndex
      vol = sum( m_dZw(idxTop:idxBtm)*m_Az(idxTop:idxBtm) )
      if (m_Hice<e8) then
         Qsifrc = m_surfData%Qsi / vol
      else
         Qsifrc = 0.0_r8
      end if
      do ii = 1, WATER_LAYER+1, 1
         SwDOLoad(ii) = 0.0_r8 
         if (ii>=idxTop .and. ii<=idxBtm) then
            SwPOCLoad(ii) = 0.5 * Qsifrc * m_surfData%POCQsi 
         else
            SwPOCLoad(ii) = 0.0_r8
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Adjust algae chlorophyll hourly caused by adaptation and inflow.
   !
   !------------------------------------------------------------------------------
   subroutine AdjustAlgaeChla(isHourNode, hindx) 
      implicit none
      logical, intent(in) :: isHourNode
      integer(i8), intent(in) :: hindx
      real(r8) :: temp, SRP
      integer :: top, ii

      top = m_lakeWaterTopIndex
      ! adjustment due to adaptation
      if (isHourNode .and. mod(hindx-1,24)==12 .and. &
            (m_Hsnow<e8 .and. m_Hgrayice<e8) .and. top<=WATER_LAYER+1) then   ! changed according to ALBM 2024 update
         do ii=top, WATER_LAYER+1, 1
         !temp = max(m_waterTemp(top), T0)
         !SRP = m_waterSubCon(Wsrp,top)
            temp = m_waterTemp(ii)
            SRP = m_waterSubCon(Wsrp,ii)
            call CalcChl2CRatio(Ipar(ii), temp, SRP, rChl2C(ii), &
                rC2Chl(ii), PCm(ii), PC(ii)) ! changed from m_rChl2C to rChl2C according to ALBM-update         
         end do
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Get O2 and CO2 production rate.
   !
   !------------------------------------------------------------------------------
   subroutine GetO2ProductionRate(o2pro)
      implicit none
      real(r8), intent(out) :: o2pro   ! units: umol/m2/s
      integer :: ii

      o2pro = 0.0_r8
      !do ii = 1, NPOC, 1
      o2pro = o2pro + sum(rGPP(:)*m_dZw)
      !end do
   end subroutine

   subroutine GetCO2ProductionRate(pco2pro, bco2pro, rco2pro, rco2ch4)
      implicit none
      real(r8), intent(out) :: pco2pro    ! units: umol/m2/s
      real(r8), intent(out) :: bco2pro    ! units: umol/m2/s
      real(r8), intent(out) :: rco2pro    ! units: umol/m2/s
      real(r8), intent(out) :: rco2ch4    ! units: umol/m2/s
      integer :: ii

      pco2pro = sum( rCDOM * m_dZw ) 
      bco2pro = 0.0_r8
      do ii = 1, NDOC, 1
         bco2pro = bco2pro + sum(rRDOC(ii,:)*m_dZw)
      end do
      rco2pro = 0.0_r8
      !do ii = 1, NPOC, 1
      rco2pro = rco2pro + sum(rRLPOC(:)*m_dZw)
      !end do
     ! rco2ch4 = sum( rOACH4 * m_dZw )
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: pCO2 equivalent mixing ratio at 0.39m.
   !
   !------------------------------------------------------------------------------
   subroutine GetpCO2MixingRatio(mpCO2, mpCH4)
      implicit none
      real(r4), intent(out) :: mpCO2(WATER_LAYER+1)   ! units: ppm
      real(r4), intent(out) :: mpCH4(WATER_LAYER+1)   ! units: ppm
      real(r8) :: zz, ps, temp, henry
      real(r8) :: pCO2, pCH4, henry2, pH
      integer :: ii, indx

      do ii = 1, WATER_LAYER+1, 1
         ps = m_surfData%pressure + Roul*G*m_Zw(ii)
         if (m_waterIce(ii)>e8) then
            mpCO2(ii) = 0.0_r4
            mpCH4(ii) = 0.0_r4
         else
            temp = m_waterTemp(ii)
            if (m_radPars%month==7 .or. m_radPars%month==8) then
                pH = SMpH(lake_info%itype)
            else
                pH = LKpH(lake_info%itype)
            end if
            henry = CalcHenrySolubility(lake_info%itype, Wco2, temp, pH)
            henry2 = CalcHenrySolubility(lake_info%itype, Wch4, temp, pH)
            pCO2 = m_waterSubCon(Wco2,ii)
            pCH4 = m_waterSubCon(Wch4,ii)
            mpCO2(ii) = REAL( 1.0d6 * pCO2 / henry / ps )
            mpCH4(ii) = REAL( 1.0d6 * pCH4 / henry2 / ps )
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: convective mixing of dissolved substances.
   !
   !------------------------------------------------------------------------------
   subroutine ConvectiveMixing()
      implicit none
      real(r8) :: cavg(NWSUB)
      real(r8) :: tzw, dzw, Vz
      integer :: top, ii

      ! convective mixing
      top = m_lakeWaterTopIndex
      if (m_mixTopIndex>top) then
         cavg = 0.0_r8
         tzw = 0.0_r8
         do ii = top, m_mixTopIndex, 1
            Vz = m_dZw(ii) * m_Az(ii)
            tzw = tzw + Vz
            cavg = cavg + m_waterSubCon(:,ii)*Vz
         end do
         cavg = cavg / tzw
         do ii = top, m_mixTopIndex, 1
            m_waterSubCon(:,ii) = cavg
         end do
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Initialize water gas profiles.
   !
   !------------------------------------------------------------------------------
   subroutine InitializeCarbonStateVariables()
      implicit none
      real(r8) :: pGas, temp, pH
      real(r8) :: mixing_ratio(NGAS)
      integer :: gas, ii

      if (m_radPars%month==7 .or. m_radPars%month==8) then
          pH = SMpH(lake_info%itype)
      else
          pH = LKpH(lake_info%itype)
      end if
      mixing_ratio = (/Xn2, Xo2, Xco2, Xch4/)
      do ii = 1, WATER_LAYER+1, 1
         temp = m_waterTemp(ii)
         do gas = Wn2, Wch4, 1
            pGas = mixing_ratio(gas) * P0   ! gas partial pressure
            m_waterSubCon(gas,ii) = CalcEQConc(lake_info%itype, gas, temp, pH, pGas)
         end do
      end do
      !m_waterSubCon(Wch4,:) = 0.0_r8
      
      m_waterSubCon(Wsrp,:) = 1.0d3 * 0.03
      m_waterSubCon(Waqdoc,:) = 0.0_r8
      m_waterSubCon(Wtrdoc,:) = 1.0d6 * 0.38
      m_chla = 3._r8 ! changed according to ALBM-update
      if (lake_info%itype == 3.or.lake_info%itype == 4) then
         m_chla = 1._r8 ! added by Lin
      end if   
      m_waterPOC = 1d3 * m_chla / 0.24 ! changed from 5d2 * 1.0 / 0.24 according to ALBM-update

      m_gasExchange = 0.0_r8
      m_burialAtCarb = 0.0_r8
      m_burialAlCarb = 0.0_r8
      m_sinkPOCPool = 0.0_r8
   end subroutine

end module carbon_cycle_mod 
