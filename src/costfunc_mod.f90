module costfunc_mod
!---------------------------------------------------------------------------------
! Purpose: this module provides cost functions for calibration and sensitivity test
!
!---------------------------------------------------------------------------------
   use phy_utilities_mod
   use read_data_mod
   use data_buffer_mod
   use math_utilities_mod

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Get cost function returns for calibration test
   !
   !------------------------------------------------------------------------------
   subroutine EvaluateSampleLikelihood(time, calib_vars, odata)
      implicit none
      type(SimTime), intent(in) :: time
      character(len=8), intent(in) :: calib_vars(:)
      real(r8), intent(out) :: odata(:)
      character(cx) :: fullname
      character(len=8) :: str
      real(r8), allocatable :: tmpArr(:)
      real(r8), allocatable :: tmpAnArr(:)
      real(r8), allocatable :: tmpArr2(:,:)
      real(r8), allocatable :: tmpfc(:)
      real(r4), allocatable :: surftemp(:) ! added by Lin
      integer, allocatable :: yearHist(:)
      integer, allocatable :: montHist(:) ! added by Lin
      integer, parameter :: ntd_w = 4     ! for weekly average
      integer, parameter :: ntd_d = 1     ! for daily average
      integer :: nvar, ii, jj, kk
      integer :: nyear, JDN0, JDN1, JDNa, JDNb
      integer :: nt, nt0, nt1, nzw, nzs ! nz is divided into nzw and nzs by Lin
      real(r8) :: fcost, std
      real(r4) :: T0_r4

      nt = size(m_timeHist)
      nzw = size(m_Zw) ! changed by Lin
      nzs = size(m_Zs) ! added by Lin
      nyear = time%year1 - time%year0 - 1 
      allocate(tmpArr(nt))
      allocate(tmpAnArr(nyear))
      allocate(yearHist(nyear))
      allocate(montHist(12)) ! added by Lin
      allocate(tmpfc(12)) ! added by Lin
      allocate(tmpArr2(nzw,nt)) ! changed by Lin
      allocate(surftemp(nt)) ! added by Lin
      do ii = 1, nt, 1 ! added by Lin
         surftemp(ii) = m_tempsHist(1,ii)
      end do
      
      do ii = 1, nyear, 1
         yearHist(ii) = time%year0 + ii
      end do

      do ii = 1, 12, 1
         montHist(ii) = ii
      end do   
      nvar = size(calib_vars) 
      odata = 0.0_r8
      do ii = 1, nvar, 1
         str = calib_vars(ii)
         if (len_trim(lakeid_file)==0) then
            write(fullname, "(A,I0,A)") trim(obs_dir) // trim(str) // &
               '_obs_', lake_info%id, '.dat'
         else
            fullname = trim(obs_dir) // trim(str) // '_obs_' // &
               trim(lake_info%name) // '.dat'
         end if 

         fcost = -9999.0_r8
         ! depth usually measured from the bottom of the ice layer
         if (trim(str)=='tw') then

            ! water temperature (K)
            !print *,'str=tw'
            tmpArr2 = -9999.0_r8
            do jj = 1, nt, 1
               tmpArr2(:,jj) = m_Zw - m_iceHist(jj)
            end do
            std = 1d-1
            T0_r4 = REAL(T0)
            call CalcCostfunc4Var(fullname, m_timeHist, ntd_d, &
               tmpArr2, m_tempwHist, std, fcost)

         else if (trim(str)=='ts')then ! added by Lin
            
            ! sediment temeprature(celcius)
            !print *,'str=ts'
                    
            std = 1d-1
            T0_r4 = REAL(T0)
            call CalcCostfunc4Var(fullname, m_timeHist, ntd_d, &
               surftemp-T0_r4, std, fcost)

         else if (trim(str)=='ice') then

            ! ice thickness (m) 
            !print *,'str=ice'
            std = 1d-2
            call CalcCostfunc4Var(fullname, m_timeHist, ntd_w, &
               m_iceHist, std, fcost)

         else if (trim(str)=='snow') then

            ! snow thickness (m) 
            std = 1d-2
            call CalcCostfunc4Var(fullname, m_timeHist, ntd_w, &
               m_snowHist, std, fcost)

         else if (trim(str)=='iceon') then

            ! ice on (DOY) 
            call Date2JDN(time%year0, time%month0, time%day0, JDN0)
            do jj = 1, nyear, 1
               call Date2JDN(time%year0+jj, 1, 1, JDN1)
               nt0 = 24*(JDN1-JDN0) + 1
               if (IsLeapYear(time%year0+jj)) then
                  nt1 = nt0 + 8783
               else
                  nt1 = nt0 + 8759
               end if
               do kk = nt1, nt0, -1
                  if (m_iceHist(kk)<=0) then
                     tmpAnArr(jj) = (kk-nt0+1.0) / 24.0
                     exit
                  end if
               end do
            end do
            std = 1d0
            call CalcCostfunc4AnnVar(fullname, yearHist, tmpAnArr, &
               std, fcost)

         else if (trim(str)=='iceoff') then

            ! ice off (DOY) 
            call Date2JDN(time%year0, time%month0, time%day0, JDN0)
            do jj = 1, nyear, 1
               call Date2JDN(time%year0+jj, 1, 1, JDN1)
               nt0 = 24*(JDN1-JDN0) + 1
               if (IsLeapYear(time%year0+jj)) then
                  nt1 = nt0 + 8783
               else
                  nt1 = nt0 + 8759
               end if
               tmpAnArr(jj) = minloc(m_iceHist(nt0:nt1), 1) / 24.0
            end do
            std = 1d0
            call CalcCostfunc4AnnVar(fullname, yearHist, tmpAnArr, &
               std, fcost)

         else if (trim(str)=='fco2') then
            call Date2JDN(time%year0, time%month0, time%day0, JDN0)
            do jj = 1, 12, 1
               call Date2JDN(2020, jj, 1, JDNa)
               call Date2JDN(2020, jj+1, 1, JDNb)
               nt0 = 24*(JDNa-JDN0) + 1
               nt1 = 24*(JDNb-JDN0) + 1
               call Mean(DBLE(m_fco2Hist(nt0:nt1)), tmpfc(jj))
               tmpfc(jj) = tmpfc(jj)/24.0
            end do   
            std = 1d-1
            call CalcCostfunc4AnnVar(fullname, montHist, tmpfc, &
               std, fcost)
         
         else if (trim(str)=='fch4d') then
            call Date2JDN(time%year0, time%month0, time%day0, JDN0)
            do jj = 1, 12, 1
               call Date2JDN(2020, jj, 1, JDNa)
               call Date2JDN(2020, jj+1, 1, JDNb)
               nt0 = 24*(JDNa-JDN0) + 1
               nt1 = 24*(JDNb-JDN0) + 1
               call Mean(DBLE(m_fch4dHist(nt0:nt1)), tmpfc(jj))
               tmpfc(jj) = tmpfc(jj)/24.0
            end do
            std = 1d-1
            call CalcCostfunc4AnnVar(fullname, montHist, tmpfc, &
               std, fcost)

         else if (trim(str)=='fch4e') then
            call Date2JDN(time%year0, time%month0, time%day0, JDN0)
            call Date2JDN(time%year1, time%month1, time%day1, JDN1)
            do jj = 1, 12, 1
               call Date2JDN(2020, jj, 1, JDNa)
               call Date2JDN(2020, jj+1, 1, JDNb)
               nt0 = 24*(JDNa-JDN0) + 1
               nt1 = 24*(JDNb-JDN0) + 1
               call Mean(DBLE(m_fch4eHist(nt0:nt1)), tmpfc(jj))
               tmpfc(jj) = tmpfc(jj)/24.0
            end do
            std = 1d-1
            call CalcCostfunc4AnnVar(fullname, montHist, tmpfc, &
               std, fcost)

         end if
         if (fcost<0) then
            odata(ii) = -fcost
         else
            odata(ii) = -0.5 * fcost
         end if
      end do

      deallocate(tmpArr)
      deallocate(tmpArr2)
      deallocate(tmpAnArr)
      deallocate(tmpfc)
      deallocate(yearHist)
      deallocate(montHist)
      deallocate(surftemp)
   end subroutine

end module costfunc_mod
