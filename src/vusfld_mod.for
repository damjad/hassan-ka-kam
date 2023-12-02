!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Subroutine VUSDFLD
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      subroutine VUSDFLD(
     + NBLOCK, NSTATEV, NFIELDV, NPROPS, NDIR, NSHR,
     + JELEM, KINTPT, KLAYER, KSECPT,
     + STEPTIME, TOTALTIME, DT, CMNAME,
     + COORDMP, DIRECT, T, CHARLENGTH, PROPS,
     + STATEOLD,STATENEW,FIELD)
      include 'vaba_param.inc'
!-----------------------------------------------------------------------
!-----Declaration ABAQUS variables
!-----------------------------------------------------------------------
      dimension JELEM(nblock),COORDMP(nblock,*),DIRECT(nblock,3,3),
     .          T(nblock,3,3),CHARLENGTH(nblock),PROPS(nprops),
     .          STATEOLD(nblock,nstatev),STATENEW(nblock,nstatev),
     .          FIELD(nblock,nfieldv)
      character*80 cmname
!-----Data from ABAQUS
      dimension stressdata(maxblk*(NDIR+NSHR)), straindata(maxblk*(NDIR+NSHR))
      integer jSData(maxblk*(NDIR+NSHR))
      character*3 cSData(maxblk*(NDIR+NSHR))
      integer jStatus
!-----------------------------------------------------------------------
!-----Declaration internal variables
!-----------------------------------------------------------------------
      integer i
      real*8 eps(NBLOCK,NDIR+NSHR)
      real*8 SMISES(NBLOCK),SIGH(NBLOCK),TRIAX(NBLOCK)
      real*8 si(NBLOCK)
      real*8 sig_l
      real*8 sig_cc
      real*8 eps_cc
      real*8 eps_l
      real*8 d_eps_c_d_eps_l
      real*8 d_sig_c_d_eps_l
      real*8 d_eps_p_l_d_eps_p_c

      real*8 sig_1c(NBLOCK)
      real*8 sig_2c(NBLOCK)
      real*8 sig_3c(NBLOCK)

      real*8 eps_1c(NBLOCK)
      real*8 eps_2c(NBLOCK)
      real*8 eps_3c(NBLOCK)

!----- Input
!----- TODO: assumptions
      real*8, parameter :: v_c = 0.2
      real*8, parameter :: E_c = 25742.960202742808
      real*8, parameter :: sig_cu = 30
      real*8, parameter :: f_co = 30
      real*8, parameter :: eps_cu = 0.001962081528017131
      real*8, parameter :: si_0 = 31
      real*8, parameter :: eps_c0 = 0.001962081528017131
!     todo: bhai dekh ly is ko
      real*8, parameter :: eps_c = 0.0015
      real*8, parameter :: eps_co = 0.001962081528017131


!-----------------------------------------------------------------------
!     Access stress and strain tensor
!-----------------------------------------------------------------------
      call vgetvrm( 'S', stressdata,jSData,cSData,jStatus)
      call vgetvrm( 'LE', straindata,jSData,cSData,jStatus)
!-----------------------------------------------------------------------
!     Extract data from straindata
!-----------------------------------------------------------------------
      if(nshr.gt.1)then
         do i=1,nblock
            sig(i,1) = stressdata(i)
            sig(i,2) = stressdata(i+nblock)
            sig(i,3) = stressdata(i+nblock*2)
            sig(i,4) = stressdata(i+nblock*3)
            sig(i,5) = stressdata(i+nblock*4)
            sig(i,6) = stressdata(i+nblock*5)
         enddo
      else
         do i=1,nblock
            sig(i,1) = stressdata(i)
            sig(i,2) = stressdata(i+nblock)
            sig(i,3) = stressdata(i+nblock*2)
!           sig_33 = sig(i, 4)
            sig(i,4) = stressdata(i+nblock*3)
         enddo
      endif
      if(nshr.gt.1)then
         do i=1,nblock
            eps(i,1) = straindata(i)
            eps(i,2) = straindata(i+nblock)
            eps(i,3) = straindata(i+nblock*2)
            eps(i,4) = straindata(i+nblock*3)
            eps(i,5) = straindata(i+nblock*4)
            eps(i,6) = straindata(i+nblock*5)
         enddo
      else
         do i=1,nblock
            eps(i,1) = straindata(i)
            eps(i,2) = straindata(i+nblock)
            eps(i,3) = straindata(i+nblock*2)
!           eps_33 = eps(i, 4)
            eps(i,4) = straindata(i+nblock*3)
         enddo
      endif
!-----------------------------------------------------------------------
!     Main Logic
!-----------------------------------------------------------------------
      do i=1,nblock
         if (-1 * eps(i, 4) < 0) then
            sig_l = 0
            si(i) = si_0
         else
!        Assumption: principal stress and strain in a direction is the direct stress in that direction
!        sig_11 = sig(i, 1)
            sig_1c(i) = sig(i, 1)
!        sig_22 = sig(i, 2)
            sig_2c(i) = sig(i, 2)
!        sig_33 = sig(i, 3)
            sig_3c(i) = sig(i, 3)

!        eps_11 = eps(i, 1)
            eps_1c(i) = eps(i, 1)
!        eps_22 = eps(i, 2)
            eps_2c(i) = eps(i, 2)
!        eps_33 = eps(i, 3)
            eps_3c(i) = eps(i, 3)

!-----------------------------------------------------------------------
!        Confinement pressure calculation
!-----------------------------------------------------------------------
            sig_leff = eq_16__sig_leff(sig_1c, sig_2c, sig_cu)
            FIELD(i, 2) = sig_leff

!-----------------------------------------------------------------------
!        Dilation angle calculation
!-----------------------------------------------------------------------
            if (-1 * eps_3c <= eps_cu) then
                si(i) = si_0
                FIELD(i, 1) = si(i)
            else
                eps_l = eq_33__eps_l(eps_1c(i), eps_2c(i))
                sig_l = eq_32__sig_l(sig_1c(i), sig_2c(i))
                sig_cc = eq_04__sig_cc(sig_cu, sig_l(i))
                eps_cc = eq_05__eps_cc(sig_l(i), sig_cu, eps_c0)
                r = eq_09__r(E_c, sig(i), eps_cc)
                sig_c(i) = eq_08__sig_c(sig_cc, eps(i), eps_c, eps_cc, r)

                d_eps_c_d_eps_l = eq_21__d_eps_c_d_eps_l(sig_l, f_co, eps_l, eps_c0)
                d_sig_c_d_eps_l = eq_19__d_sig_c_d_eps_l(r, sig_cc, eps_l, eps_cc, d_eps_c_d_eps_l)
                d_eps_p_l_d_eps_p_c = eq_30__d_eps_p_l_d_eps_p_c(v_c, E_c, d_sig_c_d_eps_l)

                si(i) = eq_31__si(d_eps_p_l_d_eps_p_c)
                FIELD(i, 1) = si(i)
            endif

         endif
      enddo
!-----------------------------------------------------------------------
!     Compute stress triaxiality
!-----------------------------------------------------------------------
      do i=1,nblock
         if(SMISES(i).gt.0.0)then
            TRIAX(i) = SIGH(i)/SMISES(i)
         else
            TRIAX(i) = 0.0
         endif
      enddo
!-----------------------------------------------------------------------
!     Update field variable
!-----------------------------------------------------------------------
      do i=1,nblock
         if(TRIAX(i).ge.0.0)then
            field(i,1) = 1.05
         else
            field(i,1) =-1.05
         endif
      enddo
!-----------------------------------------------------------------------
!     End of subroutine
!-----------------------------------------------------------------------
      return
      contains

      real*8 function eq_16__sig_leff(sig_1c, sig_2c, sig_cu)
          implicit none
          real*8, intent(in) :: sig_1c, sig_2c, sig_cu

          eq_16__sig_leff = (2.0d0 * (sig_1c + 0.039d0 * sig_cu) * (sig_2c + 0.039d0 * sig_cu) / (sig_1c + sig_2c + 0.078d0 * sig_cu)) - 0.039d0 * sig_cu

      end function eq_16__sig_leff

      real*8 function eq_33__eps_l(eps_1c, eps_2c)
          implicit none
          real*8, intent(in) :: eps_1c, eps_2c

          eq_33__eps_l = min(eps_1c, eps_2c)

      end function eq_33__eps_l

      real*8 function eq_32__sig_l(sig_1c, sig_2c)
          implicit none
          real*8, intent(in) :: sig_1c, sig_2c

          eq_32__sig_l = max(sig_1c, sig_2c)

      end function eq_32__sig_l

      real*8 function eq_04__sig_cc(sig_cu, sig_l)
          implicit none
          real*8, intent(in) :: sig_cu, sig_l

          eq_04__sig_cc = sig_cu * (1.0d0 + 3.5d0 * (sig_l / sig_cu))

      end function eq_04__sig_cc

      real*8 function eq_09__r(E_c, sig, eps_cc)
          implicit none
          real*8, intent(in) :: E_c, sig, eps_cc

          ! Check for division by zero
          if (eps_cc == 0.0d0) then
              print *, 'Division by zero encountered in eq_09__r'
              stop
          else
              eq_09__r = E_c / (E_c - (sig / eps_cc))
          endif

      end function eq_09__r

      real*8 function eq_08__sig_c(sig_cc, eps, eps_c, eps_cc, r)
          implicit none
          real*8, intent(in) :: sig_cc, eps, eps_c, eps_cc, r

          eq_08__sig_c = sig_cc * ((eps / eps_cc) * r) / (r - 1.0d0 + (eps_c / eps_cc)**r)
      end function eq_08__sig_c

      real*8 function eq_21__d_eps_c_d_eps_l(sig_l, f_co, eps_l, eps_co)
          implicit none
          real*8, intent(in) :: sig_l, f_co, eps_l, eps_co
          real*8 :: term1, term2, term3

          term1 = (0.85d0 + 6.8d0 * (sig_l / f_co))
          term2 = (-0.525d0 * (1.0d0 + 0.75d0 * (-eps_l / eps_co)**0.3d0))
          term3 = (-7.0d0 * exp(-7.0d0 * (-eps_l / eps_co)))

          eq_21__d_eps_c_d_eps_l = term1 * (term2 + term3)

      end function eq_21__d_eps_c_d_eps_l

      real*8 function eq_19__d_sig_c_d_eps_l(r, sig_cc, eps_l, eps_cc, d_eps_c_d_eps_l)
            implicit none
            real*8, intent(in) :: r, sig_cc, eps, eps_c, eps_cc, d_eps_c_d_eps_l

            eq_19__d_sig_c_d_eps_l = r * (sig_cc/eps_cc) * (r - 1.0d0) * (1.0d0 - (eps / eps_cc)**r) /(r - 1.0d0 + (eps_c / eps_cc)**r)**2 * d_eps_c_d_eps_l
      end function eq_19__d_sig_c_d_eps_l

      real*8 function eq_30__d_eps_p_l_d_eps_p_c(v, E_c, d_sig_c_d_eps_l)
            implicit none
            real*8, intent(in) :: v_c, E_c, d_sig_c_d_eps_l

            eq_30__d_eps_p_l_d_eps_p_c = (1.0d0 + v_c/E_c) * d_sig_c_d_eps_l / (d_sig_c_d_eps_l - (1.0d0 / E_c) * d_sig_c_d_eps_l)
      end function eq_30__d_eps_p_l_d_eps_p_c


      real*8 function eq_31__si(d_eps_p_l_d_eps_p_c)
            implicit none
            real*8, intent(in) :: d_eps_p_l_d_eps_p_c
            real*8 :: tan_si
            real*8, parameter :: rad_to_deg = 180.0d0 / acos(-1.0d0) ! Conversion factor from radians to degrees

            tan_si = -3.0d0 / 2.0d0 * (1.0d0 + 2.0d0 * d_eps_p_l_d_eps_p_c) / (1.0d0 - d_eps_p_l_d_eps_p_c)
            eq_31__si = atan(tan_si) * rad_to_deg
      end function eq_31__si

      end
