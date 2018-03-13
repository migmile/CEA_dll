	subroutine InitValues()
! Инициализация:      
!     задание нач. значений переменных      
      
      use CommonData      
      IMPLICIT NONE
      INCLUDE 'cea.inc'
	integer i

      Nonly = 0
      Nomit = 0
      Nsert = 0
      Trace = 0
      
!      DO i = 1,NCOL
!        Debug(i) = .FALSE.
!      ENDDO
      Nplt = 0
      Siunit = .TRUE.
      Moles = .FALSE.

      P=0
      V=0
      T=0
!          DO i = 1,MAXPV
!            P(i) = 0.
!            V(i) = 0.
!          ENDDO
!          DO i = 1,MAXT
!            T(i) = 0.
!          ENDDO
          P(1) = 1.
          Trace = 0.
          Lsave = 0
          R = Rr/4184.D0
          S0 = 0.

          Acat = 0.
          Ma = 0.
          Nfz = 1
          Nsub = 0
          Nsup = 0
          Npp = 0
          
          Pcp=0
          Supar=0
          Subar=0
          Mach1=0
          U1=0
!          DO i = 1,NCOL
!            Pcp(i) = 0.
!            Pcp(i+NCOL) = 0.
!            Supar(i) = 0.
!            Subar(i) = 0.
!            Mach1(i) = 0.
!            U1(i) = 0.
!          ENDDO

          Gamma1 = 0.
          Np = 0
		Nfz=1
          Nt = 1

	!	Hp = .TRUE.
		Pecwt = -1.
		Hsub0=1.e+30
		R = Rr/1000.
	
	    Iplt = 0
          Nplt = 0

          Nreac = 0
          Pecwt = -1.
!          DO i = 1,MAXR
!            Pecwt(i) = -1.
!          ENDDO
      end

!______________________________________________________________________________________________________________________
      subroutine SetInitFlags()
      IMPLICIT NONE
      INCLUDE 'cea.inc'
      
      Short = .TRUE.
      Massf = .TRUE.
      Moles = .FALSE.
      Debug = .FALSE.
      Tp = .FALSE.    
      Hp = .FALSE.
      Sp = .FALSE.
      Rkt = .FALSE.
      Shock = .FALSE.
      Detn = .FALSE.
      Vol = .FALSE.
      Ions = .FALSE.
      Eql = .FALSE.
      Froz = .FALSE.
      Fac = .FALSE.
      Debugf = .FALSE.
      Trnspt = .FALSE.
      Shkdbg = .FALSE.
      Incdeq = .FALSE.
      Incdfz = .FALSE.
      Refleq = .FALSE.
      Reflfz = .FALSE.
      
      end
      
!______________________________________________________________________________________________________________________
	subroutine Input(Nof_,Nreac_,Fox_,Rname_,Rtemp_,Pecwt_)
      IMPLICIT NONE
      INCLUDE 'cea.inc'
	integer i
	real*8 xyz,denmtr
      INTEGER Nof_,Nreac_
      CHARACTER*8 Fox_(MAXR)
      CHARACTER*15 Rname_(MAXR)
      REAL*8  Pecwt_(MAXR),Rtemp_(MAXR)      

      real eratio
      
!Исходные компоненты
      i_hok=1
      Nof=Nof_
      Nreac = Nreac_		! кол-во компонентов
      Fox=Fox_            !'fuel    '
      Rname=Rname_        !'Jet-A(g)       '
      Rtemp=Rtemp_        ! начальная температура
      Pecwt=Pecwt_        ! w=1. - доля компонента 
      Energy='lib'
      Nfla = 0
      Enth = 0.
      Jray = 0
      Rnum = 0.
      Rmw = 0.

! параметры расчета
	Nt=1		!кол-во точек по температуре
	Np=1		!кол-во точек по давлению

!______________________________________________________________________________________________________________________
      CALL REACT()

      DO i = 1,Nreac
          Jray(i) = 0
      ENDDO
      
	end subroutine Input

! ___________________________________________________________________________
	subroutine Init2()
      IMPLICIT NONE
      INCLUDE 'cea.inc'
      REAL*8 xi,xln
      REAL*8 DLOG
      INTEGER inc,j


C INITIAL ESTIMATES

        Npr = 0
        Gonly = .TRUE.
        Enn = .1D0
        Ennl = -2.3025851
        Sumn = Enn
        xi = Ng
        IF ( xi.EQ.0. ) xi = 1.
        xi = Enn/xi
        xln = DLOG(xi)

        DO inc = 1,Nc
          j = Ng + inc
          En(j,1) = 0.D0
          Enln(j) = 0.D0
        ENDDO

        DO j = 1,Ng
          En(j,1) = xi
          Enln(j) = xln
        ENDDO
	end
