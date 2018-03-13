!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! все функции инициализации
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer function CEA_INIT_in(Nof_,
     +     Nreac_,
     +     Fox_,
     +     Rname_,
     +     Symbols,
     +     Rtemp_,
     +     Pecwt_,
     +     Ent_)
      
      use CommonData
      IMPLICIT NONE
      INCLUDE 'cea.inc'
	integer i
      INTEGER Nof_,Nreac_
      CHARACTER*8 Fox_(100)
      CHARACTER*15 Rname_(100)
      character(15), intent(in) :: Symbols(10)
      REAL*8 Pecwt_(100),Rtemp_(100),Ent_(100)      
      integer(kind=4) :: ret
      character (100) fname
      
! ___________________________________________________________________________
	
!      cmData_dir_name=Data_dir_name
      CEA_INIT_in=-1
      
!     задание нач. значений  
      call SetInitFlags()
	call InitValues()     
      
!     открытие файлов баз данных       
      OPEN (IOSCH,STATUS='scratch',FORM='unformatted',
     *      iostat=ret)
      if (ret/=0) then
          call WrtError('scratch')
          ErrorFlag=10
          return 
      end if    
      fname=trim(cmData_dir_name)//'thermo.lib'
      OPEN (IOTHM,FILE=fname,status='OLD',
     *      FORM='unformatted',iostat=ret)     
      if (ret/=0) then
          call WrtError(fname)
          ErrorFlag=10
          return 
      end if    
      fname=trim(cmData_dir_name)//'trans.lib'
      OPEN (IOTRN,FILE=fname,status='OLD',
     *      FORM='unformatted',iostat=ret)
      if (ret/=0) then
          call WrtError(fname)
          ErrorFlag=10
          return 
      end if    

      if (ErrorFlag/=0) return

	call Input(Nof_,Nreac_,Fox_,Rname_,Symbols,Rtemp_,Pecwt_,Ent_)
      CALL SEARCH()
      CLOSE (IOSCH)
      CLOSE (IOTHM)
      CLOSE (IOTRN)
	Call Init2()
      CEA_INIT_in=i_hok
      end function CEA_INIT_in
      
!______________________________________________________________________________________________________________________
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
	subroutine Input(Nof_,Nreac_,Fox_,Rname_,Symbols,Rtemp_,Pecwt_,Ent_)
      IMPLICIT NONE
      INCLUDE 'cea.inc'
	integer i,j
	real*8 xyz,denmtr
      INTEGER Nof_,Nreac_
      CHARACTER*8 Fox_(MAXR)
      CHARACTER*15 Rname_(MAXR),Symbols(10)
      REAL*8  Pecwt_(MAXR),Rtemp_(MAXR)  ,Ent_(MAXR)
      character*4 RAtom_(12)

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
      RAtom_=""
      Rnum=0
      Enth=0
      Nfla=0
      do i=1,Nreac_
         RAtom_=""
         call SplitStr(Symbols(i),RAtom_)
         call GetAtomNumbers(RAtom_,Rnum(i,:))
         if (Rnum(i,1)>0) then
             Energy(i)=""
             Enth(i) = Ent_(i)*1.0D6/8314.51
             do j=1,10
              if (Rnum(i,j)==0) exit
              RAtom(i,j)=RAtom_(j)
              Nfla(i)=j               ! кол-во элементов в компоненте
             end do
          end if   
      end do   

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

! ___________________________________________________________________________
      subroutine SplitStr(str,RAtoms)
      implicit none
      character*15 str
      integer i,j,k,itmp
      character*4 RAtoms(12),tmp*3
      i=1
      j=0
      k=1
      if (str=="" .or. str(1:1)==" " .or. ichar(str(1:1))==0 ) return
      itmp=ichar(str(1:1))
!      itmp=str(1:1)      
      tmp=""
      do while(i<=len(str) .and. str(i:i)/=" ")
        j=j+1
        do while(i<=len(str) .and. str(i:i)/=";" .and. str(i:i)/=" "
     *       .and. ichar(str(i:i))/=0)
          itmp=ichar(str(i:i))           
             tmp(k:k)=str(i:i)
            k=k+1
            i=i+1
        enddo
        k=1
        RAtoms(j)=tmp
        tmp=""
        i=i+1
      end do
      end
 
    
      subroutine GetAtomNumbers(RAtoms,Rnum)
    ! split Atoms charaters to a name and a number
      implicit none
      character*4 RAtoms(12)
      REAL*8 Rnum(12)
      integer i,j,stat
      integer tmp

        do i=1,12
         if (RAtoms(i)=="" .or. RAtoms(i)(1:1)==" " ) then
              return   
         endif
         tmp=10;   
         do j=1,len(RAtoms(i))
             if (index("0123456789",RAtoms(i)(j:j))/=0) exit 
         end do
         if (j==13) return
         read(RAtoms(i)(j:),*,iostat=stat)  Rnum(i) 
         RAtoms(i)(j:)=""
        end do
      end
      