      integer(4) function DllMain(hInst, ul_reason_being_called, 
     * lpReserved)
!DEC$ IF DEFINED(_X86_)
!DEC$ ATTRIBUTES STDCALL, ALIAS : '_DllMain@12' :: DllMain
!DEC$ ELSE
!DEC$ ATTRIBUTES STDCALL, ALIAS : 'DllMain' :: DllMain
!DEC$ ENDIF
      use DFWINA
      USE KERNEL32
      implicit none
      integer(4) hInst
      integer(4) ul_reason_being_called
      integer(4) lpReserved
      character(len=MAX_PATH) DLL_Name
      integer(HANDLE) hDLL
      integer n
      common /DLL_Handle/ hDLL
    
      hDLL = hInst

      ul_reason_being_called = ul_reason_being_called
      lpReserved = lpReserved

#if defined (_DEBUG_DLL_)
	open(26,file="c:\temp\ffffff")

	DLL_Name=""
      n = GetModuleFileName(hDLL, DLL_Name, sizeof(DLL_Name))
      if (n /= 0) then
          write(26,*)("DLL_main dll name:")
          write(26,*)(DLL_Name)
      else
          write(26,*)("can't get dll file name")
      end if   


      hDll = GetModuleHandle("NTDC_Dll.dll"//char(0))          ! почему-то дает другой хандл!!!!!
      n = GetModuleFileName(hDLL, DLL_Name, sizeof(DLL_Name))
      if (n /= 0) then
          write(26,*)("DLL_main dll name:")
          write(26,*)(DLL_Name)
      else
          write(26,*)("can't get dll file name")
      end if   


	write(26,*)"sdfgsdfgsdfgsdfg"
	close(26)
#endif

      DllMain = 1
      return

      end

!____________________________________________________________________      
      subroutine NTDC_INIT(log_path,names,symbols,params,reslt)
!     инициализация       
#if defined (_DLL_V_)
!DEC$ ATTRIBUTES DLLEXPORT,STDCALL,REFERENCE::NTDC_INIT
!DEC$ ATTRIBUTES  ALIAS:"ntdc_init" :: NTDC_INIT
!DEC$ ATTRIBUTES REFERENCE :: names
!DEC$ ATTRIBUTES REFERENCE :: symbols
!DEC$ ATTRIBUTES REFERENCE :: log_path
!DEC$ ATTRIBUTES REFERENCE :: params,reslt
#endif
      USE KERNEL32
	USE DFLIB 
      use CommonData
      implicit none
      real, intent(inout) :: params(100)
      integer  :: reslt   ! признак успешности завершения инициализации 
      character(MAX_PATH), intent(inout) :: log_path(1)      
      character(15), intent(inout) :: names(100),symbols(10)
      character(len=MAX_PATH) dir_name
  
      integer CEA_INIT_in

      real*8 Res(100),phi
      INTEGER Nof,Nreac
      CHARACTER*8 Fox(100)
      CHARACTER*15 Rname(100)
      REAL*8 Pecwt(100),Rtemp(100),Ent_(100)
      integer i
      CHARACTER*266 curdir

      integer(HANDLE) hDLL
      common /DLL_Handle/ hDLL

      character(len=256) EXE_name,DLL_name, cea_libs_path_name
      INTEGER DATE_TIME (8)
      CHARACTER (LEN = 12) REAL_CLOCK (3)
      integer(4) n
      logical Lexists

      dir_name='' 
      cmLog_dir_name=''      
      
      CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), 
     *                 REAL_CLOCK (3), DATE_TIME)
     
      ! обрубается запуск после 2017 года
      DATE_TIME(1)=DATE_TIME(1)-18
      if (DATE_TIME(1)>2000) then
          reslt=-1
          return
      endif 
      
      ! определим путь к файлу log      
      if (Len(log_path(1))>0) then
          cmLog_dir_name=trim(log_path(1))
          n=len(trim(cmLog_dir_name))+1
          cmLog_dir_name(n:n)='\'
      endif
      call WrtError(REAL_CLOCK(1))      
      call WrtError(REAL_CLOCK(2))      
      call WrtError(REAL_CLOCK(3))            

      
! Соберем статистику и определим пути
      IF (GetCurrentDirectory(LEN(curdir), curdir)>0) THEN 
          curdir = curdir(1 : INDEX(curdir, CHAR(0)))
          call WrtError("cur dir:")
          call WrtError(curdir)
       else
          call WrtError("can't get cur dir")
       end if
       
! name of exe
      n = GetModuleFileName(0, EXE_Name, sizeof(EXE_Name))
      if (n /= 0) then
          call WrtError("exe name:")
          call WrtError(EXE_Name)
      else
          call WrtError("can't get exe file name")
      end if    
       
!name of dll       
      hDll = GetModuleHandle("NTDC_Dll.dll"//char(0))          ! почему-то дает другой хандл!!!!!
      n = GetModuleFileName(hDLL, DLL_Name, sizeof(DLL_Name))
      if (n /= 0) then
          call WrtError("dll name:")
          call WrtError(DLL_Name)
      else
          call WrtError("can't get dll file name")
      end if    
      
      n=index(DLL_name,"\",.true.)
      dir_name=DLL_name(1:n)    
      
      inquire(FILE=trim(dir_name)//'thermo.lib',exist = Lexists )
      if (.NOT.Lexists) then
         dir_name=trim(curdir)
         n=LEN(trim(curdir))
         dir_name(n:n)='\'
         inquire(FILE=trim(dir_name)//'thermo.lib',exist=Lexists)
         if (.NOT.Lexists) then
            n=GETENVQQ("NTDC_DATA_PATH",cea_libs_path_name)
            if (n>0) then
             dir_name=cea_libs_path_name(1:n)
             call WrtError("env var NTDC_DATA_PATH:")
             call WrtError(dir_name)
            endif
         end if
      end if

      inquire(FILE=trim(dir_name)//'thermo.lib',exist=Lexists)
      if (.NOT.Lexists) then
          ErrorFlag=-100
          reslt=-100
          call WrtError("can't find file thermo.lib")
          return
      endif
      cmData_dir_name=dir_name

      call WrtError("lib dir:")
      call WrtError(dir_name)
      
      ErrorFlag=0
      call WrtError("____________________NTDC_INIT")
      open(unit=20,file=trim(cmLog_dir_name)//'NTDC_debug',
     *  access='APPEND')
       write(20,*) names(1),names(2),params
      close(20)

! собственно, начало инициализации
        Nof=1			        ! кол-во точек по er
        Nreac = params(1)		! кол-во компонентов
        Pecwt=-1
	  do i=1, Nreac
 	    if (params(1+(i-1)*4+1)==1) then
	    	Fox(i)='fuel    '
	    else  	    if (params(1+(i-1)*4+1)==0) then
	        Fox(i)='oxid    '
          else
              Fox(i)='name    '
	    end if
	    Rname(i) = names(i)
          Rtemp(i) = params(2+(i-1)*4+1)	  ! начальная температура
          Pecwt(i) = params(3+(i-1)*4+1)    !  wt=1. - доля компонента 
          Ent_(i) = params(4+(i-1)*4+1)     !  энтальпия 
	  end do

        reslt= CEA_INIT_in(Nof,Nreac,Fox,Rname,symbols,Rtemp,Pecwt,Ent_)
	  call WrtError("end NTDC_INIT")

      end subroutine NTDC_INIT

!_______________________________________________________________________________
      subroutine NTDC_HP(phi,Res)
! расчет параметров при заданном давлении и постоянной энтальпии      
#if defined (_DLL_V_)      
!DEC$ ATTRIBUTES DLLEXPORT,STDCALL,REFERENCE::NTDC_HP
!DEC$ ATTRIBUTES REFERENCE ::Res
!DEC$ ATTRIBUTES REFERENCE ::phi
#endif
      use CommonData
        IMPLICIT NONE
        interface
        subroutine CEA_calculate(c_t,phi,Res)
      	real*8 Res(100),phi      
          integer c_t
        end subroutine CEA_calculate
        end interface
        
        real*8 Res(100),phi
 
        open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access='APPEND')
        write(20,*) phi,Res
        close(20)

        call CEA_calculate(0,phi,Res)

      end subroutine NTDC_HP

!_______________________________________________________________________________
      subroutine NTDC_TP(phi,Res)
! расчет параметров при заданном давлении и температуре            
#if defined (_DLL_V_)      
!DEC$ ATTRIBUTES DLLEXPORT,STDCALL, REFERENCE::NTDC_TP
!DEC$ ATTRIBUTES REFERENCE ::Res
!DEC$ ATTRIBUTES REFERENCE ::phi
#endif
      use CommonData
        IMPLICIT NONE

        interface
        subroutine CEA_calculate(c_t,phi,Res)
      	real*8 Res(100),phi      
          integer c_t
        end subroutine CEA_calculate
        end interface
                
        real*8 Res(100),phi
 
        open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access="Append")
        write(20,*) phi,Res
        close(20)

        call CEA_calculate(1,phi,Res)

      end subroutine NTDC_TP

! ___________________________________________________________________________      
	
#if defined (_DLL_V_)      
      subroutine MixFractions (strings,x_) 
!DEC$ ATTRIBUTES DLLEXPORT,STDCALL, REFERENCE::MixFractions
!DEC$ ATTRIBUTES REFERENCE :: strings
!DEC$ ATTRIBUTES REFERENCE :: x_
	real, intent(inout) :: x_(100)
	character(15), intent(inout) :: strings(100)

#elif defined (_EXE_Intel_)
      subroutine MixFractions (strings,x_) 
      implicit none
      real, intent(inout) :: x_(100)
	character(15), intent(inout) :: strings(100)
#else
!      subroutine MixFractions (strings,x_) bind(C,name="MixFractions")
 !     use ISO_C_BINDING
 !     implicit none
  !    real, intent(inout) :: x_(100)
!	character(15), intent(inout) :: strings(100)
#endif

      integer n
      call CalculateFractions(strings,x_,n)
      end subroutine MixFractions
    
      
! ___________________________________________________________________________    
      subroutine NTDC_ROCKET(phi,Acat_,Supar_,res,strings,x_)
#if defined (_DLL_V_)      
!DEC$ ATTRIBUTES DLLEXPORT,STDCALL, REFERENCE::NTDC_ROCKET
!DEC$ ATTRIBUTES REFERENCE :: phi
!DEC$ ATTRIBUTES REFERENCE :: res,Acat_,Supar_  
!DEC$ ATTRIBUTES REFERENCE :: strings
!DEC$ ATTRIBUTES REFERENCE :: x_
#endif
      implicit none
       interface
        subroutine Rocket1(phi,Acat_,Supar_,Res,strings,x_)
      	real*8 Res(*),Acat_,Supar_,phi 
          real, intent(inout) :: x_(100,10)
	    character(15), intent(inout) :: strings(100)          
        end subroutine Rocket1
       end interface
       
       real*8 Res(200),Acat_,Supar_,phi 
       real, intent(inout) :: x_(100,10)
	 character(15), intent(inout) :: strings(100)
      
       call Rocket1(phi,Acat_,Supar_,res,strings,x_)
      end subroutine NTDC_ROCKET

! ___________________________________________________________________________    
      subroutine NTDC_IE_T(phi,T1,P1,T2,Res,IsOk)   !CEA_Isoentr_T
#if defined (_DLL_V_)      
! Подбирает давление под заданные Т и S      
!DEC$ ATTRIBUTES DLLEXPORT,STDCALL, REFERENCE::NTDC_IE_T
!DEC$ ATTRIBUTES  ALIAS:"NTDC_IE_T" :: NTDC_IE_T
!DEC$ ATTRIBUTES REFERENCE ::Res
!DEC$ ATTRIBUTES REFERENCE ::phi
!DEC$ ATTRIBUTES REFERENCE ::T1
!DEC$ ATTRIBUTES REFERENCE ::P1
!DEC$ ATTRIBUTES REFERENCE ::T2
!DEC$ ATTRIBUTES REFERENCE ::IsOk
#endif
      use CommonData
        IMPLICIT NONE

        real*8::  Res(100),phi,T1,P1,T2
        integer   IsOk
        
        real*8::  EPS=1e-4
        real*8::  Ptmp,S1,k1,S2,dS,R
        integer:: i
 
        open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access="Append")
        write(20,*) "NTDC_Isoentr_T"
        write(20,*) phi,Res
        close(20)
! TP - initial S&k        
      Res(1)=T1
      Res(2)=P1      
      call CEA_calculate(1,phi,Res)       
      S1=Res(6)
      k1=Res(9)
      
      k1=1.4
      Ptmp=P1*(T2/T1)**(k1/(k1-1))    ! first assumption
      ! find p for state S=const
      S2=S1      
      do i=1,500
            Res(1)=T2
            Res(2)=Ptmp      
            call CEA_calculate(1,phi,Res)       ! TP 
            dS=S1-Res(6)
            if (abs(dS)<EPS) exit
            k1=Res(9)
            R=Res(8)*(1-1/k1)
            Ptmp=Ptmp*exp(-ds/R)
      end do
      
      if (i>500) then
        open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access="Append")
         write(20,*) "NTDC_Isoentr_T - ERROR!"
         write(20,*) "STOPED after 10 steps without desision!"
        close(20)
        IsOk=0
        return
      end if          
      IsOk=1
      end subroutine NTDC_IE_T
      

! ___________________________________________________________________________    
      subroutine NTDC_Stat_Total(phi,T,P,M,Res,IsOk)   
!   идея взята у Д.Давиденко
! Numerical method:
!         Iterative process is intended to satisfy energy conservation. It
!         is devided into two stages: streight forward iteration, when total
!         pressure is approximated using gasdynamic function; iteration
!         between limits provided by the first stage, when enthalpy is
!         linearized versus temperature and pressure is determined using
!         adiabatic relation. Entropy conservation condition is provided
!         by 'ISENTA' subroutine.
!----------------------------------------------------------------------------      
!   M->H2->T2->P2            
#if defined (_DLL_V_)      
!DEC$ ATTRIBUTES DLLEXPORT,STDCALL, REFERENCE::NTDC_Stat_Total
!DEC$ ATTRIBUTES  ALIAS:"NTDC_Stat_Total" :: NTDC_Stat_Total
!DEC$ ATTRIBUTES REFERENCE ::Res
!DEC$ ATTRIBUTES REFERENCE ::phi
!DEC$ ATTRIBUTES REFERENCE ::T
!DEC$ ATTRIBUTES REFERENCE ::P
!DEC$ ATTRIBUTES REFERENCE ::M
!DEC$ ATTRIBUTES REFERENCE ::IsOk
#endif
      use CommonData
        IMPLICIT NONE

        real*8::  Res(100),phi,T,P,M
        integer   IsOk
        
        real*8::  EPS=1e-6
        real*8::  Ptmp,S1,k1,S2,T1,P1,M1,dS,T2,H1,H2,dh,R,Hf,Ht2,Tt2
        integer:: i
 
        open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access="Append")
        write(20,*) "NTDC_Stat_Total"
        write(20,*) phi,Res
        close(20)
        
! TP - initial S&k  
        T1=T
        P1=P
        M1=M
      Res(1)=T1
      Res(2)=P1      
      call CEA_calculate(1,phi,Res) 
      H1=Res(4)                           ! начальное значение
      Hf=H1+(M1*Res(10))**2/2              ! конечное значение 
      k1=Res(9)
      Ptmp=P1 
      do i=1,50
            T2=T1+(M1*Res(10))**2/2/Res(8)
            Res(1)=T2
            Res(2)=Ptmp      
            call NTDC_IE_T(phi,T1,P1,T2,Res,IsOk)  
            if (IsOk==0) return                      
            H2=Res(4)
            if (H2>Hf) exit
            R=Res(8)*(1-1/Res(9))
            M1=sqrt(2.0*(Hf-H2)/Res(9)/T2/R)
            T1=T2
            H1=H2
      enddo
    
      if (i>50) then
        open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access="Append")
        write(20,*) "NTDC_Isoentr_T - ERROR!"
        write(20,*) "STOPED after 10 steps without desision!"
        close(20)
        IsOk=0
        return
      end if
      
      Ptmp=Res(2)      
      do i=1,50
            T2=T1+(T2-T1)*(Hf-H1)/(H2-H1)
            Res(1)=T2
            Res(2)=Ptmp      
            call NTDC_IE_T(phi,T1,P1,T2,Res,IsOk)  
            if (IsOk==0) return                      
            dH=Hf-Res(4)
            if (abs(dH)<EPS) exit
            Ptmp=Res(2)
            H2=Res(4)
      end do 
      
      if (i>50) then
        open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access="Append")
        write(20,*) "NTDC_Stat_Total - ERROR!"
        write(20,*) "STOPED after 50 steps without desision!"
        close(20)
        IsOk=0
        return
      end if
      IsOk=1
      end subroutine NTDC_Stat_Total
      

! ___________________________________________________________________________    
      subroutine NTDC_Total_Stat(phi,T,P,M,Res,IsOk)   
      !   идея взята у Д.Давиденко
      !   M->H2->T2->P2       
#if defined (_DLL_V_)      
!DEC$ ATTRIBUTES DLLEXPORT,STDCALL, REFERENCE::NTDC_Total_Stat
!DEC$ ATTRIBUTES  ALIAS:"NTDC_Total_Stat" :: NTDC_Total_Stat
!DEC$ ATTRIBUTES REFERENCE ::Res
!DEC$ ATTRIBUTES REFERENCE ::phi
!DEC$ ATTRIBUTES REFERENCE ::T
!DEC$ ATTRIBUTES REFERENCE ::P
!DEC$ ATTRIBUTES REFERENCE ::M
!DEC$ ATTRIBUTES REFERENCE ::IsOk
#endif
      use CommonData
        IMPLICIT NONE

        real*8::  Res(100),phi,T,P,M
        integer   IsOk
        
        real*8::  EPS=1e-6,Cp,a
        real*8::  Ptmp,S1,k1,S2,T1,P1,M1,dS,T2,H1,H2,dh,R,Hf,Ht2
        real*8::  Tf,M2,k,Ht
        integer:: i
 
        open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access="Append")
        write(20,*) "NTDC_Stat_Total"
        write(20,*) phi,Res
        close(20)
        
! TP - initial S&k  
      T1=T
      P1=P
      M1=M
      Res(1)=T1
      Res(2)=P1      
      call CEA_calculate(1,phi,Res) 
      Ht=Res(4)                           ! начальное значение
      H1=Ht
      k=1.4d0
      a=300.d0
      Cp=1005.d0
      Ptmp=P1 
!      do i=1,50
          ! вторая точка
            T2=T1-(M1*a)**2/2.d0/Cp
            T2=T1/(1.d0+(k-1)*M1**2/2.d0)            
            T2=max(50.d0,T2)                 ! ограничение          
            call NTDC_IE_T(phi,T1,P1,T2,Res,IsOk) 
            if (IsOk==0) return                      
            H2=Res(4)
            M2=sqrt(2.*(H1-H2))/Res(10)
            a=Res(10)
            Cp=Res(8)
            k=Res(9)
!            if (abs(M2-M)/M<0.1) exit
!      enddo
!      if (i>50) then
!        open(20,file="cea_debug",access="Append")
!        write(20,*) "CEA_Total_Stat - ERROR!"
!        write(20,*) "STOPED after 50 steps without desision!"
!        close(20)
!        IsOk=0
!        return
!      end if
 
      ! первая точка
      T1=T1/(1.d0+(k-1)*M2**2/2.)
      T1=max(50.d0,T1)                 ! ограничение          
      call NTDC_IE_T(phi,T,P,T1,Res,IsOk) 
      if (IsOk==0) return                
      H1=Res(4)
      M1=sqrt(2.*(Ht-H1))/Res(10)
      P1=Res(2)      
      
      do i=1,50
          Tf=T1+(T2-T1)/(M2-M1)*(M-M1)
          Tf=max(10.,Tf)
          call NTDC_IE_T(phi,T1,P1,Tf,Res,IsOk) 
          if (IsOk==0) return          
          H2=Res(4)
          M2=sqrt(2.*(Ht-H2))/Res(10)
          if (abs(M2-M)/M<EPS) exit
          T2=Tf
          H2=Res(4)
      enddo
      if (i>50) then
        open(20,file=trim(cmLog_dir_name)//'NTDC_debug',access="Append")
        write(20,*) "NTDC_Total_Stat - ERROR!"
        write(20,*) "STOPED after 50 steps without desision!"
        close(20)
        IsOk=0
        return
      end if
      IsOk=1
      end subroutine NTDC_Total_Stat

