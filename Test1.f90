    program Test1
     USE KERNEL32
	 USE DFLIB 
     implicit none    
  
      character(15) :: symbols(10)
    ! Variables
     real :: params(100)
     integer  :: reslt   ! признак успешности завершения инициализации 
     character(MAX_PATH) :: log_path(1)      
     character(15) :: names(100)
     real*8::Res(200),phi,T,P,M
     integer i
     real *8 ratios(4)
	 real :: x_(100)


     
     params(1)=2        ! NReac     кол-во компонентов
   
     names(1)="H2" 
     params(2)=1        ! 0-окисл/ 1-топл
     params(3)=293      ! T0     
     params(4)=1    ! 0.768        ! amount

     names(2)="O2"          !Jet-A(L)     
     params(6)=0        ! 0-окисл/ 1-топл
     params(7)=293     ! T0     
     params(8)=1.00    !-0.768  ! amount 1-0.768
     

!     names(3)="Jet-A(L)"                
!     params(8)=1       ! 0-окисл/ 1-топл
!     params(9)=283     ! T0     
!    params(10)=1  ! amount 1-0.768

!     names(3)="O2"           
!     params(8)=0        ! 0-окисл/ 1-топл
!     params(9)=283      ! T0     
!     params(10)=1  ! amount 1-0.768
     
     
!     params(8)=-10000       ! pressure  ???
    log_path(1)=""!c:\cea_dll"
!    call CEA_INIT(log_path,names,params,reslt)
     call    NTDC_INIT(log_path,names,symbols,params,reslt)    
    if (reslt<0) then
        stop "Error INIT"
    end if    

    P=10
    Res(2)=P      
    phi=0.5
    call NTDC_HP(phi,Res)
    call MixFractions (names,x_) 


    phi=1
!    call CEA_ROCKET(phi,res)
  
!    stop
    
!    T=2300
!    P=4d
!    Res(1)=T
    Res(2)=P      
!    phi=0
!    call CEA_HP(phi,Res)
! call MixFractions (names,x_) 
    
    write(*,*) 'P1=',Res(2),'T1=',Res(1),'H=',Res(4),'S=',Res(6)

    
    T=2300
    P=10
    Res(1)=T
    Res(2)=P      
    phi=0
!    call CEA_HP(phi,Res)
!    call CEA_calculate(0,phi,Res) 
    write(*,*) 'P1=',Res(2),'T1=',Res(1),'H=',Res(4),'S=',Res(6)

!    call CEA_calculate(0,phi,Res) 
    write(*,*) 'P1=',Res(2),'T1=',Res(1),'H=',Res(4),'S=',Res(6)
    
!    call CEA_INIT(log_path,names,params,reslt)
!    call CEA_calculate(0,phi,Res) 
    write(*,*) 'P1=',Res(2),'T1=',Res(1),'H=',Res(4),'S=',Res(6)

!    call CEA_calculate(0,phi,Res) 
    write(*,*) 'P1=',Res(2),'T1=',Res(1),'H=',Res(4),'S=',Res(6)

    pause
    stop
 
    T=3000
    P=10
    M=6
 
    
    do i=1,10
        M=i
        call NTDC_Stat_Total(phi,T,P,M,Res,reslt)   
            write(*,*) 'P2=',Res(2),'T2=',Res(1),'H=',Res(4),'S=',Res(6)    
        call NTDC_Total_Stat(phi,Res(1),Res(2),M,Res,reslt)   
            write(*,*) 'P1=',Res(2),'T1=',Res(1),'H=',Res(4),'S=',Res(6)
        if (abs(Res(1)-T)>0.001) print *,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        if (abs(Res(2)-p)>0.001) print *,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"        
    enddo
    
    ! Ambient Parameters 
 !     call RATU (287.d0,31261.d0,T,PH,DP,ROH,CS)

    pause
    end program Test1

