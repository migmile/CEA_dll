      subroutine CalculateFractions(ProductNames,X_,n_pr)
      IMPLICIT NONE
      INCLUDE 'cea.inc'      
      character*15 ProductNames(100)      ! products names
      real X_(100)                        ! fractions of product
      integer n_pr                        ! quantaty of products with sugnificant fractions
      
      integer i,k,m,im,kin,ione,notuse
      real*8 tra,tem
      LOGICAL kok
      
      n_pr=0
      tra = 5.D-6
      IF ( Trace.NE.0. ) tra = Trace
C MASS OR MOLE FRACTIONS 
      IF ( Eql ) THEN
       notuse = 0
        DO k = 1,Ngc
          kok = .TRUE.
          IF ( k.GT.Ng.AND.k.LT.Ngc.AND.Prod(k).EQ.Prod(k+1) ) THEN
            kok = .FALSE.
            im = 0
            GOTO 120
          ENDIF
          DO m = 1,Nplt
            im = 0
            IF ( Pltvar(m).EQ.Prod(k).OR.'*'//Pltvar(m).EQ.Prod(k) )
     &           THEN
              im = m
              GOTO 120
            ENDIF
          ENDDO
 120      kin = 0
          DO i = 1,Npt
            IF ( Massf ) THEN
              tem = Mw(k)
            ELSE
              tem = 1.D0/Totn(i)
            ENDIF
            IF ( k.LE.Ng ) THEN
              X(i) = En(k,i)*tem
            ELSE
              IF ( Prod(k).NE.Prod(k-1) ) X(i) = 0.D0
              IF ( En(k,i).GT.0.D0 ) X(i) = En(k,i)*tem
            ENDIF
            IF ( Nplt.NE.0.AND.i.GT.ione.AND.im.GT.0 )
     &           Pltout(i+Iplt-ione,im) = X(i)
            IF ( kok.AND.X(i).GE.tra ) kin = 1
          ENDDO
          IF ( kin.EQ.1 ) THEN

              n_pr=n_pr+1  
              ProductNames(n_pr) = Prod(k)
              X_(n_pr)= X(1)

            IF ( Prod(k).EQ.Omit(notuse) ) notuse = notuse - 1
          ELSEIF ( Prod(k).NE.Prod(k-1) ) THEN
            notuse = notuse + 1
            Omit(notuse) = Prod(k)
          ENDIF
        ENDDO
      ENDIF
      end subroutine

!____________________________________________________________________________
       subroutine  MixParams(arr )
!      
!     Подпрограмма возвращает соотношения компонентов смеси после последнего расчета равновесия:
!      O/F;   %FUEL;  R,EQ.RATIO;  PHI,EQ.RATIO;
!      
      IMPLICIT NONE
      real*8 :: arr(100)
      INCLUDE 'cea.inc'
      
      real*8 ::phi,tem,pfuel,rho
      integer n
      !      DO n = 1,Nreac
!        WRITE (IOOUT,99007) Fox(n),Rname(n),Pecwt(n),Enth(n)*R,Rtemp(n)
!      ENDDO
      phi = 0.
      tem = (Vpls(1)+Vmin(1))*Oxfl
      IF ( ABS(tem).GE.1.D-3 ) phi = -(Vmin(2)+Vpls(2))/tem
      IF ( Fox(1).EQ.'NAME' ) THEN
        pfuel = 0.
      ELSE
        pfuel = 100.D0/(1.D0+Oxfl)
      ENDIF
      IF ( Rh(1).NE.0..OR.Rh(2).NE.0. ) THEN
        IF ( Rh(1).EQ.0..OR.Rh(2).EQ.0. ) THEN
          rho = MAX(Rh(1),Rh(2))
        ELSE
          rho = (Oxfl+1.)*Rh(1)*Rh(2)/(Rh(1)+Oxfl*Rh(2))
        ENDIF
        rho = rho*1000.D0
      ENDIF
      arr(1)=Oxfl
      arr(2)=pfuel
      arr(3)=Eqrat
      arr(4)=phi
      
      DO n = 1,Nreac
          arr(4+n)=Enth(n)*R
      end do
      
      end
      