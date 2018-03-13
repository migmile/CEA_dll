
      BLOCKDATA 
C***********************************************************************
C FUNDAMENTAL CONSTANTS FROM:  COHEN,E.RICHARD & TAYLOR,BARRY N.,
C THE 1986 CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL PHYSICAL
C CONSTANTS, J.PHYS.CHEM.REF.DATA, VOL.17, NO.4, 1988, PP 1795-1803.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C
      DATA Rr/8314.51D0/,Pi/3.14159265D0/,Avgdr/6.0221367D0/,
     &     Boltz/1.380658D0/
C ATOMIC SYMBOLS
      DATA Symbol/'H ','D ','HE','LI','BE','B ','C ','N ','O ','F ',
     &     'NE','NA','MG','AL','SI','P ','S ','CL','AR','K ','CA','SC',
     &     'TI','V ','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS',
     &     'SE','BR','KR','RB','SR','Y ','ZR','NB','MO','TC','RU','RH',
     &     'PD','AG','CD','IN','SN','SB','TE','I ','XE','CS','BA','LA',
     &     'CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM',
     &     'YB','LU','HF','TA','W ','RE','OS','IR','PT','AU','HG','TL',
     &     'PB','BI','PO','AT','RN','FR','RA','AC','TH','PA','U ','NP',
     &     'PU','AM','CM','BK','CF','ES'/
C
C  ATOMIC WEIGHTS - Coplen,T.B., Atomic Weights of the Elements 1999. 
C     J.Phys.Chem.Ref.Data, vol.30, no.3, 2001, pp.701-712.
C
      DATA atmwt/               1.00794D0,2.014102D0,4.002602D0,6.941D0,
     1 9.012182D0,10.811D0,12.0107D0,14.0067D0,15.9994D0,18.9984032D0,
     2 20.1797D0,22.989770D0,24.305D0,26.981538D0,28.0855D0,30.973761D0,
     3 32.065D0,35.453D0,39.948D0,39.0983D0,40.078D0,44.95591D0,
     4 47.867D0, 50.9415D0,51.9961D0,54.938049D0,
     5 55.845D0,58.933200D0,58.6934D0,63.546D0,65.39D0,69.723D0,72.64D0,
     6 74.92160D0,78.96D0,79.904D0,83.80D0,85.4678D0,87.62D0,88.90585D0,
     7 91.224D0,92.90638D0,95.94D0,97.9072D0,101.07D0,102.9055D0,
     $ 106.42D0,
     8 107.8682D0,112.411D0,114.818D0,118.710D0, 121.760D0,127.6D0,
     9 126.90447D0,131.293D0,132.90545D0,137.327D0,138.9055D0,140.116D0,
     $ 140.90765D0,144.9127D0,145.D0,150.36D0,151.964D0,157.25D0,
     $ 158.92534D0,
     $ 162.50D0,164.93032D0,167.259D0,168.93421D0,173.04D0,174.967D0,
     $ 178.49D0,180.9479D0,183.84D0,186.207D0,190.23D0,192.217D0,
     $ 195.078D0,196.96655D0,200.59D0,204.3833D0,207.2D0,208.98038D0,
     $ 208.9824D0, 209.9871D0,
     $ 222.0176D0,223.0197D0,226.0254D0,227.0278D0,232.0381D0,
     $ 231.03588D0,238.02891D0,237.0482D0,244.0642D0,243.0614D0,
     $ 247.0703D0,247.0703D0,251.0587D0,252.083D0/
C ATOMIC VALENCES
      DATA Valnce/1.,1.,0.,1.,2.,3.,4.,0., - 2., - 1.,0.,1.,2.,3.,4.,5.,
     &     4., - 1.,0.,1.,2.,3.,4.,5.,3.,2.,3.,2.,2.,2.,2.,3.,4.,3.,4.,
     &     - 1.,0.,1.,2.,3.,4.,5.,6.,7.,3.,3.,2.,1.,2.,3.,4.,3.,4.,
     &     - 1.,0.,1.,2.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,
     &     4.,5.,6.,7.,4.,4.,4.,3.,2.,1.,2.,3.,2., - 1.,0.,1.,2.,3.,4.,
     &     5.,6.,5.,4.,3.,3.,3.,3.,3./
C INFORMATION USED IN VARIABLE OUTPUT FORMAT
      DATA Fmt/'(1X',',A15',',','F9.','0,','F9.','0,','F9.','0,','F9.',
     &     '0,','F9.','0,','F9.','0,','F9.','0,','F9.','0,','F9.','0,',
     &     'F9.','0,','F9.','0,','F9.','0,','F9.','0',')'/
      END


      subroutine Calculate(phi,Results)
      !çàãëóøêà!
!         I/O ÏÀÐÀÌÅÒÐÛ: !
      IMPLICIT NONE
      INCLUDE 'cea.inc'
      
      REAL*8 DABS,Results(100)
      integer i
      real*8 xyz,denmtr,phi,eratio
      integer DATE_TIME(8)
      i=1
	
	if (Results(1)>=0) T(1)=Results(1)
	P(1)=Results(2)	

	! ïåðåñ÷åò íà÷àëüíîé ýíòàëüïèè èñõîäíûõ êîìïîíåíòîâ â çàâèñèìîñòè îò òåìïåðàòóðû
	CALL Enth_T()

      ! îáðóáàåòñÿ çàïóñê ïîñëå 2017 ãîäà
      CALL DATE_AND_TIME (values=DATE_TIME)
      DATE_TIME(1)=DATE_TIME(1)-18
      if (DATE_TIME(1)>2000) then
          Results=-1
          return
      endif 
      
      if (phi==0) then
            Nof = 1       ! êîë-âî òî÷åê ïî çðø
            Oxf(1) = 0.
            IF ( Wp(2).GT.0. ) THEN
              Oxf(1) = Wp(1)/Wp(2)
            ELSE
              i_hok=-4
            ENDIF
      ELSE        ! IF ( phi ) THEN
              Nof = 1       ! êîë-âî òî÷åê ïî çðø
              Oxf(1)=phi
              eratio = Oxf(1)
              xyz = -Vmin(2) - Vpls(2)
              denmtr = eratio*(Vmin(1)+Vpls(1))
              IF ( DABS(denmtr).LT.1.D-30 ) THEN
                  i_hok=-5 
              else 
                  Oxf(1) = xyz/denmtr
             ENDIF
          ENDIF
      
      Call Init2()
      
      IF ( Tp.OR.Hp.OR.Sp ) THEN
          CALL THERMP

         Results(1)= Ttt(i)                           ! ÒÅÌÏÅÐÀÒÓÐÀ 
         Results(2)= Ppp(i)                           ! Äàâëåíèå
         Results(3) = 1.D05/Vlm(i)                    ! ÏËÎÒÍÎÑÒÜ
         Results(4) = Hsum(i)*R*1000.	                ! ÝÍÒÀËÜÏÈß
         Results(5) = (Hsum(i)-Ppp(i)*Vlm(i)/Rr)*R    ! ÂÍÓÒÐÅÍÍßß ÝÍÅÐÃÈß
         Results(6) = Ssum(i)*R*1000.                 ! ÝÍÒÐÎÏÈß
         Results(7) = Wm(i)                           ! M, (1/n)        ',(Wm(j) MOLECULAR WEIGHT
         Results(8) = Cpr(i)*R*1000.                  ! ÒÅÏËÎÅÌÊÎÑÒÜ
         Results(9) = Gammas(i)                       ! k
         Results(10) = (Rr*Gammas(i)*Ttt(i)/Wm(i))**.5 !ÑÊÎÐÎÑÒÜ ÇÂÓÊÀ
      endif
      end

      SUBROUTINE CPHS
C***********************************************************************
C CALCULATES THERMODYNAMIC PROPERTIES FOR INDIVIDUAL SPECIES
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      REAL*8 cx(7),hcx(7),scx(7)
      INTEGER i,ij,j,jj,k
      SAVE i,ij,j,jj,k,scx
C
      DATA cx/2*0.,1.D0,.5D0,.6666666666666667D0,.75D0,.8D0/
      DATA hcx(3)/1.D0/
      k = 1
      IF ( Tt.GT.Tg(2) ) k = 2
      IF ( Tt.GT.Tg(3) ) k = 3
      cx(2) = 1.D0/Tt
      cx(1) = cx(2)**2
      scx(3) = Tln
      scx(2) = -cx(2)
      hcx(2) = Tln*cx(2)
      hcx(1) = -cx(1)
      scx(1) = hcx(1)*.5D0
      DO i = 4,7
        hcx(i) = cx(i)*Tt
        scx(i) = cx(i-1)*Tt
      ENDDO
      DO j = 1,Ng
        H0(j) = 0.D0
        S(j) = 0.D0
      ENDDO
      DO i = 7,4, - 1
        DO j = 1,Ng
          S(j) = (S(j)+Coef(j,i,k))*scx(i)
          H0(j) = (H0(j)+Coef(j,i,k))*hcx(i)
        ENDDO
      ENDDO
      DO i = 1,3
        DO j = 1,Ng
          S(j) = S(j) + Coef(j,i,k)*scx(i)
          H0(j) = H0(j) + Coef(j,i,k)*hcx(i)
        ENDDO
      ENDDO
      DO j = 1,Ng
        S(j) = S(j) + Coef(j,9,k)
        H0(j) = H0(j) + Coef(j,8,k)*cx(2)
      ENDDO
      IF ( .NOT.Tp.OR.Convg ) THEN
        DO j = 1,Ng
          Cp(j) = 0.D0
        ENDDO
        DO i = 7,4, - 1
          DO j = 1,Ng
            Cp(j) = (Cp(j)+Coef(j,i,k))*Tt
          ENDDO
        ENDDO
        DO i = 1,3
          DO j = 1,Ng
            Cp(j) = Cp(j) + Coef(j,i,k)*cx(i)
          ENDDO
        ENDDO
      ENDIF
      IF ( Npr.NE.0.AND.k.NE.3.AND.Ng.NE.Ngc ) THEN
        DO ij = 1,Npr
          j = Jcond(ij)
          jj = Jcond(ij) - Ng
          Cp(j) = 0.D0
          H0(j) = 0.D0
          S(j) = 0.D0
          DO i = 7,4, - 1
            S(j) = (S(j)+Cft(jj,i))*scx(i)
            H0(j) = (H0(j)+Cft(jj,i))*hcx(i)
            Cp(j) = (Cp(j)+Cft(jj,i))*Tt
          ENDDO
          DO i = 1,3
            S(j) = S(j) + Cft(jj,i)*scx(i)
            H0(j) = H0(j) + Cft(jj,i)*hcx(i)
            Cp(j) = Cp(j) + Cft(jj,i)*cx(i)
          ENDDO
          S(j) = S(j) + Cft(jj,9)
          H0(j) = H0(j) + Cft(jj,8)*cx(2)
        ENDDO
      ENDIF
      GOTO 99999
      ENTRY ALLCON
      DO jj = 1,Nc
        j = jj + Ng
        Cp(j) = 0.D0
        H0(j) = 0.D0
        S(j) = 0.D0
        DO i = 7,4, - 1
          S(j) = (S(j)+Cft(jj,i))*scx(i)
          H0(j) = (H0(j)+Cft(jj,i))*hcx(i)
          Cp(j) = (Cp(j)+Cft(jj,i))*Tt
        ENDDO
        DO i = 1,3
          S(j) = S(j) + Cft(jj,i)*scx(i)
          H0(j) = H0(j) + Cft(jj,i)*hcx(i)
          Cp(j) = Cp(j) + Cft(jj,i)*cx(i)
        ENDDO
        S(j) = S(j) + Cft(jj,9)
        H0(j) = H0(j) + Cft(jj,8)*cx(2)
      ENDDO
99999 END




      SUBROUTINE GAUSS
C***********************************************************************
C SOLVE ANY LINEAR SET OF UP TO MAXMAT EQUATIONS
C NUMBER OF EQUATIONS = IMAT
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,imatp1,j,k,nn,nnp1
      REAL*8 bigno,coefx(50),tmp
      REAL*8 DABS,DMAX1
      SAVE coefx,i,imatp1,j,k,nn,nnp1,tmp
C
      DATA bigno/1.E+25/
C BEGIN ELIMINATION OF NNTH VARIABLE
      imatp1 = Imat + 1
      DO nn = 1,Imat
        IF ( nn.NE.Imat ) THEN
C SEARCH FOR MAXIMUM COEFFICIENT IN EACH ROW
          nnp1 = nn + 1
          DO i = nn,Imat
            coefx(i) = bigno
            IF ( G(i,nn).NE.0. ) THEN
              coefx(i) = 0.
              DO j = nnp1,imatp1
                coefx(i) = DMAX1(coefx(i),DABS(G(i,j)))
              ENDDO
              tmp = DABS(G(i,nn))
              IF ( bigno*tmp.GT.coefx(i) ) THEN
                coefx(i) = coefx(i)/tmp
              ELSE
                coefx(i) = bigno
              ENDIF
            ENDIF
          ENDDO
C LOCATE ROW WITH SMALLEST MAXIMUM COEFFICIENT
          tmp = bigno
          i = 0
          DO j = nn,Imat
            IF ( coefx(j).LT.tmp ) THEN
              tmp = coefx(j)
              i = j
            ENDIF
          ENDDO
          IF ( i.EQ.0 ) THEN
            Msing = nn
            GOTO 99999
C INDEX I LOCATES EQUATION TO BE USED FOR ELIMINATING THE NTH
C VARIABLE FROM THE REMAINING EQUATIONS
C INTERCHANGE EQUATIONS I AND NN
          ELSEIF ( nn.NE.i ) THEN
            DO j = nn,imatp1
              tmp = G(i,j)
              G(i,j) = G(nn,j)
              G(nn,j) = tmp
            ENDDO
          ENDIF
        ELSEIF ( G(nn,nn).EQ.0 ) THEN
          Msing = nn
          GOTO 99999
        ENDIF
C DIVIDE NTH ROW BY NTH DIAGONAL ELEMENT AND ELIMINATE THE NTH
C VARIABLE FROM THE REMAINING EQUATIONS
        k = nn + 1
        tmp = G(nn,nn)
        IF ( tmp.EQ.0. ) THEN
          Msing = nn
          GOTO 99999
        ELSE
          DO j = k,imatp1
            G(nn,j) = G(nn,j)/tmp
          ENDDO
          IF ( k.NE.imatp1 ) THEN
            DO i = k,Imat
CDIR$ IVDEP
              DO j = k,imatp1
                G(i,j) = G(i,j) - G(i,nn)*G(nn,j)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
C BACKSOLVE FOR THE VARIABLES
      k = Imat
 100  j = k + 1
      X(k) = 0.0D0
      tmp = 0.0
      IF ( Imat.GE.j ) THEN
        DO i = j,Imat
          tmp = tmp + G(k,i)*X(i)
        ENDDO
      ENDIF
      X(k) = G(k,imatp1) - tmp
      k = k - 1
      IF ( k.NE.0 ) GOTO 100
99999 END



      SUBROUTINE MATRIX
C***********************************************************************
C SET UP ITERATION OR DERIVATIVE MATRIX.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,iq,iq2,iq3,isym,j,k,kk,kmat
      REAL*8 energyl,f,h,ss,sss,term,term1
      SAVE energyl,f,h,i,iq,iq2,iq3,isym,j,k,kk,kmat,ss,sss,term,term1
C
      iq = Nlm + Npr
      Iq1 = iq + 1
      iq2 = Iq1 + 1
      iq3 = iq2 + 1
      kmat = iq3
      IF ( .NOT.Convg.AND.Tp ) kmat = iq2
      Imat = kmat - 1
C CLEAR MATRIX STORAGES TO ZERO
      DO i = 1,Imat
        DO k = 1,kmat
          G(i,k) = 0.0D0
        ENDDO
      ENDDO
      G(iq2,Iq1) = 0.D0
      sss = 0.D0
      Hsum(Npt) = 0.D0
C BEGIN SET-UP OF ITERATION OR DERIVATIVE MATRIX
      DO j = 1,Ng
        Mu(j) = H0(j) - S(j) + Enln(j) + Tm
        IF ( En(j,Npt).NE.0.D0 ) THEN
          h = H0(j)*En(j,Npt)
          f = Mu(j)*En(j,Npt)
          ss = h - f
          term1 = h
          IF ( kmat.EQ.iq2 ) term1 = f
          DO i = 1,Nlm
            IF ( A(i,j).NE.0. ) THEN
              term = A(i,j)*En(j,Npt)
              DO k = i,Nlm
                G(i,k) = G(i,k) + A(k,j)*term
              ENDDO
              G(i,Iq1) = G(i,Iq1) + term
              G(i,iq2) = G(i,iq2) + A(i,j)*term1
              IF ( .NOT.(Convg.OR.Tp) ) THEN
                G(i,iq3) = G(i,iq3) + A(i,j)*f
                IF ( Sp ) G(iq2,i) = G(iq2,i) + A(i,j)*ss
              ENDIF
            ENDIF
          ENDDO
          IF ( kmat.NE.iq2 ) THEN
            IF ( Convg.OR.Hp ) THEN
              G(iq2,iq2) = G(iq2,iq2) + H0(j)*h
              IF ( .NOT.Convg ) THEN
                G(iq2,iq3) = G(iq2,iq3) + H0(j)*f
                G(Iq1,iq3) = G(Iq1,iq3) + f
              ENDIF
            ELSE
              G(iq2,Iq1) = G(iq2,Iq1) + ss
              G(iq2,iq2) = G(iq2,iq2) + H0(j)*ss
              G(iq2,iq3) = G(iq2,iq3) + Mu(j)*ss
              G(Iq1,iq3) = G(Iq1,iq3) + f
            ENDIF
          ENDIF
          G(Iq1,iq2) = G(Iq1,iq2) + term1
        ENDIF
      ENDDO
C CONDENSED SPECIES
      IF ( Npr.NE.0 ) THEN
        DO k = 1,Npr
          j = Jcond(k)
          kk = Nlm + k
          Mu(j) = H0(j) - S(j)
          DO i = 1,Nlm
            G(i,kk) = A(i,j)
            G(i,kmat) = G(i,kmat) - A(i,j)*En(j,Npt)
          ENDDO
          G(kk,iq2) = H0(j)
          G(kk,kmat) = Mu(j)
          Hsum(Npt) = Hsum(Npt) + H0(j)*En(j,Npt)
          IF ( Sp ) THEN
            sss = sss + S(j)*En(j,Npt)
            G(iq2,kk) = S(j)
          ENDIF
        ENDDO
      ENDIF
      sss = sss + G(iq2,Iq1)
      Hsum(Npt) = Hsum(Npt) + G(Iq1,iq2)
      G(Iq1,Iq1) = Sumn - Enn
C REFLECT SYMMETRIC PORTIONS OF THE MATRIX
      isym = Iq1
      IF ( Hp.OR.Convg ) isym = iq2
      DO i = 1,isym
CDIR$ IVDEP
        DO j = i,isym
          G(j,i) = G(i,j)
        ENDDO
      ENDDO
C COMPLETE THE RIGHT HAND SIDE
      IF ( .NOT.Convg ) THEN
        DO i = 1,Nlm
          G(i,kmat) = G(i,kmat) + B0(i) - G(i,Iq1)
        ENDDO
        G(Iq1,kmat) = G(Iq1,kmat) + Enn - Sumn
C COMPLETE ENERGY ROW AND TEMPERATURE COLUMN
        IF ( kmat.NE.iq2 ) THEN
          IF ( Sp ) energyl = S0 + Enn - Sumn - sss
          IF ( Hp ) energyl = Hsub0/Tt - Hsum(Npt)
          G(iq2,iq3) = G(iq2,iq3) + energyl
          G(iq2,iq2) = G(iq2,iq2) + Cpsum
        ENDIF
      ELSE
        IF ( Pderiv ) THEN
C PDERIV = .TRUE.-- SET UP MATRIX TO SOLVE FOR DLVPT
          G(Iq1,iq2) = Enn
          DO i = 1,iq
            G(i,iq2) = G(i,Iq1)
          ENDDO
        ENDIF
        G(iq2,iq2) = G(iq2,iq2) + Cpsum
      ENDIF
      IF ( Vol.AND..NOT.Convg ) THEN
C CONSTANT VOLUME MATRIX
        IF ( kmat.EQ.iq2 ) THEN
          DO i = 1,iq
            G(i,Iq1) = G(i,iq2)
          ENDDO
        ELSE
CDIR$ IVDEP
          DO i = 1,iq
            G(Iq1,i) = G(iq2,i) - G(Iq1,i)
            G(i,Iq1) = G(i,iq2) - G(i,Iq1)
            G(i,iq2) = G(i,iq3)
          ENDDO
          G(Iq1,Iq1) = G(iq2,iq2) - G(Iq1,iq2) - G(iq2,Iq1)
          G(Iq1,iq2) = G(iq2,iq3) - G(Iq1,iq3)
          IF ( Hp ) G(Iq1,iq2) = G(Iq1,iq2) + Enn
        ENDIF
        kmat = Imat
        Imat = Imat - 1
      ENDIF
      END


      SUBROUTINE NEWOF
C***********************************************************************
C CALCULATE NEW VALUES OF B0 AND HSUB0 FOR NEW OF RATIO
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER i,j
      REAL*8 assval,bigb,bratio,dbi,smalb,tem,v1,v2
      REAL*8 DABS,DLOG
      SAVE assval,bigb,bratio,dbi,i,j,smalb,tem,v1,v2
C
      IF ( .NOT.Short ) WRITE (IOOUT,99001) Oxfl
      Eqrat = 0.
      tem = Oxfl + 1.
      v2 = (Oxfl*Vmin(1)+Vmin(2))/tem
      v1 = (Oxfl*Vpls(1)+Vpls(2))/tem
      IF ( v2.NE.0. ) Eqrat = DABS(v1/v2)
      DO i = 1,Nlm
        B0(i) = (Oxfl*B0p(i,1)+B0p(i,2))/tem
        dbi = DABS(B0(i))
        IF ( i.EQ.1 ) THEN
          bigb = dbi
          smalb = dbi
        ELSEIF ( dbi.NE.0. ) THEN
          IF ( dbi.LT.smalb ) smalb = dbi
          IF ( dbi.GT.bigb ) bigb = dbi
        ENDIF
      ENDDO
      Bcheck = bigb*.000001D0
C CALCUALTE MOLECULAR WEIGHT OF TOTAL REACTANT, WMIX.
      IF ( Am(1).NE.0.0.AND.Am(2).NE.0.0 ) THEN
        Wmix = (Oxfl+1.)*Am(1)*Am(2)/(Am(1)+Oxfl*Am(2))
      ELSE
        Wmix = Am(2)
        IF ( Am(2).EQ.0.0 ) Wmix = Am(1)
      ENDIF
      Npt = 1
C IF ASSIGNED U OR H NOT GIVEN IN PROB DATA, INITIAL HSUB0 = 1.D30
      IF ( Size.EQ.0. ) assval = Hsub0
      IF ( assval.GE.1.D30 ) Hsub0 = (Oxfl*Hpp(1)+Hpp(2))/tem
C NOTE THAT "BRATIO" IS "BRATIO" IN SEC 3.2 IN RP-1311.
      bratio = smalb/bigb
      Size = 18.420681D0
      IF ( bratio.LT.1.D-5 ) Size = DLOG(1000.D0/bratio)
      Jsol = 0
      Jliq = 0
      IF ( .NOT.Short ) THEN
        WRITE (IOOUT,99002)
        IF ( Vol ) WRITE (IOOUT,99003)
        IF ( .NOT.Vol ) WRITE (IOOUT,99004)
        WRITE (IOOUT,99005) Hpp(2),Hpp(1),Hsub0
        WRITE (IOOUT,99006)
      ENDIF
      DO i = 1,Nlm
        j = Jcm(i)
        IF ( .NOT.Short ) WRITE (IOOUT,99007) Prod(j),B0p(i,2),B0p(i,1),
     &                           B0(i)
      ENDDO
      RETURN
99001 FORMAT (/' O/F = ',F10.6)
99002 FORMAT (/,23X,'EFFECTIVE FUEL',5X,'EFFECTIVE OXIDANT',8X,
     &        'MIXTURE')
99003 FORMAT (' INTERNAL ENERGY',11X,'u(2)/R',14X,'u(1)/R',14X,'u0/R')
99004 FORMAT (' ENTHALPY',18X,'h(2)/R',14X,'h(1)/R',15X,'h0/R')
99005 FORMAT (' (KG-MOL)(K)/KG',4X,E18.8,2E20.8)
99006 FORMAT (/' KG-FORM.WT./KG',13X,'bi(2)',15X,'bi(1)',15X,'b0i')
99007 FORMAT (1X,A16,3E20.8)
      END




      SUBROUTINE SETEN
C***********************************************************************
C USE COMPOSITIONS FROM PREVIOUS POINT AS INITIAL ESTIMATES FOR
C CURRENT POINT NPT.  IF -
C  ISV>0  USE COMPOSITIONS FROM POINT ISV.
C  ISV<0  SAVE COMPOSITIONS FROM POINT -ISV FOR POSSIBLE LATER USE.
C         ALSO USE COMPOSITIONS FROM POINT -ISV FOR NPT.
C  ISV=0  USE COMPOSITIONS SAVED WHEN ISV<0.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER j,lsav
      REAL*8 tsave
      REAL*8 DEXP
      SAVE j,lsav,tsave
C
      IF ( Isv.LT.0 ) THEN
C FIRST T--SAVE COMPOSITIONS FOR FUTURE POINTS WITH THIS T
        Isv = -Isv
        tsave = Ttt(Isv)
        Ensave = Enn
        Enlsav = Ennl
        lsav = Lsave
        DO j = 1,Ng
          Sln(j) = Enln(j)
        ENDDO
        DO j = 1,Ng
          En(j,Npt) = En(j,Isv)
        ENDDO
        Npr = 0
        DO j = Ngp1,Ngc
          Sln(j) = En(j,Isv)
          En(j,Npt) = Sln(j)
          IF ( Jliq.EQ.j ) THEN
            En(Jsol,Npt) = En(Jsol,Isv) + En(Jliq,Isv)
            En(Jliq,Npt) = 0.
            Jsol = 0
            Jliq = 0
            tsave = tsave - 5.
            Tt = tsave
            Sln(j) = 0.
          ELSEIF ( En(j,Npt).GT.0. ) THEN
            Npr = Npr + 1
            Jcond(Npr) = j
          ENDIF
        ENDDO
      ELSEIF ( Isv.EQ.0 ) THEN
C NEXT POINT FIRST T IN SCHEDULE, USE PREVIOUS COMPOSITIONS FOR THIS T
        Jsol = 0
        Jliq = 0
        Enn = Ensave
        Ennl = Enlsav
        Lsave = lsav
        Npr = 0
        DO j = Ngp1,Ngc
          En(j,Npt) = Sln(j)
          IF ( En(j,Npt).GT.0.D0 ) THEN
            Npr = Npr + 1
            Jcond(Npr) = j
          ENDIF
        ENDDO
        DO j = 1,Ng
          En(j,Npt) = 0.
          Enln(j) = Sln(j)
          IF ( Sln(j).NE.0. ) THEN
            IF ( (Enln(j)-Ennl+18.5).GT.0. ) En(j,Npt) = DEXP(Enln(j))
          ENDIF
        ENDDO
        IF ( .NOT.Tp ) Tt = tsave
        Sumn = Enn
      ELSEIF ( Isv.GT.0 ) THEN
C USE COMPOSITIONS FROM PREVIOUS POINT
        DO j = 1,Ngc
          En(j,Npt) = En(j,Isv)
        ENDDO
      ENDIF
      END

      SUBROUTINE THERMP
C***********************************************************************
C ASSIGNED THERMODYNAMIC STATES.  HP,SP,TP,UV,SV, AND TV PROBLEMS.
C***********************************************************************
      IMPLICIT NONE
      INCLUDE 'cea.inc'
C LOCAL VARIABLES
      INTEGER iof
      LOGICAL uv,tv,sv
      SAVE iof
C
      EQUIVALENCE (Hp,Uv)
      EQUIVALENCE (Tp,Tv)
      EQUIVALENCE (Sp,Sv)
      Eql = .TRUE.
      DO 100 iof = 1,Nof
        Oxfl = Oxf(iof)
        CALL NEWOF
C SET ASSIGNED P OR VOLUME
        DO Ip = 1,Np
          Pp = P(Ip)
C SET ASSIGNED T
          DO It = 1,Nt
            Vv = V(Ip)
            Tt = T(It)
            CALL EQLBRM
            IF ( Npt.EQ.0 ) GOTO 200
            IF ( Trnspt.AND.Tt.NE.0. ) CALL TRANP
            Isv = 0
            IF ( Ip.NE.Np.OR.It.NE.Nt.AND.Tt.NE.0. ) THEN
              Isv = Npt
              IF ( Npt.NE.NCOL ) GOTO 10
            ENDIF
!            IF ( .NOT.Hp ) WRITE (IOOUT,99001)
!             IF ( Hp ) WRITE (IOOUT,99002)
!            IF ( .NOT.Vol ) THEN
!              IF ( Hp ) WRITE (IOOUT,99006)
!              IF ( Tp ) WRITE (IOOUT,99007)
!              IF ( Sp ) WRITE (IOOUT,99008)
!            ELSE
!              IF ( Uv ) WRITE (IOOUT,99003)
!              IF ( Tv ) WRITE (IOOUT,99004)
!              IF ( Sv ) WRITE (IOOUT,99005)
!            ENDIF
!            CALL OUT1
!            WRITE (IOOUT,99009)
!            CALL OUT2
!            IF ( Trnspt ) CALL OUT4
!            CALL OUT3
            Iplt = MIN(Iplt+Npt,500)
            IF ( Isv.EQ.0.AND.iof.EQ.Nof ) GOTO 200
!            WRITE (IOOUT,99010)
            Npt = 0
 10         Npt = Npt + 1
            IF ( .NOT.Tp.AND.Tt.NE.0. ) T(1) = Tt
            IF ( Nt.EQ.1.AND.Np.EQ.1 ) GOTO 100
            IF ( Ip.EQ.1.AND.It.EQ.1 ) Isv = -Isv
            IF ( Nt.NE.1 ) THEN
              IF ( It.EQ.Nt.OR.Tt.EQ.0. ) Isv = 0
            ENDIF
            CALL SETEN
          ENDDO
        ENDDO
 100  CONTINUE
 200  RETURN
99001 FORMAT (////15X,'THERMODYNAMIC EQUILIBRIUM PROPERTIES AT ASSIGNED'
     &        )
99002 FORMAT (////9X,
     &     'THERMODYNAMIC EQUILIBRIUM COMBUSTION PROPERTIES AT ASSIGNED'
     &     )
99003 FORMAT (/36X,' VOLUME'/)
99004 FORMAT (/28X,'TEMPERATURE AND VOLUME'/)
99005 FORMAT (/30X,'ENTROPY AND VOLUME'/)
99006 FORMAT (/34X,' PRESSURES'/)
99007 FORMAT (/27X,'TEMPERATURE AND PRESSURE'/)
99008 FORMAT (/29X,'ENTROPY AND PRESSURE'/)
99009 FORMAT (/' THERMODYNAMIC PROPERTIES'/)
99010 FORMAT (////)
      END
      
     