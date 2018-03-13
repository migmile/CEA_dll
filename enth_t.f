	!integer l(10)	!!!!!- в глобальные

	Subroutine Enth_T()
	IMPLICIT NONE
	INCLUDE 'cea.inc'
	
	integer n,kr,l


	Hpp=0
		
	DO n = 1,Nreac

        Tt = Rtemp(n)

                Tln = DLOG(Rtemp(n))
		l=l_t(n)
      if (l>0) then          
        Enth(n) = (((((rcf(n,7,l)/5.D0)*Tt+rcf(n,6,l)/4.D0)*Tt+
     &  rcf(n,5,l)/3.D0)*Tt+rcf(n,4,l)/2.D0)*Tt+rcf(n,3,l))
     &  *Tt - rcf(n,1,l)/Tt + rcf(n,2,l)*Tln + rcf(n,8,l)
      end if

        IF ( Fox(n)(:1).EQ.'F' ) THEN
          kr = 2
        ELSEIF ( Fox(n)(:4).EQ.'NAME' ) THEN
          kr = 2
        ELSE
          kr = 1
        ENDIF


        Hpp(kr) = Hpp(kr) + Enth(n)*Pecwt(n)/Rmw(n) 	     
	end do
	end