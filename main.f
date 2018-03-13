! ___________________________________________________________________________
      subroutine CEA_calculate(c_t,phi,Res)
      IMPLICIT NONE
      INCLUDE 'cea.inc'
      
      real*8 Res(100),MPars(100),phi      
	integer c_t,i

      Sp = .false.

	if (c_t==0) 	then
		tp = .FALSE.
		Hp = .TRUE.
	elseif	(c_t==1) then
		TP = .TRUE.
		Hp = .FALSE.
	end if

	call Calculate(phi,Res)
      
      call  MixParams(MPars)
      
      DO i = 1,Nreac+4
          Res(10+i)=MPars(i)
      end do
      
	tp = .FALSE.
	Hp = .FALSE.      
      
      END 
