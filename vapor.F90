#include <misc.h>
#include <params.h>

module vapor
!------------------------------------------------------------------------------------------------
! Purpose: Set the consts or variables related to water vapor 
! Original version : QPDATA.f & MPDATA.f & AVNEGQ.f (IAP 9L)
! Reconstructed  : ZhangHe
! Completed : 2005.9.12
! Update : 2007.8.20, ZhangHe, 
!          2007.12.21, ZhangHe
!          2007.12.23, ZhangHe, add dumb parameter 'IORD'
!          2008.5, WuJianping, parallel version
!          2008.6.11, ZhangHe, available for both serial & parallel
!------------------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use IAP_grid,  only : NX, NY, NM, NL, NZ, IB, IE, IM, JB, JE, period
   use mathconst, only : ZERO, HALF, ONE, TWO
   use pmgrid,    only : beglatdyn, endlatdyn, loc_JB, loc_JE,    &
                         beglatdynex, endlatdynex
!!   use Dyn_const, only: IORD, ISOR

   implicit none
   save
   public

   integer  :: IQVADV   ! switch of vertical advection computation of Q(water vapor)
   integer  :: IDQADV   ! switch on choosing the scheme of Q advection
   integer  :: IQHADV   ! switch of horizontal vertical advection computation of Q
   integer  :: IVMIDT   ! swich of using wind of middle time layer in computing Q advection 
   integer  :: NONOS    ! NON-OSSCILATORY OPTION
!!   integer  :: IORD     ! times of iterations in water vapor scheme
   integer  :: ISOR     ! switch of third order modification of water vapor scheme
   real(r8) :: EP       ! 
   real(r8) :: UVW0(NX,NL,NY,3)   ! wind velocity used for Q advection
   real(r8) :: WW(NX,NY)
!!   real(r8) :: DONOR,A,Y1,Y2,VDYF,R,X1,X2,VCOR,B,YY,VC31,VC32,PP,Y,PN  !2007.8.20

! Set data  (Initialize)  
   data IQVADV / -1 /
   data IQHADV / 2  /
   data IDQADV / 1  /
   data IVMIDT / 1  /
!   data NONOS  / 0  /
   data NONOS  / 1  /
!!   data IORD   / 3  /
!!   data IORD   / 1  /
   data ISOR   / 3  /
!!   data ISOR   / 1  /
   data EP     / 1.E-10 /

!================================================================================================
CONTAINS
!================================================================================================

!================================================================================================
   SUBROUTINE QPDATA( IFILTQ,IBCFFT )
!------------------------------------------------------------------------------------------------
!  Purpose: PREDICT POSITIVE DEFINITE FIELD Q  DUE TO 3-D ADVECTION
!           BY USING THE SCHEME PROPOSED BY P.K.Smolarkiewicz
!------------------------------------------------------------------------------------------------
      use IAP_prog   !! , only: U, V, WS, QT
	  use tendency, only: DQ, DQliq, DQice
	  implicit none
!----------------------------------------Arguments-----------------------------------------------
      integer, intent(in) :: IFILTQ   ! identification if SMOOTHING & FILTER are REQUESTED
	  integer, intent(in) :: IBCFFT   ! switch of filter
!-------------------------------------Local workspace-------------------------------------------
      real(r8) :: UQ(NX,NL,NY),VQ(NX,NL,NY),WQ(NX,NL,NY),PQ(NX,NY)
      integer  :: I, J, K
!=============================== zhh ====================================
      real(r8) :: DQmax, DQm
      real(r8) :: DQh(NX,NL,NY),DQv(NX,NL,NY)
      integer  :: Imax, Jmax, Kmax
!============================ 2007.12.5 =================================
!------------------------------------------------------------------------------------------------
!=============================== zhh ====================================
      DQmax = ZERO
!!      write (6,*) 'At beginning of sub. QPDATA: IORD = ', IORD, 'NONOS =', NONOS, 'ISOR =', ISOR
!!      write (6,*) 'At beginning of sub. QPDATA: '
!!      write (6,*) 'QT(1,128,1) =', QT(1,128,1), 'QT(160,127,1) =', QT(160,127,1)
!!      write (6,*) 'QTliq(1,128,1) =', QTliq(1,128,1), 'QTliq(160,127,1) =', QTliq(160,127,1)
!!      write (6,*) 'QTice(1,128,1) =', QTice(1,128,1), 'QTice(160,127,1) =', QTice(160,127,1)
!============================ 2007.12.8 =================================
!     GET THE ADVECTION VELOCITY
!$DOACROSS LOCAL(K,J,I)
      IF ( IVMIDT.EQ.+1 ) THEN
         DO J = beglatdyn, endlatdyn
            DO K = 1 ,NL
               DO I = 1 ,NX
                  UQ(I,K,J)      = HALF*(U (I,K,J)+UVW0(I,K,J,1))
                  VQ(I,K,J)      = HALF*(V (I,K,J)+UVW0(I,K,J,2))
                  WQ(I,K,J)      = HALF*(WS(I,K,J)+UVW0(I,K,J,3))
                  UVW0(I,K,J,1)  = U (I,K,J)
                  UVW0(I,K,J,2)  = V (I,K,J)
                  UVW0(I,K,J,3)  = WS(I,K,J)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO J = beglatdyn, endlatdyn
            DO K = 1 ,NL
               DO I = 1 ,NX
                  UQ(I,K,J)  = U (I,K,J)
                  VQ(I,K,J)  = V (I,K,J)
                  WQ(I,K,J)  = WS(I,K,J)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!     SAVE THE FIELD ON INPUT
      DO J = beglatdyn, endlatdyn
         DO K = 1 ,NL
            DO I = 1 ,NX
               DQ(I,K,J) = QT(I,K,J)
!======================== added by zhh ============================
               DQliq(I,K,J) = QTliq(I,K,J)
               DQice(I,K,J) = QTice(I,K,J)
!========================= 2007.12.21 =============================
            ENDDO
         ENDDO
      ENDDO
!
!     PERFORM HORIZONTAL ADVECTION IN SPHERICAL GEOMETRY
!$DOACROSS LOCAL(J,I)
      DO J = beglatdyn, endlatdyn
         IF ( IDQADV.EQ.+1 ) THEN
            DO I = 1 ,NX
               PQ(I,J)    = ONE
            ENDDO
         ELSE IF( IDQADV.EQ.+2 ) THEN
            DO I = 1 ,NX
               PQ(I,J)    = PT(I,J)
            ENDDO
         ELSE
            DO I = 1 ,NX
               PQ(I,J)    = P(I,J)
            ENDDO
         ENDIF
      ENDDO
!     
!------------------------ DO THE 2-D ADVECTION BY MPDATA ------------------------
      CALL MPDATA( QT,UQ,VQ,PQ,3 )
!======================== added by zhh ============================
      CALL MPDATA( QTliq,UQ,VQ,PQ,1 )
      CALL MPDATA( QTice,UQ,VQ,PQ,1 )
!!      CALL MPDATA( QTliq,UQ,VQ,PQ,3 )
!!      CALL MPDATA( QTice,UQ,VQ,PQ,3 )
!========================= 2007.12.21 =============================
!--------------------------------------------------------------------------------
!     PERFORM THE VERTICAL ADVECTION
      IF ( IQVADV.EQ.-1 ) THEN
!        BY  P.K.Smolarkiewicz
!--------------------------------------------------------------------------------
         CALL VPDATA( QT,WQ,3 )  
!======================== added by zhh ============================
         CALL VPDATA( QTliq,WQ,1 )
         CALL VPDATA( QTice,WQ,1 )
!!         CALL VPDATA( QTliq,WQ,3 )
!!         CALL VPDATA( QTice,WQ,3 )
!========================= 2007.12.21 =============================
!--------------------------------------------------------------------------------
      ENDIF
!
!     CALCULATE THE TENDENCE
! filter DQ
      DO J = beglatdyn, endlatdyn
         DO K = 1 ,NL
            DO I = 1 ,NX
               UQ(I,K,J) = QT(I,K,J)
               VQ(I,K,J) = DQ(I,K,J)
            ENDDO
         ENDDO
      ENDDO
      DO J = beglatdyn, endlatdyn
         DO K = 1 ,NL
            DO I = 1 ,NX
               DQ(I,K,J) = UQ(I,K,J) - VQ(I,K,J)
            ENDDO
         ENDDO
      ENDDO
!
!=============================== zhh ====================================
!!    write (6,*) 'befrore filter: '
!!    write (6,*) 'DQ(1,128,1) =', DQ(1,128,1), 'DQ(160,127,1) =', DQ(160,127,1)
!============================ 2007.12.8 =================================
!     PERFORM SMOOTHING & FILTER IF REQUESTED
      IF ( IFILTQ.NE.0 ) THEN
!$DOACROSS LOCAL(K,J,I)
         DO K = 1 ,NL
            WW(:,beglatdyn:endlatdyn)=DQ(:,K,beglatdyn:endlatdyn)
            CALL FILT2D( WW,0,1,IBCFFT )
            DQ(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
            DO J = beglatdyn, endlatdyn
               DO I = 1 ,NX
                  QT(I,K,J) = VQ(I,K,J) + DQ(I,K,J)
!=============================== zhh =====================================
                  if ( abs(DQ(I,K,J)) > DQmax) then
                     DQmax = abs(DQ(I,K,J))
                     DQm   = DQ(I,K,J)
                     Imax = I
                     Jmax = J
                     Kmax = K
                  end if
!============================ 2007.12.5 ==================================
               ENDDO
            ENDDO
            DO J = beglatdyn, endlatdyn
               DO I = 1 ,NX
                  IF ( UQ(I,K,J).GE.ZERO.AND.QT(I,K,J).LT.ZERO ) THEN
                     QT(I,K,J)= ZERO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!=============================== zhh ====================================
!!    write (6,*) 'after filter: '
!!    write (6,*) 'DQ(1,128,1) =', DQ(1,128,1), 'DQ(160,127,1) =', DQ(160,127,1)
!!    write (6,*) 'QT(1,128,1) =', QT(1,128,1), 'QT(160,127,1) =', QT(160,127,1)
!============================ 2007.12.8 =================================
!
! filter DQliq
      DO J = beglatdyn, endlatdyn
         DO K = 1 ,NL
            DO I = 1 ,NX
               UQ(I,K,J) = QTliq(I,K,J)
               VQ(I,K,J) = DQliq(I,K,J)
            ENDDO
         ENDDO
      ENDDO
      DO J = beglatdyn, endlatdyn
         DO K = 1 ,NL
            DO I = 1 ,NX
               DQliq(I,K,J) = UQ(I,K,J) - VQ(I,K,J)
            ENDDO
         ENDDO
      ENDDO
!=============================== zhh ====================================
!!    write (6,*) 'befrore filter: '
!!    write (6,*) 'DQliq(1,128,1) =', DQliq(1,128,1), 'DQliq(160,127,1) =', DQliq(160,127,1)
!============================ 2007.12.8 =================================
!
!     PERFORM SMOOTHING & FILTER IF REQUESTED
      IF ( IFILTQ.NE.0 ) THEN
         DO K = 1 ,NL
            WW(:,beglatdyn:endlatdyn)=DQliq(:,K,beglatdyn:endlatdyn)
            CALL FILT2D( WW,0,1,IBCFFT )
            DQliq(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
            DO J = beglatdyn, endlatdyn
               DO I = 1 ,NX
                  QTliq(I,K,J) = VQ(I,K,J) + DQliq(I,K,J)
               ENDDO
            ENDDO
            DO J = beglatdyn, endlatdyn
               DO I = 1 ,NX
                  IF ( UQ(I,K,J).GE.ZERO.AND.QTliq(I,K,J).LT.ZERO ) THEN
                     QTliq(I,K,J)= ZERO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!     
! filter DQice
      DO J = beglatdyn, endlatdyn
         DO K = 1 ,NL
            DO I = 1 ,NX
               UQ(I,K,J) = QTice(I,K,J)
               VQ(I,K,J) = DQice(I,K,J)
            ENDDO
         ENDDO
      ENDDO
      DO J = beglatdyn, endlatdyn
         DO K = 1 ,NL
            DO I = 1 ,NX
               DQice(I,K,J) = UQ(I,K,J) - VQ(I,K,J)
            ENDDO
         ENDDO
      ENDDO
!=============================== zhh ====================================
!!    write (6,*) 'befrore filter: '
!!    write (6,*) 'DQice(1,128,1) =', DQice(1,128,1), 'DQice(160,127,1) =', DQice(160,127,1)
!============================ 2007.12.8 =================================
!
!     PERFORM SMOOTHING & FILTER IF REQUESTED
      IF ( IFILTQ.NE.0 ) THEN
         DO K = 1 ,NL
            WW(:,beglatdyn:endlatdyn)=DQice(:,K,beglatdyn:endlatdyn)
            CALL FILT2D( WW,0,1,IBCFFT )
            DQice(:,K,beglatdyn:endlatdyn)=WW(:,beglatdyn:endlatdyn)
            DO J = beglatdyn, endlatdyn
               DO I = 1 ,NX
                  QTice(I,K,J) = VQ(I,K,J) + DQice(I,K,J)
               ENDDO
            ENDDO
            DO J = beglatdyn, endlatdyn
               DO I = 1 ,NX
                  IF ( UQ(I,K,J).GE.ZERO.AND.QTice(I,K,J).LT.ZERO ) THEN
                     QTice(I,K,J)= ZERO
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!------------------- PERFORM VERTICAL REDISTRIBUTION TO AVOID NEGATIVE Q-H2O -------------------
      CALL AVNEGQ( QT )
      CALL AVNEGQ( QTliq )
      CALL AVNEGQ( QTice )
!-----------------------------------------------------------------------------------------------
!=============================== zhh =====================================
      if ( abs(DQm) > 1E1 ) then
         write (6,*) 'DQm =', DQm, ' at: I =',Imax, ' J =',Jmax, ' K =',Kmax  
         write (6,*) 'UQ(',Imax,Kmax,Jmax, ') =', UQ(Imax,Kmax,Jmax), 'UQ(',Imax-1,Kmax,Jmax, ') =', UQ(Imax-1,Kmax,Jmax)  
         write (6,*) 'VQ(',Imax,Kmax,Jmax, ') =', VQ(Imax,Kmax,Jmax), 'VQ(',Imax,Kmax,Jmax-1, ') =', VQ(Imax,Kmax,Jmax-1)  
      end if
!!      write (6,*) 'after AVNEGQ: '
!!      write (6,*) 'QT(1,128,1) =', QT(1,128,1), 'QT(160,127,1) =', QT(160,127,1)
!!      write (6,*) 'QTliq(1,128,1) =', QTliq(1,128,1), 'QTliq(160,127,1) =', QTliq(160,127,1)
!!      write (6,*) 'QTice(1,128,1) =', QTice(1,128,1), 'QTice(160,127,1) =', QTice(160,127,1)
!============================ 2007.12.5 ==================================

      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE MPDATA (Q,U,V,P,IORD)
!------------------------------------------------------------------------------------------------
!  Purpose: PERFORM 2-D ADVECTION IN GLOBAL SPHERICAL GEOMETRY
!           WITH THE UNIFORM LAT-LON C-GRID MESH
!------------------------------------------------------------------------------------------------
      use Dyn_const, only: IP, DSNP, DSSP, DTDLN, DTDLT, GC
	  implicit none
!----------------------------------------Arguments-----------------------------------------------
      real(r8), intent(inout) :: Q(NX,NL,beglatdynex:endlatdynex)
	  real(r8), intent(in) :: U(NX,NL,NY)
	  real(r8), intent(in) :: V(NX,NL,NY)
	  real(r8), intent(in) :: P(NX,NY)
      integer,  intent(in) :: IORD     ! times of iterations in water vapor scheme
!----------------------------------Local workspace-----------------------------------------------
	  real(r8) :: VU(NX,NY), VV(NX,NY), FU(NX,NY), FV(NX,NY)
	  real(r8) :: CU(NX,NY), CV(NX,NY), HU(NX,NY), HV(NX,NY)
	  real(r8) :: UU(NX,NY), UV(NX,NY), BU(NX,NY), BD(NX,NY)
	  real(r8) :: X (NX,NY), H (NX,NY), XM(NX,NY), XN(NX,NY)
      real(r8) :: C (NX,2 )
	  real(r8) :: XNP, XSP 
      integer  :: IO, I, J, K
!      
      real(r8) :: DONOR,A,Y1,Y2,VDYF,R,X1,X2,VCOR,B,YY,VC31,VC32,PP,Y,PN
!      real(r8) :: CVMGP
!------------------------------------------------------------------------------------------------
!
!============================== define functions ==============================
!      DONOR(Y1,Y2,A)  = CVMGP(Y1,Y2,A)*A
      VDYF(R,A,X1,X2) = (ABS(A)-A**2/R)*(X2-X1)/(X2+X1+EP)
      VCOR(R,A,B,YY)  = -0.125*A*B*YY/R
      VC31(R,A,YY)    = -(1.0-3.0*ABS(A)/R+2.0*(A/R)**2)*A*YY/3.0
      VC32(R,A,B,YY)  = 0.25*B/R*(ABS(A)-2.0*A**2/R)*YY
!!      PP(Y)           = MAX(dble(0.0),Y)
!!      PN(Y)           = -MIN(dble(0.0),Y)
      PP(Y)           = MAX(0.0,Y)
      PN(Y)           = -MIN(0.0,Y)
!==============================================================================
!
!     GET WORKSHOPS INDEPENDENT OF VERTICAL COORDINATE
      DO J = beglatdyn, endlatdyn
         IF (J.GE.JB.AND.J.LE.JE) THEN
            DO I = IB,IE
               H(I,J) = GC(J) * P(I,J)
            ENDDO
            call period( H(1,J) )
         ELSE   ! at polar
            DO I = 1 ,NX
               H(I,J) = ZERO
            ENDDO
         ENDIF
!
         IF (J.GE.JB.AND.J.LE.JE) THEN
            DO I = IB,IE
               HU(I,J) = HALF*(H(I-1,J)+H(I,J))
            ENDDO
!            call period( H(1,J) )
            call period( HU(1,J) )
            DO I = 1 ,NX
               CU(I,J) = DTDLN(J) * HU(I,J)
            ENDDO
         ELSE   ! at polar
            DO I = 1 ,NX
               HU(I,J) = ZERO
               CU(I,J) = ZERO
            ENDDO
         ENDIF
      ENDDO

#if (defined SPMD)
!jjr      call mpi_move_right( H(1,beglatdyn), H(1,endlatdyn+1),NX)
#endif
      DO J = beglatdyn ,loc_JE
         DO I = IB,IE
            HV(I,J) = HALF*( H(I,J)+H(I,J+1) )
         ENDDO
         call period( HV(1,J) )
         DO I = 1 ,NX
            CV(I,J) = DTDLT(J) * HV(I,J)
         ENDDO
      ENDDO
!
!     START VERTICAL LAYER DO LOOP
!$DOACROSS LOCAL(K,J,I,UU,UV,X,XM,XN,FU,FV,C,XNP,XSP,VU,VV,BU,BD,IO)
      DO K = 1 ,NL
         DO J = 1 ,NY,JE
            DO I = 1 ,NX
               VU(I,J) = ZERO
               FU(I,J) = ZERO
            ENDDO
         ENDDO
!
!     CALCULATE COURANT NUMBERS MULTIPLIED BY H
         DO J = loc_JB, loc_JE
            DO I = 1 ,NX
               UU(I,J)  = CU(I,J) * U(I,K,J)
            END DO
		 END DO
         DO J = beglatdyn ,loc_JE
            DO I = 1 ,NX
               UV(I,J)  = CV(I,J) * V(I,K,J)
            END DO
		 END DO
         DO J = loc_JB,loc_JE
            call period( UU(1,J) )
         END DO
         DO J = beglatdyn ,loc_JE
            call period( UV(1,J) )
         END DO
!
         DO J = loc_JB, loc_JE
            DO I = 1 ,NX
               VU(I,J)  = UU(I,J)
            END DO
		 END DO
         DO J = beglatdyn ,loc_JE
            DO I = 1 ,NX
               VV(I,J)  = UV(I,J)
            END DO
		 END DO
         DO J = beglatdyn ,endlatdyn
            DO I = 1 ,NX
               X (I,J)  = Q (I,K,J)
            END DO
		 END DO
!
!     PREPARE FOR NON-OSSCILATORY OPTION
         IF ( NONOS.EQ.1 ) THEN
#if (defined SPMD)
!jjr            call mpi_move_left(X(1,endlatdyn), X(1,beglatdyn-1),NX)
!jjr            call mpi_move_right( X(1,beglatdyn), X(1,endlatdyn+1),NX)
#endif
            DO J = loc_JB,loc_JE
               DO I = IB,IE
                  XM(I,J) = MAX( X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1) )
                  XN(I,J) = MIN( X(I-1,J),X(I,J),X(I+1,J),X(I,J-1),X(I,J+1) )
               END DO
			END DO
         ENDIF
!
         DO IO = 1 ,IORD
!     ++++++++++++++++++++++++++++++++
!     PREDICTOR STEP : UPSTREAM SCHEME
!     ++++++++++++++++++++++++++++++++
            DO J = loc_JB, loc_JE
               DO I = IB, IE
!b    FU(I,J)  = DONOR( X(I-1,J),X(I,J),VU(I,J) )
                  IF (VU(I,J).GE.0.0E0) THEN
                     FU(I,J) = X(I-1,J)*VU(I,J)
                  ELSE
                     FU(I,J) = X(I,J)*VU(I,J)
                  ENDIF
               END DO
			END DO
            DO J = loc_JB,loc_JE
               call period( FU(1,J) )
            END DO
#if (defined SPMD)
!jjr            call mpi_move_right(GC(  beglatdyn), GC(  endlatdyn+1),1)
!jjr            call mpi_move_right(X (1,beglatdyn), X (1,endlatdyn+1),NX)
#endif
            DO J = beglatdyn ,loc_JE
               DO I = IB, IE
!b    FV(I,J)  = DONOR( X(I,J),X(I,J+1),VV(I,J) )
                  IF (VV(I,J).GE.0.0E0) THEN
                     FV(I,J) = X(I,J)*VV(I,J)
                  ELSE
                     FV(I,J) = X(I,J+1)*VV(I,J)
                  ENDIF
               END DO
               call period( FV(1,J) )   !zhh 2007.12.6
			END DO
!=============================== zhh ====================================
            DO I = 1, NX
               FU(I,1)  = ZERO
               FU(I,NY) = ZERO
               FV(I,NY) = ZERO
			END DO
!============================ 2007.12.6 =================================
#if (defined SPMD)
!jjr            call mpi_move_left(FV(1,endlatdyn), FV(1,beglatdyn-1),NX)
#endif
            DO J = loc_JB, loc_JE
               DO I = IB, IE
                  X(I,J) = X(I,J)-( FU(I+1,J)-FU(I,J)+FV(I,J)-FV(I,J-1) )/H(I,J)
               END DO
			END DO
!     B.C. BY THE SPHERICAL CYCLICITY & MASS CONSERVATION RESTRICT
            DO I = IB, IE
               C(I,1)    = FV(I,1) / H(I,2)
               C(I,2)    = FV(I,JE) / H(I,JE)
			END DO
            XNP      = ZERO
            XSP      = ZERO
            DO I = IB, IE
               XNP   = XNP + C(I,1)
               XSP   = XSP + C(I,2)
            END DO
            XNP      = X(IB,1) - XNP*DSNP
            XSP      = X(IB,NY) + XSP*DSSP
            DO I = IB, IE
               X(I,1)   = XNP
               X(I,NY)  = XSP
            END DO
            DO J = beglatdyn, endlatdyn
               call period( X(1,J) )
            END DO
!
            IF ( IO.EQ.IORD ) GOTO 700
!     ++++++++++++++++++++++++++++++++++++++
!     CORRECTOR STEP : ANTI-DIFFUSION SCHEME
!     ++++++++++++++++++++++++++++++++++++++
            DO J = loc_JB, loc_JE
               DO I = 1 ,NX
                  FU(I,J)  = VU(I,J)
               END DO
			END DO
            DO J = beglatdyn ,loc_JE
               DO I = 1 ,NX
                  FV(I,J)  = VV(I,J)
               END DO
			END DO
!
!     CALCULATE THE  PSEUDO VELOCITIES  LONGITUDINAL DIRECTION
#if (defined SPMD)
!jjr            call mpi_move_left(X(1,endlatdyn),   X(1,beglatdyn-1),NX)
!jjr            call mpi_move_left(FV(1,endlatdyn),  FV(1,beglatdyn-1),NX)
!jjr            call mpi_move_right( X(1,beglatdyn-1), X(1,endlatdyn+1),NX*2)
#endif
            DO J = loc_JB, loc_JE
               DO I = IB, IE
                  VU(I,J) = VDYF( HU(I,J), FU(I,J), X(I-1,J), X(I,J) )                         &
                          + VCOR( HU(I,J), FU(I,J), (FV(I-1,J-1)+FV(I-1,J)+FV(I,J)+FV(I,J-1)), &
                                  (X(I-1,J+1)+X(I,J+1)-X(I-1,J-1)-X(I,J-1))                    &
                                / (EP+X(I-1,J+1)+X(I,J+1)+X(I-1,J-1)+X(I,J-1)) )
               END DO
			END DO
!               LATITUDINAL  DIRECTION
!     LATITUDINAL  DIRECTION
#if (defined SPMD)
!jjr            call mpi_move_right(X(1,beglatdyn), X(1,endlatdyn+1),NX)
!jjr            call mpi_move_right(FU(1,beglatdyn),FU(1,endlatdyn+1),NX)
#endif
            DO J = beglatdyn, loc_JE
               DO I = IB, IE
                  VV(I,J) = VDYF( HV(I,J), FV(I,J), X(I,J), X(I,J+1) )                         &
                          + VCOR( HV(I,J), FV(I,J), (FU(I,J)+FU(I,J+1)+FU(I+1,J)+FU(I+1,J+1)), &      
                                  (X(I+1,J)+X(I+1,J+1)-X(I-1,J)-X(I-1,J+1))                    &
                                / (EP+X(I+1,J)+X(I+1,J+1)+X(I-1,J)+X(I-1,J+1)) )
               END DO
			END DO
!
!     ADD THE THIRD ORDER CORRECTION IF REQUESTED
            IF ( ISOR.EQ.3 ) THEN
!     LONGITUDINAL DIRECTION
#if (defined SPMD)
!jjr               call mpi_move_left(X(1,endlatdyn), X(1,beglatdyn-1),NX)
!jjr               call mpi_move_left(FV(1,endlatdyn),FV(1,beglatdyn-1),NX)
!jjr!               call mpi_move_right( X(1,beglatdyn-1), X(1,endlatdyn+1),NX*2)
#endif
               DO J = loc_JB,loc_JE
                  DO I = IB ,IE
!============================ 2007.12.8 =================================
                     VU(I,J) = VU(I,J) + VC31( HU(I,J), FU(I,J),                           &
                                               (X(I-2,J)+X(I+1,J)-X(I-1,J)-X(I,J))         &
                                             / (EP+X(I-2,J)+X(I+1,J)+X(I-1,J)+X(I,J)) )
                  END DO
			   END DO
               DO J = loc_JB,loc_JE
                  DO I = IB,IE
                     VU(I,J) = VU(I,J) + VC32( HU(I,J), FU(I,J),                           &
                                               (FV(I-1,J-1)+FV(I-1,J)+FV(I,J-1)+FV(I,J)),  &
                                               (X(I,J+1)-X(I,J-1)-X(I-1,J+1)+X(I-1,J-1))   &
                                             / (EP+X(I,J+1)+X(I,J-1)+X(I-1,J+1)+X(I-1,J-1)) )
                  END DO
			   END DO
!               LATITUDINAL  DIRECTION
               DO J = loc_JB,min(JE-1,endlatdyn)
                  DO I = IB,IE
                     VV(I,J) = VV(I,J) + VC31( HV(I,J), FV(I,J),                           &
                                               (X(I,J-1)+X(I,J+2)-X(I,J)-X(I,J+1) )        &
                                             / (EP+X(I,J-1)+X(I,J+2)+X(I,J)+X(I,J+1)) )
                  END DO
			   END DO
!     B.C. BY THE SPHERICAL CYCLICITY
               DO I = IB,IE
                  C(I,1)  = X(IP(I),2)
                  C(I,2)  = X(IP(I),JE)
               END DO
               DO I = IB,IE
                  VV(I,1) = VV(I,1) + VC31( HV(I,1), FV(I,1),                              &
				                            (C(I,1)+X(I,3)-X(I,1)-X(I,2) )                 &   
                                          / (EP+C(I,1)+X(I,3)+X(I,1)+X(I,2)) )
                  VV(I,JE)= VV(I,JE) + VC31( HV(I,JE), FV(I,JE),                           &
                                             (X(I,JE-1)+C(I,2)-X(I,JE)-X(I,NY))            &
                                           / (EP+X(I,JE-1)+C(I,2)+X(I,JE)+X(I,NY)) )
               END DO
               DO J = beglatdyn ,loc_JE
                  DO I = IB,IE
                     VV(I,J)  = VV(I,J) + VC32( HV(I,J), FV(I,J),                          &
                                                (FU(I,J)+FU(I+1,J)+FU(I+1,J+1)+FU(I,J+1)), &
                                                (X(I+1,J+1)-X(I-1,J+1)-X(I+1,J)+X(I-1,J))  &
                                              / (EP+X(I+1,J+1)+X(I-1,J+1)+X(I+1,J)+X(I-1,J)) )
                  END DO
			   END DO
            ENDIF
!
            DO J = loc_JB,loc_JE
               DO I = IB,IE
                  VU(I,J)  = SIGN(ONE,VU(I,J))*MIN(ABS(UU(I,J)),ABS(VU(I,J)))
               END DO
			END DO
            DO J = beglatdyn ,loc_JE
               DO I = IB,IE
                  VV(I,J)  = SIGN(ONE,VV(I,J))*MIN(ABS(UV(I,J)),ABS(VV(I,J)))
               END DO
			END DO
!     B.C. BY THE SPHERICAL CYCLICITY
            DO J = loc_JB,loc_JE
               call period( VU(1,J) )
            END DO
            DO J = beglatdyn ,loc_JE
               call period( VV(1,J) )
            END DO
!
!     PERFORM THE NON-OSSCILATORY OPTION
            IF ( NONOS.EQ.1 ) THEN
               DO J = loc_JB,loc_JE
                  DO I = IB,IE
                     XM(I,J) = MAX( X(I-1,J), X(I,J), X(I+1,J), X(I,J-1), X(I,J+1), XM(I,J) )      
                     XN(I,J) = MIN( X(I-1,J), X(I,J), X(I+1,J), X(I,J-1), X(I,J+1), XN(I,J) )        
                  END DO
			   END DO
!
               DO J = loc_JB,loc_JE
                  DO I = IB,IE
!b    FU(I,J)  = DONOR( X(I-1,J),X(I,J),VU(I,J) )
                     IF (VU(I,J).GE.0.0E0) THEN
                        FU(I,J) = X(I-1,J)*VU(I,J)
                     ELSE
                        FU(I,J) = X(I,J)*VU(I,J)
                     ENDIF
                  END DO
			   END DO
               DO J = loc_JB,loc_JE
                  call period( FU(1,J) )
               END DO
               DO J = beglatdyn ,loc_JE
                  DO I = IB,IE
!b    FV(I,J)  = DONOR( X(I,J),X(I,J+1),VV(I,J) )
                     IF (VV(I,J).GE.0.0E0) THEN
                        FV(I,J) = X(I,J)*VV(I,J)
                     ELSE
                        FV(I,J) = X(I,J+1)*VV(I,J)
                     ENDIF
                  END DO
			   END DO
               DO J = loc_JB,loc_JE
                  DO I = IB,IE
                     BU(I,J) = (XM(I,J)-X(I,J))*H(I,J) /                                   &
                               (PN(FU(I+1,J))+PP(FU(I,J))+PN(FV(I,J))+PP(FV(I,J-1))+EP)
                     BD(I,J) = (X(I,J)-XN(I,J))*H(I,J) /                                   &
                               (PP(FU(I+1,J))+PN(FU(I,J))+PP(FV(I,J))+PN(FV(I,J-1))+EP)
                  END DO
			   END DO
               DO J = loc_JB,loc_JE
			      call period( BU(1,J) )
				  call period( BD(1,J) )
               END DO
!     KEEP IN MIND THAT H ARE ZERO AT POLE
               DO I = 1 ,IE
                  BU(I,01) = ZERO
                  BU(I,NY) = ZERO
                  BD(I,01) = ZERO
                  BD(I,NY) = ZERO
               END DO
               DO J = loc_JB,loc_JE
                  DO I = IB,IE
                     VU(I,J)  = PP( VU(I,J) ) * MIN(ONE,BD(I-1,J),BU(I,J))                 &
                              - PN( VU(I,J) ) * MIN(ONE,BU(I-1,J),BD(I,J))
                  END DO
			   END DO
#if (defined SPMD)
!jjr               call mpi_move_right(BD(1,beglatdyn), BD(1,endlatdyn+1),NX)
!jjr               call mpi_move_right(BU(1,beglatdyn), BU(1,endlatdyn+1),NX)
#endif
               DO J = beglatdyn ,loc_JE
                  DO I = IB,IE
                     VV(I,J)  = PP( VV(I,J) ) * MIN(ONE,BD(I,J),BU(I,J+1))                 &
                              - PN( VV(I,J) ) * MIN(ONE,BU(I,J),BD(I,J+1))
                  END DO
			   END DO
!     B.C. BY THE SPHERICAL CYCLICITY
               DO J = loc_JB,loc_JE
			      call period( VU(1,J) )
               END DO
               DO J = beglatdyn ,loc_JE
			      call period( VV(1,J) )
               END DO
            ENDIF
700      END DO         ! end for IO = 1 ,IORD
!
!     UPDATE THE PREDICTED FIELD
         DO J = beglatdyn, endlatdyn
            DO I = 1 ,NX
               Q(I,K,J) = X(I,J)
            END DO
	     END DO
      END DO           ! end for K = 1 ,NL
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE VPDATA(Q,W,IORD)
!------------------------------------------------------------------------------------------------
!  Purpose: PERFORM 1-D ADVECTION IN THE VERTICAL DIRECTION
!           WITH A NON-UNIFORM SIGMA C-GRID MESH
!------------------------------------------------------------------------------------------------
      use Dyn_const, only: DTDSG 
	  implicit none
!----------------------------------------Arguments-----------------------------------------------
      real(r8), intent(inout) :: Q(NX,NL,beglatdynex:endlatdynex)
	  real(r8), intent(in) :: W(NX,NL,NY)
      integer,  intent(in) :: IORD     ! times of iterations in water vapor scheme
!----------------------------------Local workspace-----------------------------------------------
	  real(r8) :: C (NX,2 ), HW(NL) , HS(NL) , Z (NL)
	  real(r8) :: VW(NZ)   , ZM(NL) , ZN(NL) , FW(NZ)
	  real(r8) :: UW(NZ)   , BUW(NL), BDW(NL)
!!	  real(r8) :: VDYF, VC31, PN, PP, 
      real(r8) :: A,VDYF,R,X1,X2,YY,VC31,PP,Y,PN
      integer  :: IO, IS, IT
	  integer  :: I, J, K
!------------------------------------------------------------------------------------------------
!============================== define functions ==============================
!      DONOR(Y1,Y2,A)  = CVMGP(Y1,Y2,A)*A
      VDYF(R,A,X1,X2) = (ABS(A)-A**2/R)*(X2-X1)/(X2+X1+EP)
      VC31(R,A,YY)    = -(1.0-3.0*ABS(A)/R+2.0*(A/R)**2)*A*YY/3.0
!!      PP(Y)           = MAX(dble(0.0),Y)
!!      PN(Y)           = -MIN(dble(0.0),Y)
      PP(Y)           = MAX(0.0,Y)
      PN(Y)           = -MIN(0.0,Y)
!==============================================================================

!     GET WORKSHOPS INDEPENDENT OF HORIZONTAL MESH
      DO K = 1 ,NL
         HS(K) = ONE  / DTDSG(K)
      END DO
      DO K = 2 ,NL
         HW(K) = HALF * ( HS(K)+HS(K-1) )
      END DO
!
!     START HORIZONTAL GRID DO LOOP
!     KEEP IN MIND THAT Z ARE THE SAME FOR ALL I AT POLE
!$DOACROSS LOCAL(J,IS,IT,I,K,Z,UW,VW,ZM,ZN,FW,BUW,BDW,IO)
      DO J = beglatdyn, endlatdyn
         FW(01)   = ZERO
         FW(NZ)   = ZERO
         IF (J.EQ.1) THEN
            IS = NX          !!???????
            IT = NX          !!???????
         ELSE IF(J.EQ.NY) THEN
            IS = 1           !!???????
            IT = 1           !!???????
         ELSE
            IS = 1           !!???????
            IT = NX          !!???????
         ENDIF
         DO I = IS,IT
            DO K = 1 ,NL
               Z (K) = Q(I,K,J)
               UW(K) = W(I,K,J)
               VW(K) = UW(K)
            END DO
!
!     PREPARE FOR NON-OSSCILATORY OPTION
            IF ( NONOS.EQ.1 ) THEN
               DO K = 2 ,NM
                  ZM(K) = MAX( Z(K-1),Z(K),Z(K+1) )
                  ZN(K) = MIN( Z(K-1),Z(K),Z(K+1) )
               END DO
               ZM(1)    = MAX( Z(1),Z(2) )
               ZN(1)    = MIN( Z(1),Z(2) )
               ZM(NL)   = MAX( Z(NM),Z(NL) )
               ZN(NL)   = MIN( Z(NM),Z(NL) )
            ENDIF
!
            DO IO = 1 ,IORD
!     ++++++++++++++++++++++++++++++++
!     PREDICTOR STEP : UPSTREAM SCHEME
!     ++++++++++++++++++++++++++++++++
               DO K = 2 ,NL
!b    FW(K)    = DONOR( Z(K-1),Z(K),VW(K) )
                  IF (VW(K).GE.0.0D0) THEN
                      FW(K) = Z(K-1)*VW(K)
                  ELSE
                      FW(K) = Z(K)*VW(K)
                  ENDIF
               END DO
               DO K = 1 ,NL
                  Z(K)      = Z(K)-(FW(K+1)-FW(K))*DTDSG(K)
               END DO
!
               IF ( IO.EQ.IORD ) GOTO 980
!     ++++++++++++++++++++++++++++++++++++++
!     CORRECTOR STEP : ANTI-DIFFUSION SCHEME
!     ++++++++++++++++++++++++++++++++++++++
               DO K = 2 ,NL
                  FW(K)    = VW(K)
               END DO
!
!     CALCULATE THE  PSEUDO VELOCITIES
               DO K = 2 ,NL
                  VW(K)    = VDYF( HW(K),FW(K),Z(K-1),Z(K) )
               END DO
!     ADD THE THIRD ORDER CORRECTION IF REQUESTED
               IF ( ISOR.EQ.3 ) THEN
                  DO K = 3 ,NM
                     VW(K) = VW(K) + VC31( HW(K), FW(K), (Z(K-2)+Z(K+1)-Z(K-1)-Z(K))       & 
                                                       / (EP+Z(K-2)+Z(K+1)+Z(K-1)+Z(K)) )
                  END DO
!     ASSUME CONSTANT Z ABOVE K=1 & BELOW K=NL
                  VW(2)  = VW(2) + VC31( HW(2), FW(2), (Z(3)-Z(2))/(EP+Z(1)+Z(3)+Z(1)+Z(2)) )                 
                  VW(NL) = VW(NL) + VC31( HW(NL), FW(NL), (Z(NL-2)-Z(NL-1))                &
                                                       / (EP+Z(NL-2)+Z(NL)+Z(NL-1)+Z(NL)) )
               ENDIF
!
               DO K = 2 ,NL
                  VW(K)  = SIGN(ONE,VW(K))*MIN(ABS(UW(K)),ABS(VW(K)))
               END DO
!
!     PERFORM THE NON-OSSCILATORY OPTION
               IF ( NONOS.EQ.1 ) THEN
                  DO K = 2 ,NM
                     ZM(K) = MAX( Z(K-1),Z(K),Z(K+1),ZM(K) )
                     ZN(K) = MIN( Z(K-1),Z(K),Z(K+1),ZN(K) )
                  END DO
                  ZM(1)    = MAX( Z(1),Z(2),ZM(1) )
                  ZN(1)    = MIN( Z(1),Z(2),ZN(1) )
                  ZM(NL)   = MAX( Z(NM),Z(NL),ZM(NL) )
                  ZN(NL)   = MIN( Z(NM),Z(NL),ZN(NL) )
!
                  DO K = 2 ,NL
!b    FW(K)    = DONOR( Z(K-1),Z(K),VW(K) )
                     IF (VW(K).GE.0.0D0) THEN
                        FW(K) = Z(K-1)*VW(K)
                     ELSE
                        FW(K) = Z(K)*VW(K)
                     ENDIF
                  END DO
                  DO K = 1 ,NL
                     BUW(K)   = (ZM(K)-Z(K))*HS(K) / (PN(FW(K+1))+PP(FW(K))+EP)
                     BDW(K)   = (Z(K)-ZN(K))*HS(K) / (PP(FW(K+1))+PN(FW(K))+EP)
                  END DO
                  DO K = 2 ,NL
                     VW(K)    = PP( VW(K) ) * MIN(ONE,BDW(K-1),BUW(K))                     &
                              - PN( VW(K) ) * MIN(ONE,BUW(K-1),BDW(K))
                  END DO
               ENDIF
980         END DO       ! end IO = 1 ,IORD
!
!     UPDATE THE PREDICTED FIELD
            DO K = 1 ,NL
               Q(I,K,J) = Z(K)
            END DO
		 END DO      ! end I = IS,IT (1, NX)
!
         IF (J.EQ.1) THEN
            DO K = 1 ,NL
               DO I = 1 ,NX
                  Q(I,K,J)= Q(NX,K,J)        !!???????
               ENDDO
            ENDDO
         ELSE IF(J.EQ.1.OR.J.EQ.NY) THEN
            DO K = 1 ,NL
               DO I = 1 ,NX
                  Q(I,K,J) = Q(1,K,J)        !!???????
               ENDDO
            ENDDO
         ENDIF
      END DO       ! end J = 1 ,NY
      RETURN
   END SUBROUTINE

!================================================================================================
   SUBROUTINE AVNEGQ( Q )
!------------------------------------------------------------------------------------------------
!  Purpose: AVOID   NEGATIVE MIXING RATIO
!------------------------------------------------------------------------------------------------
      use Dyn_const, only: DSGHL 
	  implicit none
!----------------------------------------Arguments-----------------------------------------------
      real(r8), intent(inout) :: Q(NX,NL,beglatdynex:endlatdynex)
!--------------------------------------Local workspace-------------------------------------------
	  real(r8) :: QR(NL)
	  real(r8) :: QI
      integer  :: I, J, K
!------------------------------------------------------------------------------------------------

!!      DO J = loc_JB, loc_JE
      DO J = beglatdyn, endlatdyn
         DO I = IB, IE
            DO K = 1 ,NL
               QR(K)      = Q(I,K,J)
            END DO
            DO K = 2 ,NL
               QI         = QR(K-1)
               IF ( QI.LT.ZERO ) THEN
                  QR(K-1) = ZERO
                  QR(K )  = QR(K ) + QI*DSGHL(K)
               ENDIF
            END DO
            IF ( QR(NL).LT.ZERO ) QR(NL) = ZERO
            DO K  = 1 ,NL
               Q(I,K,J) = QR(K)
            END DO
         END DO
	  END DO
!     
!!      DO J = 1 ,NY,JE !at the north and south polar
!!         DO K = 1 ,NL
!!            QR(K)       = Q(IB,K,J)
!!		 END DO
!!         DO K = 2 ,NL
!!            QI          = QR(K-1)
!!            IF ( QI.LT.ZERO ) THEN
!!               QR(K-1)  = ZERO
!!               QR(K )   = QR(K ) + QI*DSGHL(K)
!!            ENDIF
!!         END DO
!!         IF ( QR(NL).LT.ZERO ) QR(NL) = ZERO
!!         DO K = 1 ,NL
!!            DO I = IB,IE
!!               Q(I,K,J) = QR(K)
!!            END DO
!!         END DO
!!	  END DO

      DO J = beglatdyn, endlatdyn
         DO K = 1 ,NL
            call period( Q(1,K,J) )
         END DO
	  END DO
!
      RETURN
   END SUBROUTINE

end module
