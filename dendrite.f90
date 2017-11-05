!=============================================================================
!
! Dendric solidification analysis (demonstration)
! Programed by kruho   2012/10/15
! See. 高木知弘・山中晃徳, "フェーズフィールド法", 養賢堂, (2012), pp.61-70.
!
!=============================================================================
!=============================================================================
! Output format module
!=============================================================================

MODULE FORMAT_FNC

  IMPLICIT NONE

  CONTAINS

    INTEGER FUNCTION DIGIT(NUM)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NUM
      IF ( NUM == 0 ) THEN
         DIGIT = 1
      ELSE
         DIGIT = INT(LOG10(DBLE(NUM))) + 1
      END IF
      RETURN
    END FUNCTION DIGIT

    FUNCTION FMT(HEAD,NUM)
      IMPLICIT NONE
      CHARACTER :: FMT*30
      CHARACTER, INTENT(IN) :: HEAD*(*)
      INTEGER, INTENT(IN) :: NUM
      FMT = '("'//TRIM(HEAD)//'_",I1,".res")'
      WRITE(FMT(LEN(TRIM(HEAD))+7:LEN(TRIM(HEAD))+7),'(I1)') DIGIT(NUM)
      RETURN
    END FUNCTION FMT

END MODULE FORMAT_FNC

!=============================================================================
!=============================================================================

MODULE OUT_INFO

  IMPLICIT NONE
  LOGICAL :: OUTPUT_FLAG
  CHARACTER :: FNOUT1*20, FNOUT2*20
  INTEGER :: JGRAPH=-1, NOUT=100, JINCRE, NINCRE=100000
  CHARACTER(LEN=10), PARAMETER :: FNRES = './result/'
  INTEGER, PARAMETER :: IFRED = 5
  INTEGER, PARAMETER :: IFOUT1 = 10
  INTEGER, PARAMETER :: IFOUT2 = 11

END MODULE OUT_INFO

!=============================================================================
!=============================================================================

MODULE CONSTANT

  IMPLICIT NONE
  REAL(8) :: THEDIF, DELTAI, CNST_L, CNST_B, A_REFR, WALL_E, PF_MOB, DELTAT
  INTEGER, PARAMETER :: NXSIZE = 501
  INTEGER, PARAMETER :: NYSIZE = 751
  REAL(8), PARAMETER :: H_COND = 84.01D0
  REAL(8), PARAMETER :: SPHEAT = 5.42D6
  REAL(8), PARAMETER :: HEAT_L = 2.350D9
  REAL(8), PARAMETER :: TEMP_M = 1728.D0
  REAL(8), PARAMETER :: COEF_M = 2.0D0
  REAL(8), PARAMETER :: ENGY_I = 0.37D0
  REAL(8), PARAMETER :: THETAZ = 0.D0
  REAL(8), PARAMETER :: DELTAX = 20.D-9
  REAL(8), PARAMETER :: DELTAY = 20.D-9
  REAL(8), PARAMETER :: AISOMD = 4.D0
  REAL(8), PARAMETER :: X_SIZE = 10.D-6
  REAL(8), PARAMETER :: Y_SIZE = 15.D-6
  REAL(8), PARAMETER :: AISO_I = 0.03D0
  REAL(8), PARAMETER :: TEMP_0 = 1511.2D0

END MODULE CONSTANT

!=============================================================================
!=============================================================================

PROGRAM DENDRITE

  !$ USE OMP_LIB
  USE OUT_INFO
  USE CONSTANT
  IMPLICIT NONE
  INTEGER :: LPREV=1, LNEXT=2, IX, IY, IXM, IXP, IYM, IYP
  REAL(8) ::  X1, X2, F_A, F_DA, F_B, P_M0, P_0M, P_0P, P_P0, THETA, &
       & P_00, T_00, T_M0, T_P0, T_0M, T_0P, TMP, DXINV, DYINV, DXINV2, &
       & DYINV2, DX2INV, DY2INV, &
       & DPDX, DPDY, DPDX2, DPDY2, A_00, A_M0, A_P0, A_0M, A_0P, &
       & A2_00, A2_M0, A2_P0, A2_0M, A2_0P, DADX, DADY, DA2DX, DA2DY, DA_00, &
       & DA_M0, DA_P0, DA_0M, DA_0P, DAQDX, DAQDY, DTDX2, DTDY2, PRAT, TRAT
  REAL(8), ALLOCATABLE :: VPHASE(:,:,:), VPRATE(:,:), TEMPRT(:,:,:), COEF_A(:,:), DA_DTH(:,:)
  REAL(8), PARAMETER :: PIH = 1.57079632679490D0

  CALL OMP_SET_NUM_THREADS ( 4 )

  THEDIF = H_COND/SPHEAT
  DELTAI = 4.D0*DELTAX
  CNST_L = 0.1D0
  CNST_B = 2.D0*ATANH(1.D0-2.D0*CNST_L)
  A_REFR = SQRT(3.D0*DELTAI*ENGY_I/CNST_B)
  WALL_E = 6.D0*ENGY_I*CNST_B/DELTAI
  PF_MOB = CNST_B*TEMP_M*COEF_M/(3.D0*DELTAI*HEAT_L)
  X1 = DELTAX**2.D0/(5.D0*PF_MOB*A_REFR**2.D0)
  X2 = DELTAX**2.D0/(5.D0*THEDIF)
  DELTAT = MIN(X1,X2)

  DXINV = 1.D0/DELTAX
  DYINV = 1.D0/DELTAY
  DXINV2 = 0.5D0*DXINV
  DYINV2 = 0.5D0*DYINV
  DX2INV = DXINV**2.D0
  DY2INV = DYINV**2.D0

  ALLOCATE ( VPHASE(NXSIZE,NYSIZE,2), &
           & VPRATE(NXSIZE,NYSIZE), &
           & COEF_A(NXSIZE,NYSIZE), &
           & DA_DTH(NXSIZE,NYSIZE), &
           & TEMPRT(NXSIZE,NYSIZE,2) )

  DO IY = 1, NYSIZE
     DO IX = 1, NXSIZE
        TMP = (DBLE(IX)*DELTAX-0.5D0*X_SIZE)**2.D0 &
          & + (DBLE(IY)*DELTAY-0.0D0*Y_SIZE)**2.D0
        TMP = SQRT(TMP)
        P_00 = 0.5D0*(1.D0-TANH(SQRT(2.D0*WALL_E)/(2.D0*A_REFR)*(TMP-2.D0*DELTAX)))
        VPHASE(IX,IY,LPREV) = P_00
        TEMPRT(IX,IY,LPREV) = TEMP_0 + P_00*(TEMP_M-TEMP_0)
        VPRATE(IX,IY) = 0.D0
     END DO
  END DO

  CALL OUTPUT_VALUE ( VPHASE(:,:,LPREV), TEMPRT(:,:,LPREV) )

  DO JINCRE = 1, NINCRE

     OUTPUT_FLAG = .FALSE.
     IF ( MOD(JINCRE,NINCRE/NOUT) == 0 ) OUTPUT_FLAG = .TRUE.

     !$OMP PARALLEL PRIVATE ( IY, IYM, IYP, IX, IXM, IXP, &
     !$OMP & P_M0, P_0M, P_0P, P_P0, DPDX, DPDY, THETA )
     !$OMP DO
     DO IY = 1, NYSIZE
        IYM = IY - 1
        IYP = IY + 1
        IF ( IYM == 0 ) IYM = 1
        IF ( IYP == NYSIZE+1 ) IYP = NYSIZE
        DO IX = 1, NXSIZE
           IXM = IX - 1
           IXP = IX + 1
           IF ( IXM == 0 ) IXM = 1
           IF ( IXP == NXSIZE+1 ) IXP = NXSIZE
           P_0M = VPHASE(IX ,IYM,LPREV)
           P_M0 = VPHASE(IXM,IY ,LPREV)
           P_0P = VPHASE(IX ,IYP,LPREV)
           P_P0 = VPHASE(IXP,IY ,LPREV)
           DPDX = DXINV2*(P_P0-P_M0)
           DPDY = DYINV2*(P_0P-P_0M)
           IF ( DPDX == 0.D0 ) THEN
              THETA = PIH
           ELSE
              THETA = ATAN(DPDY/DPDX)
           END IF
           COEF_A(IX,IY) = F_A(THETA)
           DA_DTH(IX,IY) = F_DA(THETA)
        END DO
     END DO
     !$OMP END DO
     !$OMP END PARALLEL

     !$OMP PARALLEL PRIVATE ( IY, IYM, IYP, IX, IXM, IXP, &
     !$OMP & P_00, P_0M, P_M0, P_0P, P_P0, DPDX, DPDY, DPDX2, DPDY2, &
     !$OMP & A_00, A_M0, A_P0, A_0M, A_0P, A2_00, &
     !$OMP & A2_M0, A2_P0, A2_0M, A2_0P, DADX, DADY, DA2DX, DA2DY, DA_00, &
     !$OMP & DA_M0, DA_P0, DA_0M, DA_0P, DAQDX, DAQDY, PRAT, &
     !$OMP & T_00, T_M0, T_P0, T_0M, T_0P, DTDX2, DTDY2, TRAT )
     !$OMP DO
     DO IY = 1, NYSIZE
        IYM = IY - 1
        IYP = IY + 1
        IF ( IYM == 0 ) IYM = 1
        IF ( IYP == NYSIZE+1 ) IYP = NYSIZE
        DO IX = 1, NXSIZE
           IXM = IX - 1
           IXP = IX + 1
           IF ( IXM == 0 ) IXM = 1
           IF ( IXP == NXSIZE+1 ) IXP = NXSIZE

           P_00 = VPHASE(IX ,IY ,LPREV)
           P_0M = VPHASE(IX ,IYM,LPREV)
           P_M0 = VPHASE(IXM,IY ,LPREV)
           P_0P = VPHASE(IX ,IYP,LPREV)
           P_P0 = VPHASE(IXP,IY ,LPREV)

           A_00 = COEF_A(IX ,IY )
           A_M0 = COEF_A(IXM,IY )
           A_P0 = COEF_A(IXP,IY )
           A_0M = COEF_A(IX ,IYM)
           A_0P = COEF_A(IX ,IYP)

           DA_00 = DA_DTH(IX ,IY )
           DA_M0 = DA_DTH(IXM,IY )
           DA_P0 = DA_DTH(IXP,IY )
           DA_0M = DA_DTH(IX ,IYM)
           DA_0P = DA_DTH(IX ,IYP)

           T_00 = TEMPRT(IX ,IY ,LPREV)
           T_M0 = TEMPRT(IXM,IY ,LPREV)
           T_P0 = TEMPRT(IXP,IY ,LPREV)
           T_0M = TEMPRT(IX ,IYM,LPREV)
           T_0P = TEMPRT(IX ,IYP,LPREV)

           DPDX = DXINV2*(P_P0-P_M0)
           DPDY = DYINV2*(P_0P-P_0M)

           DPDX2 = DX2INV*(P_P0-2.D0*P_00+P_M0)
           DPDY2 = DY2INV*(P_0M-2.D0*P_00+P_0P)

           A2_00 = A_00*A_00
           A2_M0 = A_M0*A_M0
           A2_P0 = A_P0*A_P0
           A2_0M = A_0M*A_0M
           A2_0P = A_0P*A_0P

           DADX = DXINV2*(A_P0-A_M0)
           DADY = DYINV2*(A_0P-A_0M)

           DA2DX = DXINV2*(A2_P0-A2_M0)
           DA2DY = DYINV2*(A2_0P-A2_0M)

           DAQDX = DXINV2*(DA_P0-DA_M0)
           DAQDY = DYINV2*(DA_0P-DA_0M)

           PRAT = DA2DX*DPDX + DA2DY*DPDY        &
                & + A2_00*(DPDX2+DPDY2)          &
                & + DA_00*(DADY*DPDX-DADX*DPDY)  &
                & + A_00*(DAQDY*DPDX-DAQDX*DPDY) &
                & + F_B(P_00,T_00)

           PRAT = PF_MOB*PRAT

           DTDX2 = DX2INV*(T_P0-2.D0*T_00+T_M0)
           DTDY2 = DY2INV*(T_0P-2.D0*T_00+T_0M)

           TRAT = THEDIF*(DTDX2+DTDY2) &
                & + 30.D0*P_00**2.D0*(1.D0-P_00)**2.D0 &
                & * HEAT_L/SPHEAT*VPRATE(IX,IY)

           TEMPRT(IX,IY,LNEXT) = T_00 + TRAT*DELTAT

           VPRATE(IX,IY) = PRAT
           VPHASE(IX,IY,LNEXT) = P_00 + PRAT*DELTAT

        END DO
     END DO
     !$OMP END DO
     !$OMP END PARALLEL

     ! Output
     IF ( OUTPUT_FLAG ) THEN
        CALL OUTPUT_VALUE ( VPHASE(:,:,LNEXT), TEMPRT(:,:,LNEXT) )
     END IF

     ! Update layer
     LPREV = 3 - LPREV
     LNEXT = 3 - LNEXT

  END DO

  STOP

END PROGRAM DENDRITE

!=============================================================================

REAL(8) FUNCTION F_A ( Q )

  USE CONSTANT
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: Q

  F_A = A_REFR*(1.D0+AISO_I*COS(AISOMD*(Q-THETAZ)))

  RETURN

END FUNCTION F_A

!=============================================================================

REAL(8) FUNCTION F_DA ( Q )

  USE CONSTANT
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: Q

  F_DA = -AISOMD*A_REFR*AISO_I*SIN(AISOMD*(Q-THETAZ))

  RETURN

END FUNCTION F_DA

!=============================================================================

REAL(8) FUNCTION F_B ( P, T )

  USE CONSTANT
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: P, T
  REAL(8) :: FLUC

  CALL RANDOM_NUMBER ( FLUC )
  FLUC = 0.2D0*(FLUC-0.5D0)
  F_B = 15.D0/(2.D0*WALL_E) * HEAT_L*(T-TEMP_M)/TEMP_M * P*(1.D0-P)
  F_B = 4.D0*WALL_E*P*(1.D0-P)*( P - 0.5D0 - F_B + FLUC )

  RETURN

END FUNCTION F_B

!=============================================================================
!=============================================================================

SUBROUTINE OUTPUT_VALUE ( P, T )

  USE FORMAT_FNC
  USE OUT_INFO
  USE CONSTANT
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: P(NXSIZE,NYSIZE), T(NXSIZE,NYSIZE)
  INTEGER :: IX, IY
  REAL(8) :: TMP

  JGRAPH = JGRAPH + 1
  WRITE(FNOUT1,FMT('phi',JGRAPH)) JGRAPH
  WRITE(FNOUT2,FMT('thrm',JGRAPH)) JGRAPH

  OPEN ( IFOUT1, STATUS = 'REPLACE', FILE = TRIM(FNRES)//FNOUT1 )
  OPEN ( IFOUT2, STATUS = 'REPLACE', FILE = TRIM(FNRES)//FNOUT2 )
 
  DO IY = 1, NYSIZE
     DO IX = 1, NXSIZE
        TMP = P(IX,IY)
        IF ( ABS(P(IX,IY)) < 1.D-20 ) TMP = 0.D0
        WRITE(IFOUT1,'(2I5,ES12.3)') IX, IY, TMP
        TMP = T(IX,IY)
        WRITE(IFOUT2,'(2I5,ES12.3)') IX, IY, TMP
     END DO
     WRITE(IFOUT1,*)
     WRITE(IFOUT2,*)
  END DO

  CLOSE ( IFOUT1 )
  CLOSE ( IFOUT2 )

  WRITE(*,1000) JINCRE, NINCRE, JGRAPH, DELTAT*DBLE(JINCRE)
1000 FORMAT (I7,"/",I7,2X,"No.",I3,F6.2,"s")

  RETURN

END SUBROUTINE OUTPUT_VALUE

!=============================================================================
!=============================================================================
