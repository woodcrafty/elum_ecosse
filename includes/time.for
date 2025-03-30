      COMMON /TIME/ N_STEPS, CONVER_F, N_REMAIN,
     &              T_GERM

	INTEGER N_STEPS, N_REMAIN,T_GERM
	
	REAL CONVER_F
C
C TIMESTEP 1 = Daily
C          2 = Weekly
C          3 = Monthly
C
C N_STEPS	: Number of timesteps in one year 
C
C CONVER_F :  Conversion factor
C                1 / 7 for Daily 
C                1     for Weekly
C                52 / 12 or 365 / (7 * 12) for monthly