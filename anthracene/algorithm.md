1. Included:
  - anthracene.top, the same top file for all simulations.
  - dg.gro, which can be modified by each run by replacing `init-lambda-state` with the lambda state to be run at.
  - lambda0X0.gro to initialize simulations at states 0, 10, 20, . . . , 100.

  - A sample analysis script `dg.py` which takes in all of the .dhdl files, and runs MBAR.

2. Algorithm would be after each simulation comes in, run the
analysis on the simulations that have come in, and decide which lambda
the next simulations should be run at.

3. Initial algorithm: Identify where the largest estimated uncertainty
is, and run simulations at state that is in between the ones with the
largest estimated uncertainty is, i.e. if the largest error is between
states 60 and 70, run simulations at state 65.  If there is no state
between them, run at one of the two states (whichever has highest
uncertainty in the neighboringing state).