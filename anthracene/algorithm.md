1. Included:
  - anthracene.top, the same top file for all simulations.
  - dg.gro, which can be modified by each run by replacing `init-lambda-state` with the lambda state to be run at.
  - lambda0X0.gro to initialize simulations at states 0, 10, 20, . . . , 100.

  - A sample analysis script `dg.py` which takes in all of the .dhdl files, and runs MBAR.

2. Algorithm would be after each simulation comes in, run the analysis
(`basic_analysis.py`) on the `dhdl.xvg` simulation outputs that have
come in, and decide which lambda the next simulation (or simulations
should be run at) should be run at.

3. Initial algorithm: Identify where the largest estimated uncertainty
is between two neighboring states, and run simulations at state that
is in between the ones with the largest estimated uncertainty is,
i.e. if the largest error is between states 60 and 70, run a
simulation at state 65.  If there is no state between them, run at one
of the two states (whichever has highest uncertainty in the next state
over.  We will probably have to iterate.

This information is stored in `ti.d_delta_f_` and `mbar.d_delta_f_`,
in the off-diagonal elements.  Note that ti.d_delta_f_ only shows data
at states that HAVE samples, whereas `mbar.d_delta_f_` shows data at
all states, regardless of whether they have samples. So if you use
`mbar.d_delta_f_`, you don't have to bisect, you could just run a
simulation at the states that have the highest uncertainty, which is
easier.

Note that if you are running multiple simulations, the states with the
lowest uncertainty might be next to each other, so you might want to
include some criteria that the next simulation is at least delta
lambda = 0.1 from currently running simulations or something, as the
two neigboring simulations would be duplicative.

To edit the .mdp for the next run, you just need to fill in the value
for `init-lambda-state`. 

To pick which .gro to use, if there does not exist a simulation with
the same lambda value as the last .gro. pick whichever one is closest
in the LOWER direction (i.e. for 0.69, pick 0.6).  This is because the
particle is disappearing in the increasing lambda direction, and it's
easier to start with a simulation that is at a higher level of
insertion.

The length of the current simulation seems reasonable.  The analysis
completed without error for 1 dhdl.xvg this length. On 4 cores, it
took about 10 min for me. per itartion

A sample `sample_output.txt` of the analysis is included which ran from the .dhdl.xvg files included.

Might be nice to keep the graphs/data over time to see how they evolve
over time. We will come up with better ways to show this graphically.


We will come up with better versions to show this eventually.

`choosing_lambda.ipynb` is a notebook that describes the algorithm to select new states
based on the variance measured so far.

