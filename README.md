# dkmfr
Synthesis of reversible circuits.

## Installing required modules:
`pip install sympy`
`pip install argparse` 

or simply use:

`pip install -r requirements.txt`


## Usage:
`dkmfr.py [-h] [--bi BI] [--functions FUNCTIONS [FUNCTIONS ...]] [--samples SAMPLES] [--nbits NBITS]`

### Example:
`$ DKMFR.py --bi 2 --functions hwb4 4b15g_2`

### Optional arguments:
  `-h`, `--help` show this help message and exit.
  
  `--bi BI` choose 1 for forward only, and choose 2 for forward/reverse.
  
  `--functions FUNCTIONS [FUNCTIONS ...]` Choose benchmark functions from `4b15g_1`, `4b15g_2`, `4b15g_3`, `4b15g_4`, `4b15g_5`, `hwb4`, `4_49`, `toffoli_double_2`, `nth_prime4_inc`, `ham3_complete_47(28)`, `miller_complete_5`, `toffoli_1`, `ex-1_82`, `3_17`, `nth_prime3_inc`, `RCFK`, `ex1Miller`, `ex2Miller`, `ex3Miller`, `ex4Miller`, `ex5Miller`, `ex6Miller`, `ex7Miller`, `aj-e11_complete_74(81)`, `mod5mils_complete_26(18)`, `nth_prime5_inc`, `hwb5_13`, `mod5adder`, `hwb6`, `nth_prime6_inc`, `graycode6_complete_19`, `nth_prime7_inc`, `ham7`, `hwb7_15`, `random`. You can choose multiple functions separating them by spaces.
  
  `--samples SAMPLES` how many random permutations will be generated.
  
  `--nbits NBITS` how many bits in random permutations.

