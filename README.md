# CVP Infrastructure

## Examples & Tracks

See Simulator options:

`./cvp`

Running the simulator on `trace.gz`:

`./cvp trace.gz`

Running with value prediction (`-v`) and predict all instructions (first track, `-t 0`):

`./cvp -v -t 0 trace.gz`

Running with value prediction (`-v`) and predict only loads (second track, `-t 1`):

`./cvp -v -t 1 trace.gz`

Running with value prediction (`-v`) and predict only loads but with hit/miss information (third track, `-t 2`):

`./cvp -v -t 2 trace.gz`

Baseline (no arguments) is equivalent to (prefetcher (`-P`), 512-window (`-w 512`), 8 LDST (`-M 8`), 16 ALU (`-A 16`), 16-wide fetch (`16`) with 16 branches max per cycle (`16`), stop fetch at cond taken (`1`), stop fetch at indirect (`1`), model ICache (`1`)):

`./cvp -P -w 512 -M 8 -A 16 -F 16,16,1,1,1`

## Run bash file

'bash run.sh'

## Notes

Run `make clean && make` to ensure your changes are taken into account.

## Value Predictor Interface

See [cvp.h](./cvp.h) header.

## Getting Traces

135 30M Traces @ [TAMU ](http://hpca23.cse.tamu.edu/CVP-1/public_traces/)

2013 100M Traces @ [TAMU ](http://hpca23.cse.tamu.edu/CVP-1/secret_traces/)


## Sample Output

```VP_ENABLE = 1
VP_PERFECT = 0
VP_TRACK = LoadsOnly
WINDOW_SIZE = 256
FETCH_WIDTH = 16
FETCH_NUM_BRANCH = 0
FETCH_STOP_AT_INDIRECT = 0
FETCH_STOP_AT_TAKEN = 0
FETCH_MODEL_ICACHE = 0
PERFECT_BRANCH_PRED = 0
PERFECT_INDIRECT_PRED = 0
PIPELINE_FILL_LATENCY = 5
NUM_LDST_LANES = 0 (unbounded)
NUM_ALU_LANES = 0 (unbounded)
MEMORY HIERARCHY CONFIGURATION---------------------
STRIDE Prefetcher = 1
PERFECT_CACHE = 0
WRITE_ALLOCATE = 1
Within-pipeline factors:
        AGEN latency = 1 cycle
        Store Queue (SQ): SQ size = window size, oracle memory disambiguation, store-load forwarding = 1 cycle after store's or load's agen.
        * Note: A store searches the L1$ at commit. The store is released
        * from the SQ and window, whether it hits or misses. Store misses
        * are buffered until the block is allocated and the store is
        * performed in the L1$. While buffered, conflicting loads get
        * the store's data as they would from the SQ.
L1$: 32 KB, 4-way set-assoc., 64B block size, 2-cycle search latency
L2$: 1 MB, 8-way set-assoc., 64B block size, 12-cycle search latency
L3$: 8 MB, 16-way set-assoc., 128B block size, 60-cycle search latency
Main Memory: 150-cycle fixed search time
STORE QUEUE MEASUREMENTS---------------------------
Number of loads: 59570
Number of loads that miss in SQ: 49200 (82.59%)
Number of PFs issued to the memory system 4487
MEMORY HIERARCHY MEASUREMENTS----------------------
L1$:
        accesses   = 280836
        misses     = 200283
        miss ratio = 71.32%
        pf accesses   = 4487
        pf misses     = 249
        pf miss ratio = 5.55%
L2$:
        accesses   = 200283
        misses     = 11547
        miss ratio = 5.77%
        pf accesses   = 249
        pf misses     = 35
        pf miss ratio = 14.06%
L3$:
        accesses   = 11547
        misses     = 5857
        miss ratio = 50.72%
        pf accesses   = 35
        pf misses     = 30
        pf miss ratio = 85.71%
BRANCH PREDICTION MEASUREMENTS---------------------
Type                      n          m     mr  mpki
All                 1020870       1487  0.15%  1.46
Branch               227459       1172  0.52%  1.15
Jump: Direct           6301          0  0.00%  0.00
Jump: Indirect         4647        315  6.78%  0.31
Jump: Return              0          0  -nan%  0.00
Not control          782463          0  0.00%  0.00
ILP LIMIT STUDY------------------------------------
instructions = 1020870
cycles       = 326967
IPC          = 3.122
Prefetcher------------------------------------------
Num Trainings :59570
Num Prefetches generated :4487
Num Prefetches issued :4487
Num Prefetches filtered by PF queue :0
Num untimely prefetches dropped from PF queue :0
Num prefetches not issued LDST contention :0
Num prefetches not issued stride 0 :10391
CVP STUDY------------------------------------------
prediction-eligible instructions = 59570
correct predictions              = 20575 (34.54%)
incorrect predictions            = 3053 (5.13%)
NP_correct predictions              = 2396 (4.02%)
NP_incorrect predictions              = 33546 (56.31%)
----------------
Number of aluInstClass : 0
Number of condBranchInstClass : 0
Number of uncondDirectBranchInstClass : 0
Number of uncondIndirectBranchInstClass : 0
Number of fpInstClass : 0
Number of slowAluInstClass : 0
Number of undefInstClass : 0
Other------------------------------------------
Misclassified LDP to LD + BU: 4
Misclassified STP to STR (once, total): 96
Misclassified STP to STR (multiple times, total): 0
Misclassified STR to STP (once times, total): 0
Misclassify load vector lane: 0
 Read 1003520 instrs
```
