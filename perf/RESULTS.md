# METIS Bit-Identical Optimization Campaign — Results

Branch: `perf-bit-identical`. Driven by `PERF-REVIEW.md`. Every change here is **bit-identical**
to baseline output (verified: all 10 configs `cmp`-equal to `perf/ref/` for `seed=12345`).

## Method
- Build: `make config cc=gcc-15 && make` (gcc-15, `-O3 -march=native`, single-threaded → deterministic).
- Gate: `perf/harness.sh verify perf/ref` — must show `pass=10 fail=0` (byte-identical output).
- Timing: `perf/harness.sh bench <label> 3` — serial, min-of-3, phase timers via `-dbglvl=2`.
- Machine: Apple Silicon (10 cores), Darwin 25.5.0. Runs serial, machine otherwise idle.
- Raw rows accumulate in `perf/RESULTS.tsv`.

## Baseline (serial, min of 3) — seconds

| config | metis | coarsen | contract | initpart | refine |
|---|---|---|---|---|---|
| kway_cit_10 | 8.108 | 2.476 | 1.585 | 3.055 | 1.796 |
| kway_cit_50 | 11.202 | 2.508 | 1.617 | 2.672 | 5.008 |
| kway_cit_100 | 13.148 | 2.370 | 1.532 | 2.593 | 7.007 |
| kway_mdual_10 | 0.064 | 0.039 | 0.023 | 0.016 | 0.002 |
| kway_mdual_50 | 0.080 | 0.038 | 0.022 | 0.024 | 0.009 |
| kway_mdual_100 | 0.116 | 0.037 | 0.022 | 0.051 | 0.017 |
| rb_mdual_10 | 0.195 | 0.148 | 0.088 | 0.000 | 0.020 |
| rb_mdual_50 | 0.277 | 0.197 | 0.114 | 0.002 | 0.038 |
| rb_mdual_100 | 0.332 | 0.230 | 0.133 | 0.004 | 0.049 |
| nd_mdual | 2.895 | 0.465 | 0.266 | 0.424 | 1.864 |

Reference edgecuts (seed=12345): cit k10/50/100 = 1673921 / 2741107 / 3213730;
mdual kway 11157/24100/32113; rb 10313/22831/30689; nd opcount 5.964e+10.

## Optimizations (each verified bit-identical before measuring)

| ID | Description | Status | Key effect |
|---|---|---|---|
| T2.2 | Delete dead cmap stores in matching | — | coarsening/matching |
| T2.6 | SHEM condition reorder | — | matching |
| T2.8 | MMDOrder dead restore pass | — | ndmetis |
| T2.1 | Contraction hash-probe (htkeys) | — | contraction |
| T2.3 | Fuse rename into split extraction | — | rb/nd splitting |
| T2.5 | cnbrpool presizing | — | kway projection |
| T2.7 | sqrt table (gain priority) | — | kway refine |
| T2.9 | CompressGraph early abort | — | ndmetis front-end |
| T2.4 | Reuse pqueues/arrays across levels | — | all refine |
| T2.10 | SoA cnbr_t pool | — | kway refine body |
| T3.1 | LTO / static-inline pqueue | — | refine call overhead |

(Detailed per-opt notes appended below as they land.)
