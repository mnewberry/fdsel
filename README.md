# fdsel

`fdsel` is a command-line interface utility for inferring frequency-dependent
selection from timeseries data on abundance of types, written in OCaml.

## Building from source

`fdsel` requires OCaml and libraries for PCRE, Batteries and the GNU Scientific
Library (GSL). This requires both OCaml and system libraries. PCRE and GSL
system libraries can typically be found through a distribution package manager
such as apt, pacman, or homebrew. OCaml and the OCaml libraries are easiest to
install using OPAM. Once OPAM is installed and configured, the libraries can be
installed with:

```opam install batteries pcre gsl
```

Once the dependencies are installed, building is easiest using findlib and
ocamlbuild, which are also included as part of OPAM:

```ocamlbuild -pkgs gsl,pcre,batteries fdsel.native
```

This command is included for convenience in a file called make.sh. Hence `sh
make.sh` is also an acceptable way to build the software.

The build will create the file `fdsel.native` which is the self-contained fdsel
executable. You may rename it `fdsel` and put it in your executable path as you
see fit, or invoke the software as `./fdsel.native`.

## `fdsel` usage overview

The reference usage of `fdsel` can be found by running `fdsel --help`. `fdsel`
has three modes, which can be invoked as `fdsel infer` (for inference from
timeseries), `fdsel simulate` (for simulating timeseries), and `fdsel
timeseries` (for examining properties of timeseries without inferring or
simulating anything). Each mode has a list of options and flags, with some
features overlapping between modes.

### File formats 

`fdsel` inputs and outputs various files in tab-separated values (TSV) formats.
The output formats are designed to be self-explanatory, whereas the few input
formats (timeseries, breaks and params) are also output formats. The input
timeseries for `fdsel infer` and `fdsel timeseries` should have a header with
columns `gen`, `type`, and `count` in any order, which is the same format as
the output `fdsel simulate -o`. Only the bin boundaries files (`-B` arguments)
have no header. The parameters file (`-o` argument to `fdsel infer` and `-i`
argument to `fdsel simulate`) is complex due to the many different output
parameters, but it is also a TSV file with headers for the parameter `param`,
its index `ind` (if applicable), the value `val`, its variance `var` (if
applicable), and 95% confidence interval bounds `lci` and `uci`.

### Example

A simple session illustrates the basic features. First, we simulate a neutral
timeseries with 80 generations and a population size of 5000, per capita
mutation rate 0.001 and 5000 generations of burn-in initialized from a
monomorphic population, storing the output in `timeseries.tsv`:

```fdsel simulate -K -g 80 -n 5000 -m 0.001 -b 10000 -o timeseries.tsv
```

Now `timeseries.tsv` looks like this:

```$ head timeseries.tsv 
gen  type   count
0    24944  620
0    22545  434
0    24074  322
0    25597  312
0    25580  298
0    24508  295
0    20888  265
0    22892  263
0    21411  226
```

We see that the initial type (1) is gone, and that the largest counts are
comparable in magnitude, so we may assume that 10,000 burn-in generations was
enough to reach equilibrium.  We can infer frequency-dependent selection from
this timeseries using 5 log-evenly spaced bins, storing the result in
`params.tsv` and recording the bin boundaries in `breaks.nlsv`:

```fdsel infer -i timeseries.tsv -l 5 -o params.tsv -b breaks.nlsv
```

The resulting parameters file is consistent with neutrality, in the sense that
the confidence intervals (`lci` and `uci`) always include s=0. The inferences
for mutation rate and population size also include the true value. (Here tabs
have been replaced by spaces and some decimal places removed for clarity of
presentation.)

```$ cat params.tsv
param     ind  val        var       lci  uci
mu        NA   0.0010225  2.55e-09  0.0009234  0.0011215
ne        NA   4873.9416  8560.471  4692.5969  5055.2863
maxf      NA   2.152e-01  NA        2.152e-01  2.152e-01
dminf     NA   2.000e-04  NA        2.000e-04  2.000e-04
srep      NA   0.0026665  NA        0.0026665  0.0026665
s1        0    0.0129530  0.000181  -0.013467  0.0393736
s2        1    -0.007467  7.347-05  -0.024268  0.0093332
s3        2    -0.008886  3.286-05  -0.020122  0.0023494
s4        3    -0.001443  2.219-05  -0.010677  0.0077909
s5        4    0.0048440  1.897-05  -0.003693  0.0133818
ll        NA   -14648.50  NA        NA         NA
s0ll      NA   -14652.90  NA        NA         NA
s0ne      NA   4869.5538  8545.064  4688.2090  5050.8985
s0lrpval  NA   6.615e-02  NA        NA         NA
```

This file indicates the inferred mutation rate `mu` and parameter values
`s1`...`s5` of the different bins, the replacement fitness `srep`, the
log-likelihood of the data `ll`, the log-likelihood assuming neutrality `s0ll`,
and the p-value in the likelihood ratio test between the two `s0lrpval`. Here,
the p-value 0.07 does not quite reject neutrality. `maxf` and `dminf` refer to
the maximum and minimum frequencies present in the input.

These parameters deviate somewhat from neutrality just due to statistical
error. But we can use those slight deviations to produce a frequency-dependent
simulation, by replacing the `-K` (neutrality) and `-m` (mutation) arguments
with the parameters input file. If we simulate for a long time (g = 10,000), we
should be able to recover parameters very close to the inputs. This time, since
we know the bin boundaries used in the simulation, we can specify that the
inference use those exactly.

```fdsel simulate -i params.tsv -B breaks.nlsv \
  -g 10000 -n 5000 -b 10000 -o timeseries-selected.tsv
fdsel infer -i timeseries-selected.tsv -B breaks.nlsv -o params-selected.tsv
```

After a few seconds, we get:

```$ cat params-selected.tsv
param     ind  val        var        lci        uci
mu        NA   0.0010222  2.042e-11  0.0010133  0.0010310
ne        NA   4981.2651  109.43420  4960.7614  5001.7689
maxf      NA   9.694e-01  NA         9.694e-01  9.694e-01
dminf     NA   2.000e-04  NA         2.000e-04  2.000e-04
srep      NA   0.0055515  NA         0.0055515  0.0055515
s1        0    0.0105344  1.587e-06  0.0080648  0.0130040
s2        1    -0.006567  7.226e-07  -0.008233  -0.004901
s3        2    -0.008844  4.607e-07  -0.010175  -0.007514
s4        3    -0.000683  6.128e-07  -0.002218  0.0008505
s5        4    0.0055618  1.923e-07  0.0047022  0.0064214
ll        NA   -952663.9  NA         NA         NA
s0ll      NA   -952973.3  NA         NA         NA
s0ne      NA   4982.5432  109.49037  4962.0395  5003.0469
s0lrpval  NA   1.30e-132  NA         NA         NA
```

Although it's a little harder to check this time, the CIs include the true
value, and the inferred parameters closely resemble the inputs. This time,
however, the likelihood ratio test overwhelmingly rejects neutrality.

## Advanced uses

More advanced uses of `fdsel infer` include other binning schemes (`-k`, `-l`,
`-q` and `-B`, `-f`), strategies for dealing with censorship or types at unknown
frequency (`-w`, `-c`, `-f` and `-C`), or communicating assumptions about datasets
that do not distinguish missing data from types at count 0 (`-u` and `-U`).
Similarly, `fdsel simulate` also includes options for producing incoming
mutants in clusters (`-l`), specifying variable population sizes, simulating
novelty bias models, specifying the initial population (`-T`) or continuing
simulations from a previous timeseries (`-F`), or interspersing each timestep
with neutral timesteps in order to simulate reduced effective population size
relative to the census population. `fdsel infer` and `fdsel timeseries` also
support several other outputs besides the inferred parameters, including the
population size history (`-p`), frequency distribution (`-d`) or rank abundance
distribution a.k.a. Zipf's law distribution (`-Z`).
