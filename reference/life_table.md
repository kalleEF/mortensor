# Life tables from death rates on any age axis

Period life table for single or grouped ages with an open last interval.
Interval widths come from the age axis, so abridged and single-year
schedules are handled by the same code.

## Usage

``` r
life_table(
  mx,
  age = NULL,
  ax = NULL,
  open = c("reciprocal", "gompertz", "kannisto"),
  k = 10,
  radix = 1e+05
)
```

## Arguments

- mx:

  Death rates, one per interval of `age`.

- age:

  An [mt_axis](mt_axis.md) (default: single ages from 0 with an open
  last interval).

- ax:

  Optional average years lived in the interval by those who die.

- open:

  `"reciprocal"`, `"gompertz"` or `"kannisto"`.

- k:

  Number of ages used for the Gompertz fit.

- radix:

  Starting number of survivors.

## Value

A `data.table` with `age_lo`, `age_hi`, `n`, `mx`, `ax`, `qx`, `lx`,
`dx`, `Lx`, `Tx`, `ex`.

## Details

- `q_x = n m / (1 + (n - a) m)`; where this exceeds 1 (rates above about
  `2 / n`, far beyond observed values), the constant-hazard value
  `1 - exp(-n m)` is used, so that survivors stay positive;

- `a_x` defaults to `n / 2` for every closed interval (supply `ax` to
  override, for example for age 0, where `n / 2` overstates the average
  time lived by infants who die; the effect on e0 at low infant
  mortality is below 0.01 years);

- open interval: `"reciprocal"` uses `L = l / m`; `"gompertz"` fits a
  Gompertz line to the log rates of the last `k` closed single-year ages
  and integrates it from the start of the open interval to age 130;
  `"kannisto"` does the same with the Kannisto model, a line for the
  logit of the rates, whose hazard levels off at old ages (Thatcher,
  Kannisto and Vaupel 1998);

- a rate of 0 in the open interval (no deaths, as in small populations)
  would give an infinite life expectancy: the last positive rate of the
  closed intervals is used for the open interval instead (column `mx`
  keeps the given rate). Only a schedule without any positive rate gives
  an infinite value.

## References

Thatcher AR, Kannisto V, Vaupel JW (1998). *The Force of Mortality at
Ages 80 to 120*. Odense University Press.

## Examples

``` r
m <- 3e-4 + 2e-5 * exp(0.1 * (30:100))
lt <- life_table(m, single_ages(30, 100))
lt[age_lo %in% c(30, 65, 100)]
#>    age_lo age_hi     n           mx       ax           qx         lx         dx
#>     <num>  <num> <num>        <num>    <num>        <num>      <num>      <num>
#> 1:     30     31     1 0.0007017107 0.500000 0.0007014646 100000.000   70.14646
#> 2:     65     66     1 0.0136028327 0.500000 0.0135109391  87531.735 1182.63594
#> 3:    100    Inf   Inf 0.4408293159 2.268452 1.0000000000   1460.141 1460.14118
#>          Lx         Tx        ex
#>       <num>      <num>     <num>
#> 1: 99964.93 4973452.39 49.734524
#> 2: 86940.42 1596839.74 18.242981
#> 3:  3312.26    3312.26  2.268452
```
