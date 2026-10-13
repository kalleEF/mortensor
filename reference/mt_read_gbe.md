# Read tables exported from GBE-Bund

Parses the CSV export of tables of the Federal Health Monitoring
(Gesundheitsberichterstattung des Bundes, <https://www.gbe-bund.de>): a
few title lines, two header lines (for example years and causes of
death) and data rows whose first columns hold the row variables (for
example region, age and sex), each written only where it changes. The
file is Latin-1 encoded with `;` as separator.

## Usage

``` r
mt_read_gbe(
  file,
  row_names = c("region", "age", "sex"),
  col_names = c("year", "cause")
)
```

## Arguments

- file:

  Path of the exported CSV file.

- row_names:

  Names of the row variables, in their order.

- col_names:

  Names of the two column variables (outer, inner).

## Value

A `data.table` with the row and column variables, `value` and `symbol`.

## Details

Symbols of official statistics are kept in column `symbol`: `"-"`
(nothing, value 0), `"."` (unknown or confidential, value `NA`), `"X"`
(not applicable, value `NA`), `"..."` (not yet available, value `NA`).
