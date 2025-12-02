# Large‐scale fixed sampling design

Large‐scale fixed sampling design

## Usage

``` r
sampling_design_large_scale_fixed(data_reefs_pts_sf, config_lrge)
```

## Arguments

- data_reefs_pts_sf:

  An `sf` object representing the full field. It must contain:

  - `Year`: numeric year

  - `Reef`: unique reef ID

  - `HCC`, `SC`, `MA`: cover values (logit scale)

  - `geometry`: spatial geometry

- config_lrge:

  A list with:

  - `n_locs`: number of reefs (locations) to select

  - `n_sites`: number of sites per selected reef

  - `seed`: random seed

## Value

An `sf` object representing the large-scale sampling design. Cover
values remain on the logit scale.

## Details

Selects a set number of reef locations and then randomly selects a fixed
number of sites within each selected location.

## Author

Murray

## Examples

``` r
config_lrge <- list(n_locs = 25, n_sites = 2, seed = 123)
benthos_fixed_locs_sf <- sampling_design_large_scale_fixed(
  benthos_reefs_pts, config_lrge
)
#> Error: object 'benthos_reefs_pts' not found
```
