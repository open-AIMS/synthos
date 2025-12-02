# Generate Settings for Synthetic Reef Simulations

Generate Simulation Settings for Synthetic Reef Landscapes

## Usage

``` r
generateSettings(nreefs, nsites, nyears)
```

## Arguments

- nreefs:

  Integer. Number of reef locations to simulate.

- nyears:

  Integer. Number of years to simulate.

## Value

Invisibly returns a list containing:

- `config_sp` – spatio-temporal model configuration

- `config_lrge` – large-scale sampling configuration

- `config_fine` – fine-scale sampling configuration

- `config_pt` – point-based sampling configuration

Each configuration list is also assigned to the global environment.

## Details

Creates and saves all configuration lists required to simulate synthetic
reef landscapes, including spatial, large-scale, fine-scale, and
point-based sampling parameters. Also assigns each configuration to the
global environment for downstream modelling functions.

The function prepares:

- Spatio-temporal parameters for the synthetic field and disturbances

- Large-scale sampling parameters (reef-level)

- Fine-scale sampling parameters (site and transect-level)

- Point-based sampling parameters (photo-quadrats)

All configurations are saved to disk as `lists_of_parameters.RData`.

## Author

Murray
