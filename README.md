# wildfire-influence
Repo supporting Burke et al 2023 "Wildfire influence on recent US pollution trends".

Results from the paper are in the `figures/clean` and `tables/clean` folders. Code to replicate results are in the `scripts` folder. Data are in [Dropbox](https://www.dropbox.com/sh/3zz7ri3uzc5uf6t/AAAcwLegWlEkA31EkDXuEPZna?dl=0).

## How to replicate results
1. Download this repository.
2. Download the [Dropbox](https://www.dropbox.com/sh/3zz7ri3uzc5uf6t/AAAcwLegWlEkA31EkDXuEPZna?dl=0) folder. Place files downloaded from Dropbox in the same folder as the downloaded GitHub repository.
3. Change settings in `scripts/setup/00_03_load_settings.R`:
    1. Set `path_dropbox` to the location of the data downloaded from Dropbox.
    2. Set `path_github` to the location of this downloaded repository's root.
4. Install packages listed in `scripts/setup/00_01_load_packages.R`.
5. Set working directory to this downloaded repository's root.
6. Run scripts in `scripts/main`.
