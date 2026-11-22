# Semi-miscible CO₂–brine flow in stratified formations

MATLAB code for the analytical solution and upscaling workflow in:

> Youssef, A.A. and El-Husseiny, A. *Semi-Miscible CO₂–Brine Flow in Stratified Formations: An Analytical Framework for Upscaled Relative Permeability and Dissolved-Flux Correction.* Submitted to *Advances in Water Resources* (under review). DOI to be added on publication.

The code computes viscous-dominated displacement of CO₂-free brine by supercritical CO₂ in **non-communicating layered media**, accounting for CO₂ dissolution in the aqueous phase. It combines:

- the semi-miscible fractional-flow theory of Noh et al. (2007), where the Welge tangent starts at (1 + D, 1 + D), with D = χρ_aq / (ρ_sc − χρ_aq); and
- the exact layered (non-communicating) solution of Youssef and Matthäi (2022).

It returns frontal advance in each layer, outlet fractional flow, pressure drop, free and dissolved CO₂ volumes, and the dissolved-flux correction factor α. A second script history-matches an equivalent homogeneous model to obtain upscaled relative permeability.

## Requirements

| Script | Requirement |
|---|---|
| `MultiLayer.m`, `soluble_flux.m` | MATLAB R2018b or later (base MATLAB only) |
| `kr_inverse.m`, `kr_inverse_BC.m` | as above, plus the **Optimization Toolbox** (`fmincon`) |

No other toolboxes or external files are needed. All functions are in the repository root.

## Quick start: reproduce the Otway CRC-3 case

1. Clone or download the repository and open the folder in MATLAB (make it the current folder):
   ```
   git clone https://github.com/Abdallahyoussef/semi-miscible-layered-upscaling.git
   ```
   or use **Code → Download ZIP** on GitHub.
2. Run the analytical solution:
   ```matlab
   MultiLayer
   ```
   The inputs in this file are the Otway CRC-3 values used in the manuscript (5 grouped flowing layers, Table 2; thermodynamic conditions, Table 3; Brooks–Corey relative permeability with λ = 2, s_wr = 0.2, k_rnw end-point = 0.9). The script prints the front saturation of each layer, plots the storage volumes, outlet fractional flow and α, and saves all results to `Otway_miscible_results.mat`.
3. Estimate the upscaled relative permeability (Chierici model) by history matching pressure drop and outlet fractional flow (Eq. 7):
   ```matlab
   kr_inverse
   ```
   This reads `Otway_miscible_results.mat`. `kr_inverse_BC.m` does the same with a single-parameter Brooks–Corey model.
4. Compare dissolved-to-free CO₂ flux with and without the α correction:
   ```matlab
   soluble_flux
   ```

Expected checks for the Otway case:
- Welge tangent starting point 1.0732; flood-front aqueous saturation 0.71 and fractional flow 0.17 (Fig. 4);
- breakthrough of the fastest layer at about 0.24 PVI, with free and dissolved CO₂ in a ratio of 8.12 : 1 and total storage below 25% of the pore volume (Fig. 5a);
- stepwise outlet fractional flow and α as later layers break through (Fig. 5c).

`MultiLayer.m` uses `parfor` loops. Without the Parallel Computing Toolbox these run as ordinary `for` loops, so no extra toolbox is needed.

## Running your own case

All inputs are at the top of `MultiLayer.m`, in three blocks:

| Block | Inputs |
|---|---|
| **Patch 1** – geometry and single-phase properties | `layers(i).thickness`, `layers(i).porosity`, `layers(i).perm_abs`. Add or remove layers freely. |
| **Patch 2** – two-phase properties | `mu_w`, `mu_nw`; per-layer `models` (`'BC'`, `'VG'`, `'LET'`, `'Chierici'`), `params` (model parameters, one cell per layer), `Swr`, `Snwr`, `Krw_end`, `Krnw_end`, `SwInitial` (one entry per layer). See `relativePermeability.m` for the parameters of each model. |
| **Patch 3** – process and controls | `Process` (`'drainage'` or `'imbibition'`); `solubility_into_displaced_phase` = χρ_aq (set to 0 for immiscible flow); `solubility_into_displacing_phase` = ρ_sc; `Resolution`; `case_name`; `results_file`. |

The script checks that every per-layer input has one entry per layer and stops with a clear message if not.

**Units.** Only ratios enter the solution, so any consistent units can be used. The Otway inputs use thickness in cm, permeability in mD, viscosity in cP and mass/volume in kg/m³.

## Assumptions and limitations

As stated in Section 3.1 of the manuscript:

- non-communicating layers (no cross-flow), connected only at inlet and outlet;
- viscous-dominated flow: capillary and gravity forces are neglected;
- incompressible fluids and rock; isothermal flow;
- instantaneous local equilibrium between contacted CO₂ and brine, with a constant CO₂ solubility evaluated at a representative pressure;
- H₂O vaporization into the CO₂-rich phase is neglected;
- constant total injection rate at the inlet, fixed pressure at the outlet.

The miscible (semi-miscible) option is implemented for **drainage** (CO₂ displacing brine, initial aqueous saturation of 1). Immiscible flow can be run as drainage or imbibition.

## Files

| File | Purpose |
|---|---|
| `MultiLayer.m` | **Main script.** Inputs, workflow, plots, saves results |
| `relativePermeability.m` | Relative permeability, mobility and fractional-flow functions (Brooks–Corey, van Genuchten, LET, Chierici) |
| `findShockFront.m` | Welge tangent / shock-front saturation, including the semi-miscible tangent origin |
| `computeGammaParameter3.m` | Integral parameter γ (effective total mobility) of the layered solution |
| `computeDoubleIntegralTerm.m` | Double-integral term used after breakthrough of the fastest layer |
| `simulateFirstStage.m` | Solution before breakthrough of the fastest layer |
| `simulateSecondStage3.m` | Solution after breakthrough of the fastest layer, including α |
| `kr_inverse.m` / `kr_inverse_BC.m` | History matching of upscaled relative permeability (Chierici / Brooks–Corey) |
| `soluble_flux.m` | Dissolved-to-free CO₂ flux ratio with and without α |
| `bisection.m` | Bisection root finder used by `findShockFront.m` |
| `isPositiveIntegerArray.m` | Input-validation helper |

## How to cite

If you use this code, please cite the article above (citation details will be updated on publication) and the archived software release (see `CITATION.cff`).

## References

- Chierici, G.L., 1984. Novel Relations for Drainage and Imbibition Relative Permeabilities. Society of Petroleum Engineers Journal 24, 275–276. https://doi.org/10.2118/10165-PA
- Brooks, R.H., Corey, A.T., 1963. Hydraulic properties of porous media. Colorado State University.
- Noh, M., Lake, L.W., Bryant, S.L., Araque-Martinez, A., 2007. Implications of Coupling Fractional Flow and Geochemistry for CO2 Injection in Aquifers. SPE Reservoir Evaluation & Engineering 10, 406–414. https://doi.org/10.2118/89341-PA
- Welge, H.J., 1952. A Simplified Method for Computing Oil Recovery by Gas or Water Drive. Journal of Petroleum Technology 4, 91–98. https://doi.org/10.2118/124-G
- Youssef, A.A., Matthäi, S.K., 2022. Upscaled Dynamic Relative Permeability for Unstable CO2 Flow in Stratified Porous Media. Transport in Porous Media 143, 657–680. https://doi.org/10.1007/s11242-022-01803-6

## License

See `LICENSE`.

## Contact

AbdAllah A. Youssef – abdallah_youssef@ymail.com
