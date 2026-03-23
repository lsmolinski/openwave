# PHASE 3: Forces — 3D Validation

Port the validated 1D equations and force computation to M3 Taichi 3D engine. Verify that 3D results match 1D on-axis results, then test off-axis configurations for spherical symmetry. Validate force & motion integration — particles should move correctly under computed forces.

Ports validated 1D results from Phase 1c (non-linear) and/or Phase 1d (vector) to 3D engines. Phases 1c and 1d are deeply connected and may converge here.

**Non-linear (from Phase 1c → M3 3D):**

- **Yee & Hauger**: discrete wavelength shells r_n = 2(K-n)λ — non-uniform node spacing breaks sinc periodicity (WKB/eikonal phase integral)
- **Smoliński**: r⁵ energy scaling inside soliton's Energy Domain — defines how λ(r) varies near the wave center core
- **Ψ³ cubic non-linearity**: NLS soliton stabilizer, modifies spatial function from pure sinc
- **Variable ρ(x)**: from granule velocity — ∇ρ contributes to ∇E with different spatial structure than ∇A

**Vector (from Phase 1d → M4 3D):**

- **Divergence/curl/flux** force computation in full 3D vector field
- **Elliptical rotation handedness** as charge-sign indicator (6-phasor model)
- **Toroidal flow geometry**: non-linear + vector converge in the Energy Domain

See PHASE 1 docs for full analysis.
