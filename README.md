# Community elementary conversion modes — uranium-reducing community

Computes elementary conversion modes (ECMs) for each species at a fixed growth rate, then combines them into elementary community modes that describe feasible community compositions.

**Step 1: Build ecmtool inputs** (`make_ecmtool_inputs.ipynb`)
* For each species (Geobacter, Rhodoferax) and each exchange pattern, restrict exchanges to the pattern's allowed uptakes/secretions and calculate the max growth rate.
* For a set of growth-rate fractions (0.3–0.99 of the slower-growing species' max), fix growth rate and carbon uptake, and convert the model to a bound-free form.
* Export one SBML file per species/pattern/growth-rate combination - these are the direct ecmtool inputs.

**Step 2: Compute ECMs** (`run_models_ecmtool.sh`)
Runs ecmtool on each SBML file to compute the full conversion cone.

**Step 3: Filter ECMs and compute community modes** (`analyze_ecms.ipynb`)
* Filter ECMs down to the ones matching the pattern's cross-feeding scenario.
* Combine both species' filtered ECMs into one matrix and run efmtool on it. The resulting EFMs are the elementary community modes - combinations of Geobacter and Rhodoferax ECMs - from which community composition (relative contribution of each species) is calculated.

**Step 4: Validate against PyCoMo** (`analyse_deltas.ipynb`)
Compares the composition ranges from Step 3 against independently computed PyCoMo composition ranges, and checks per-species biomass/product yields.

**Toy example** (`toy_example.ipynb`)
Shows the calculation on a small toy example where each species has just one or two ECMs.