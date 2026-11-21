JutulCO2Lab reimplements MRST-co2lab's vertical-equilibrium (VE) CO2 storage
simulator on Jutul/JutulDarcy, verified against co2lab by numerical
comparison on identical input.  `README.md` describes what it covers, the
installation and the architecture; `STRUCTURES.md` lists every type and
returned object.

## This repository

- **Self-contained.**  The package is shared with colleagues and also
  developed as a submodule inside other projects.  Nothing in it may refer to
  such a project: no roadmap or stage numbers, no pointers to `analysis/` or
  `scripts/` files, no dates or attributions ("researcher, 2026-10-06"), no
  run lines with paths from outside the package.  Case studies built on the
  package live in their own repositories, not here.
- **Manifests are local.**  `Manifest.toml` and `examples/Manifest.toml` are
  ignored: each checkout resolves the registered Jutul, JutulDarcy and
  GeoEnergyIO releases within `Project.toml`'s compat bounds (see
  `README.md`).
- **Commits:** short subject; the body says what changed and why.  Push only
  when asked.

## Layout

- `src/ve` the VE simulator, `src/atlas` the CO2 Storage Atlas (and Johansen,
  Sleipner), `src/trapping` the spill-point trapping analysis,
  `src/inventory` the trapping inventory, `src/oracle.jl` the comparison
  harness; plotting is declared in `src/plotting.jl` and implemented in the
  Makie extension `ext/` (Makie optional: GLMakie, CairoMakie or WGLMakie).
- `examples/<topic>`: readable walkthroughs with their own environment, run
  from the package root as `julia --project=examples examples/<file>.jl`.
  `examples/sloping_aquifer.jl` is the model for layout and comments.
- `test/runtests.jl` and the shared `test/*_checks.jl`.

## Running

- `julia --threads=N --project=. -e 'using Pkg; Pkg.test()'`, or
  `julia --project=. test/runtests.jl` when the test environment cannot
  precompile every extension (e.g. GLMakie without a display).
- Data sets: `CO2LAB_ATLAS_DATA`, else `mrst-data` next to the package's
  directory (or next to the directory containing it).
- MRST reference dumps for the comparisons: `JUTULCO2LAB_REFERENCE_DATA`,
  else `analysis/reference-data` next to the package's directory.  Tests
  needing a missing dump or data set skip themselves.
- Reference dumps are generated under Octave, not MATLAB (MRST's
  `getAtlasGrid` gives a different active-cell mask in MATLAB), with MRST on
  branch `gigaccs-working-branch`; each dump records its `mrst_commit`.

## Code conventions

Julia has no `const` for arguments, so mutation is made explicit:

1. A function whose name ends in `!` puts the argument(s) it modifies first.
2. Every function that modifies anything it was given has a docstring line
   `Modifies: ...` naming exactly what, nested fields included
   (`Modifies: `report.results``).
3. A function without `!` modifies nothing it was given; it returns new
   objects.

Functions that must follow an external signature (Jutul's
`@jutul_secondary` updates, `update_equation_in_entity!`, ...) are exempt
from rule 1 but keep the `Modifies:` line.  Methods added to Jutul/JutulDarcy
functions keep those functions' names.

- **Function length:** one page (about 50 lines, docstring excluded); split
  longer logic into named helpers.
- **Internal helpers** (not exported) start with `_`, so the public interface
  is recognisable at a glance.
- **Struct fields:** every field has a brief aligned inline comment,
  `res_gas::Float64  # residual CO2 saturation`.
- **`STRUCTURES.md`:** update it in the same commit as any change that adds,
  removes, renames or changes the meaning of a type, type alias, constant or
  returned composite object (NamedTuples etc.), keeping dependency order.
- **Return values:** every docstring states what the function returns: type,
  meaning, aliasing ("`mesh` itself when nothing is removed").
- **Line length:** below 90 characters, unless a long string or URL should
  stay searchable.
- **Comments:** a brief comment, preferably one line, before every code line
  that is not obvious, saying what case or step follows and why
  (`_co2lab_multimodel` in `src/ve/co2lab_fluid.jl` shows the level).  File
  header: a few lines on what the file holds and how it fits in.
- **Whitespace and headings** (Julia and MATLAB/Octave): separate logical
  fragments by a blank line; separate larger steps by a heading line of
  exactly 80 characters, the comment marker, a space, and 78 hyphens with the
  heading text, padded by a space on each side, in the middle:

      # ------------------------------ 5. Open boundary ------------------------------

## Docstrings

The model is `ve_domain` in `src/ve/vedomain.jl`.  A reader learns what a
thing is, why it exists and how to use it, not how it is implemented.

- **Order:** signature; one or two sentences on what it does and where it
  fits; `# Arguments` (one bullet per argument or keyword, with units and
  defaults; related arguments may share a bullet); `# Returns` (one bullet
  per field of a returned NamedTuple); the `Modifies:` line; then optional
  topic headings (`# Keywords`, `# Fields`, `# Constructors`, `# Example`,
  or a topic such as `# Faces`) and `# Notes` last (design choices,
  differences from co2lab, pitfalls).
- **Exported functions** always have `# Arguments` (when they take any) and
  `# Returns`, complete enough to use the function without reading its body.
  Small or internal ones may use the compact form: the summary, then
  "Returns ...".  Do not inflate a short docstring into sections.
- **Types:** what the type represents and which model or step uses it;
  abstract types name the family's purpose and concrete subtypes.  Every
  field is described under `# Fields`, one bullet each (that is what `?`
  shows); the inline comments are brief reminders for readers of the code.
- **Writing:** short paragraphs (about four lines at most), else bullets or a
  heading; units on physical quantities; cross-references as
  `` [`name`](@ref) ``.  Keep every fact a user needs (defaults, units,
  conventions, warnings, aliasing, preconditions).  Without `!`, say nothing
  about arguments being only read.  No research history; keep the names of
  co2lab counterparts (`makeVEFluid`, `'endpoint scaling'`) and upstream
  issue numbers.

## Modelling practice

- **Follow Jutul/JutulDarcy conventions** (names, tags, layouts, units);
  translate MRST data at import rather than carrying MRST conventions inward.
- **Differences from co2lab:** where Jutul's way and co2lab's differ in a
  modelling choice, present both to the maintainers before choosing.  Where
  co2lab's solution is inferior, say so and implement the better one,
  keeping co2lab's behaviour behind a switch for comparisons (e.g.
  `co2lab_compatible`, `fault_flux = :co2lab`).
- **Comparisons with co2lab** are per quantity with a declared tolerance
  (`Oracle`), never against a plot, and never with a guessed ordering: pass an
  explicit permutation (Julia index → MRST index).  An empty report fails.
- **Zero entry pressure is the common case:** capillary-fringe code and tests
  must work for `pe = 0`, not rely on `pe > 0`.
- **Poor convergence** (time-step cuts, stalls) points to an error in the
  physics or the code: trace it (which cells, which equations, how the
  iterates move) before changing solver settings.
- **Branch on plain values:** decide branches (`if`, comparisons, upwinding,
  clamps) on `_plain(x)`, not on AD values: under ForwardDiff 1 two values
  that tie are compared by their derivatives, so the branch, and with it the
  Jacobian, would depend on them.
- **Test suite:** small and fast; 3D cases with at most 20 cells vertically.
  High resolutions belong in separate comparison scripts.
- **Plots on geological grids:** depth downward (`zreversed`); Makie's
  `limits!` with an ascending z range undoes it.
