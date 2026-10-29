> **AI-written documentation — non-authoritative.** Verify scientific and operational details against the linked contract and current source before implementation.

# GeoCSV migration plan

Status: implementation handoff, 2026-10-29; no source changes are implemented
by this planning task. Publication, deployment, source-archive modification,
commits, tags and releases require separate authorization. The accepted
[format decisions](geocsv-method-identifiers.md) and their
[consumer companion](../../mermaid-geocsv/docs/geocsv-method-identifiers.md)
define the replacement profile. **Legacy MERMAID GeoCSV versions are
unsupported.** No legacy reader/writer path, aliases, automatic fallback or
active legacy compatibility test matrix is retained. `mermaid-geocsv` supports
only the controlled MERMAID GeoCSV profile; a generic GeoCSV reader is outside
scope. This is a complete format-specification review, not only a label rename.

The current definition-review task permits documentation changes only.
Do not edit any automaid code, including the writer or version. The later
implementation steps and allowlist below require a separately requested
implementation task; they do not authorize code changes in this discussion.

## Coordination and instruction precedence

Implement one coordinated migration in separate reviewable steps. Update and
verify all affected repositories against one format contract before activating
new hosted data or presenting an incompatible consumer/data combination as
ready. A coordinated migration does not require simultaneous deployments.

The local repositories are `automaid`, `mermaid-geocsv`, `mermaid-publisher`
(not `mermaid-publishing`), and `bathymetrix-site`. They remain independent Git
repositories in the new project. Automaid is the project entry point.

Read global, workspace, shared MERMAID and repository AGENTS before editing.
For this migration, use the applicable scientific-integrity, public-interface,
association, verification and documentation rules in
[mermaid-geocsv/AGENTS.md](../../mermaid-geocsv/AGENTS.md) as the substantive
baseline, alongside each repository's ownership rules. Latest explicit migration
instructions override conflicting older rules. Do not import consumer-only
ownership restrictions into the producer or move responsibilities between repos.

**Automaid code allowlist:** only `scripts/geocsv.py` and the version value
returned by `get_version()` in `scripts/setup.py`. Documentation and AGENTS edits
are permitted. Stop and discuss before changing any other code, including tests,
GPS/cycle/clock calculations, callers, configuration, or output permissions.
Do not bypass this boundary with runtime patching, replacement modules or
changes to another repository that duplicate automaid algorithms.

The consumer's repository-wide no-default rule applies to every source-derived
field: parsing, association, projections and provenance must never substitute
values for missing, unknown, invalid or unavailable source metadata. Optional
values stay missing; required values raise. Explicit source zero or two-space
Location is preserved. Do not silently repair input or fill from another record.

## 1. Establish the working baseline and normative specification

- Inspect each repository's current status, staged changes, commits, installed
  versions and applicable instructions. Preserve unrelated work and staging.
  Use existing virtual environments; do not recreate or upgrade them.
- Read the compact format decisions and the publisher's canonical
  [publication/download contract](../../mermaid-publisher/docs/geocsv/contract.md).
  Replace the active historical-consumer requirement with the new consumer
  requirement. Older AGENTS baseline clauses are superseded for this migration;
  archived evidence is retained as history, not as a compatibility obligation.
- Write one normative payload specification in
  `mermaid-geocsv/docs/format-spec.md`, referenced by the producer, reader,
  publisher contract and site. Keep payload semantics distinct from the
  publisher-owned manifest/download protocol.
- Specify the exact header/order, types, units, per-row required/applicable
  fields, five method/source identifiers, missingness and six-digit timestamps.
  Document timing provenance, already-applied correction signs, the intentional
  crossing fallback, coordinate/vertical conventions, absolute pressure and the
  approximate station-depth projection. Require explicit `Elevation = 0 m`
  only on Recording rows, denoting the adopted sea-level reference, with
  Elevation absent on other row kinds. Reader parsing preserves source values
  and missingness and performs no scientific inference or default filling.
- Require exactly one schema declaration, planned as the first line:
  `# MERMAID GeoCSV schema: 1.0.0` (no trailing period). Keep the existing GeoCSV
  metadata declaration. The schema revision is independent of software versions;
  it does not authorize package 1.0.0 or imply support for any earlier revision.
- Selected column renames are `StartTime` to `Time`, `TimeDelay` to
  `ClockOffset`, and `TimeCorrection` to `AppliedClockCorrection`. The clock
  columns preserve their existing values and signs. `ClockOffset` is instrument
  minus GPS time (`gps.mseed_time_delay`); the already-applied correction remains
  `event.mseed_time_correction`. Never swap those source variables or negate
  either export during the rename. Reserve drift for GPS minus instrument time.
  Review every other existing column rather than
  automatically treating the current layout as a finished decade-long contract.
  Do not add a direction column or change scientific algorithms. Apart from
  the agreed explicit Recording Elevation, preserve scientific values. Discuss
  any additional column/name/meaning change before freezing it.
- Specify encoding, delimiter, quoting, line endings, header/comment rules,
  exact type/domain validation and unsupported-version errors. Include a complete
  new-format example. Reject missing/unsupported schema declarations clearly;
  do not infer a legacy schema or accept unrelated GeoCSV profiles.
- **Design gate:** review the complete specification and examples before making
  functional changes. The compact decision note fixes selected decisions; it is
  not the full normative specification. Record genuine remaining ambiguities
  before implementation rather than resolving them silently in code.
- Inventory active references to old headers, method identifiers and missing
  markers across all four repos. Historical reports, archives and frozen examples
  remain historical; distinguish them from active examples and requirements.

### Known contract conflicts to reconcile

The consumer AGENTS now states the selected target rather than requiring old
`StartTime`, empty-known-Location, built-in `nan` missingness or an active 0.24.0
baseline. Keep this instruction audit part of specification review.

The current `mermaid-geocsv/docs/reader-spec.md`, `docs/obspy-association.md`,
`docs/sac.md`, public docstrings and tutorial pairs still describe implemented
pre-migration behavior. Replace their conflicting timestamp, missingness and
hard-coded elevation rules with the implementation; do not use those older
rules to override the selected target or claim migration is already complete.
The publisher AGENTS and canonical contract still describe active historical
consumer support; update that policy explicitly for the new-only schema.

## 2. Update the allowed producer and consumer together

- In a separately authorized implementation, `automaid/scripts/geocsv.py`
  emits the five agreed identifiers, `Time`, `ClockOffset` and
  `AppliedClockCorrection`. Preserve the existing clock source variables,
  signs and all correction algorithms. Also emit
  six-digit UTC timestamps (including `#created`), true empty missing cells and
  two-space known blank Location. Emit `Elevation = 0 m` on Recording rows
  only; keep it absent elsewhere. Preserve source sensor/version provenance,
  other scientific values, row ordering, deduplication, pressure mapping and the
  intentional position/time fallback. Document corrected layer terminology in
  this file without refactoring off-limits GPS/cycle code.
- In `mermaid-geocsv`, update parsing, typing, metadata and association to the
  new profile. Require its declared schema revision, retain that revision in
  parsed/association provenance, and reject missing/unsupported revisions.
  Recognize recording rows using the new identifier; retain source labels,
  row order, SNCL requirements and complete association provenance.
  Keep quoted-empty absence distinct from two-space known blank Location.
  Retain both Location spaces literally in the DataFrame and raw source row;
  comparison with ObsPy's empty code never rewrites source or waveform identity.
  Absent Location stays unknown. Literal `nan` is valid only where the declared
  type permits that string. Never infer or default scientific source values.
- Source the approved coordinate elevation and SAC `stel` from the explicit
  Recording Elevation, preserving the documented destination-unit conventions.
  Remove the existing hard-coded coordinate elevation and SAC `stel` zero
  fallback. Missing required source Elevation raises; it must never silently
  become zero. Other reader values remain faithfully parsed.
- Preserve matching cardinality, freshness/atomicity, sample coverage, other SAC
  mappings and all scientific tolerances. Do not introduce aliases, automatic
  legacy fallback, earthquake projections or public API expansion.
- Replace active old-format fixtures/expectations with focused tests of
  temporary producer output or explicit new fixtures with documented provenance.
  Preserve archival originals as evidence without keeping them as compatibility
  targets for the active reader. Automaid gets no new test-code files
  without discussion; use existing checks or temporary external verification.
- Bump producer version only through the allowed setup value, and bump consumer
  version in its canonical authority. Use actual producer versions in output.
  Version bumps do not authorize tags or releases.

## 3. Update publication and serving consumers

- In `mermaid-publisher`, update recording selection/counts, `Time` access,
  scientific/export checks, active fixtures, documentation and active consumer
  requirements. Inspect `_verify_sac.py` and `validate_snapshot.py` explicitly.
  Preserve comparisons and tolerances; do not duplicate package association logic.
- Review report field names and schema versions if recording terminology changes
  externally stored fields. Keep old reports and baseline evidence unchanged;
  bind all new acceptance evidence to the actual new bytes/software. Bump the
  publisher version for changed public behavior in its existing version file.
- In `bathymetrix-site`, update format examples, consumer requirements and hosting
  references. The Worker serves exact bytes and may need no code change: verify
  the existing routing/serving contract before editing it.
- Preserve the manifest protocol, literal release provenance, two-colon timestamp
  path transformation, testing isolation, immutable bytes and rejection of old
  encoded paths. Do not add a manifest field or schema-dependent serving behavior
  merely because payload contents changed.

## 4. Verify the integrated replacement locally

- Exercise all five identifiers; exact timestamp width, UTC and chronology;
  unchanged fallback semantics; numeric corrections; literal string `nan`;
  numeric rejection; missing required fields; and known blank versus unknown
  Location through parsing and waveform association.
- Compare scientific values before/after migration, allowing only explicitly
  agreed representation changes. Higher displayed timestamp precision is a
  format decision, not evidence of better clock accuracy. Never widen tolerances
  to make newly generated output pass.
- Run focused parser/association/public-API tests, then the package suite. Update
  README, reader/association references and every affected Markdown/notebook pair;
  verify parity, cleared notebook outputs and top-to-bottom execution. Include
  a plain-language replacement-format notice and a CSV import example that
  preserves two-space Location, leading-zero codes and literal string `nan`.
  Remove current claims/examples implying legacy input remains supported.
- Run publisher unit/integration and SAC verification checks with the current
  installed consumer, plus site Worker tests. Exercise local latest/pinned fetch,
  one/all instruments, memory/disk/reuse, exact bytes and literal provenance.
  Use a local test server and temporary outputs; no remote write is needed.
- Perform fresh complete scientific acceptance for the chosen new candidate's
  original waveform-pair scope, reading source data without modifying it. Count
  the current inputs rather than borrowing historical totals or acceptance.
  Never run a processing entry point that writes protected data or refreshes
  external mirrors. Demonstrate isolated output before any producer run.

## 5. Handoff and separately authorized cutover

- Check diffs against each repository's starting state, including the automaid
  allowlist. Report changed files, actual source/installed versions, exact test
  results, unresolved failures and any locally changed but undeployed behavior.
  Suggest scoped commits; do not stage or commit automatically.
- Keep paired decisions, scientific references, publisher contract, website
  runbook and tutorial companions synchronized. Source changes and future
  software/data publication remain separate outcomes.
- End-user readiness requires one new-format artifact to pass writer, reader,
  association, publisher acceptance, download, parsing and the documented
  tutorial workflow with the same agreed values. Until the compatible data and
  consumer are available together, do not present current old-format hosted
  latest as a working example for the new-only reader.
- A new schema requires new immutable payloads, a new shared `#created`/literal
  release ID, fresh manifests/hashes and fresh acceptance. Do not rewrite existing
  published objects, historical fixtures or original waveform/source archives.
- Only after separate publication authorization: publish compatible software/data
  in a controlled order, verify immutable pinned downloads and parsing, then
  activate testing latest last. Production remains inactive. No deployment,
  upload, deletion, pointer change, archival commit, push, tag or release follows
  automatically from completing this implementation plan.

## Design status and remaining work

The five method names, `Time`, `ClockOffset`, `AppliedClockCorrection`,
six-digit display, missing-value rules, explicit schema
declaration, Recording-only zero Elevation and intentional fallback are
selected. They are not a complete format specification. Phase 1 must finish and review every field/row rule,
file/header grammar and examples before implementation. Exact software version
increments can be chosen during migration within the rules above; schema
`1.0.0` is a separate specification identifier. Any additional scientific/schema
decision, or automaid code edit outside the allowlist, requires discussion.
Do not describe the replacement as ready until design review and the integrated
end-user acceptance path have both passed.
