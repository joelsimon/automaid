> **AI-written documentation — non-authoritative.** Verify scientific and operational details against the linked contract and current source before implementation.

# GeoCSV migration plan

Status: schema 1.0.0 implemented locally on 2026-10-29 with automaid v4.6.0,
consumer 0.27.0 and publisher 0.4.0. Current producer source is v4.6.2. The
canonical processed corpus was regenerated on 2026-10-30; fresh testing-only
acceptance is required for those bytes. These are source versions, not claims
of software release or compatible hosted data. The numbered steps below retain
the migration procedure; they do not imply that implementation is still pending.
Publication, deployment, source-archive modification,
commits, tags and releases require separate authorization. The accepted
[format decisions](geocsv-method-identifiers.md) and their
[consumer companion](../../mermaid-geocsv/docs/geocsv-method-identifiers.md)
define the replacement profile. **Legacy MERMAID GeoCSV versions are
unsupported.** No legacy reader/writer path, aliases, automatic fallback or
active legacy compatibility test matrix is retained. `mermaid-geocsv` supports
only the controlled MERMAID GeoCSV profile; a generic GeoCSV reader is outside
scope. This is a complete format-specification review, not only a label rename.

The prior design discussion was documentation-only. The separately authorized
implementation used the writer/version allowlist below and focused
consumer/publisher/site changes. Dependency reinstalls were authorized for that
historical implementation; the current validation task uses existing environments
without modifying them. Protected products, mirrors and historical archives remain
read-only.
Finite signed Pressure-Absolute gauge readings are preserved; Recording pressure
and crossing proxies remain nonnegative. No scientific matching tolerance changes.

Fresh regenerated-corpus evidence is recorded in [the 2026-10-30 validation findings](../../mermaid-publisher/docs/geocsv/regenerated-corpus-validation.md):
43,365/43,365 pairs passed report schema 4 and its acceptance gate. This local
testing-only result does not authorize publication or accept different bytes.

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
- Review the selected normative [payload specification](../../mermaid-geocsv/docs/format-spec.md),
  now written and referenced by the producer, reader, publisher contract and
  site. It includes the complete field dictionary, row matrix, exact header and
  synthetic example. Keep payload semantics distinct from the publisher-owned
  manifest/download protocol; review the complete artifact before code changes.
- Apply the selected file structure: one MERMAID instrument and one CSV table
  per file, the fixed 16-column order in the paired decisions, UTF-8 encoding,
  comma delimiter, LF writer line endings and standard double-quote CSV escaping.
  Specify exact metadata/header rules, types, units, per-row required/applicable
  fields, five method/source identifiers, missingness and six-digit timestamps.
  Document timing provenance, already-applied correction signs, the intentional
  crossing fallback, coordinate/vertical conventions, absolute pressure and the
  approximate station-depth projection. Require Recording-only Elevation as the
  literal numeric cell `0`, type `float` and unit `meters`, denoting the adopted
  sea-level reference. Never emit `0 m`
  as a cell value. Elevation is absent on other row kinds. Reader parsing
  preserves source values and missingness and performs no scientific inference or default filling.
- Require empty cells for fields inapplicable to a row kind; populated
  inapplicable fields raise. Optional applicable fields may be empty, without
  defaults. Location is required only on Recording rows and must be empty
  on Position, Pressure and Crossing rows. The selected row matrix is recorded
  in the paired format decisions; no Location applicability decision remains open.
- Emit `InstrumentDescription` as `MERMAID_<model>`, preserving the
  actual source model code (for example `MERMAID_452.020`). This is
  an explicit string-format change, not a default model or hardware change.
- Require exactly one schema declaration, planned as the first line:
  `# MERMAID GeoCSV schema: 1.0.0` (no trailing period). Keep the existing GeoCSV
  metadata declaration. The schema revision is independent of software versions;
  it does not authorize package 1.0.0 or imply support for any earlier revision.
- Selected column renames are `StartTime` to `Time`, `TimeDelay` to
  `ClockOffset`, and `TimeCorrection` to `TimeCorrectionApplied`. The clock
  columns preserve their existing values and signs. `ClockOffset` is instrument
  minus GPS time (`gps.mseed_time_delay`); the already-applied correction remains
  `event.mseed_time_correction`. Never swap those source variables or negate
  either export during the rename. Reserve drift for GPS minus instrument time.
  `ClockOffset` is optional on Position rows; unknown remains absent.
  `TimeCorrectionApplied` is required on Recording rows; missing raises.
  Each is absent on other row kinds: do not move the measured GPS offset into
  Recording rows, infer a correction from another row, or fill with zero.
  Review every other existing column rather than
  automatically treating the current layout as a finished decade-long contract.
  Do not add a direction column or change scientific algorithms. Apart from
  the agreed explicit Recording Elevation, preserve scientific values. Discuss
  any additional column/name/meaning change before freezing it.
- Apply the selected header contract: required schema, dataset, created,
  delimiter, lineterminator, field_type and field_unit each occur once before
  the column header. Description/attribution are optional; other descriptive
  comments remain opaque source evidence. Reject repeated required declarations
  and field_missing. Follow the specification's types, units, field domains,
  CSV grammar and acceptance/rejection cases without defaults or legacy inference.
- Serialize sub-microsecond source epochs by rounding to the nearest microsecond,
  with exact ties to even and normal carry. This is GeoCSV formatting only;
  preserve clock algorithms. The original migration retained MiniSEED fixed-header
  truncation; the separately authorized v4.6.1 follow-up superseded it.
- **Design gate:** review the complete specification and examples before making
  functional changes. The compact decision note fixes selected decisions; it is
  not the full normative specification. Record genuine remaining ambiguities
  before implementation rather than resolving them silently in code.
- Inventory active references to old headers, method identifiers and missing
  markers across all four repos. Historical reports, archives and frozen examples
  remain historical; distinguish them from active examples and requirements.

### Reconciled pre-migration contract conflicts

The consumer AGENTS now states the selected target rather than requiring old
`StartTime`, empty-known-Location, built-in `nan` missingness or an active 0.24.0
baseline. Keep this instruction audit part of specification review.

The consumer's reader, association and SAC references, public docstrings and
tutorial pairs now describe the implemented schema, including timestamp,
missingness and source-elevation rules. Pre-migration reports and frozen files
remain historical evidence, not active compatibility requirements. Generic
GeoCSV examples are retained under the consumer's
[specification references](../../mermaid-geocsv/docs/references/geocsv-v2.0.4-examples/README.md),
not its data fixtures. The [fixture inventory](../../mermaid-geocsv/data/fixtures/README.md)
distinguishes active schema conversions from frozen historical evidence.
The publisher AGENTS makes the replacement's new-only policy explicit;
its contract and the site link to the locally implemented payload specification.
Compatible hosted payloads remain a separate pending publication step. Publisher 0.4.0 removes active historical-consumer and encoded-path fallbacks
while retaining frozen evidence; report schema 4 uses `recording_records`.
The older implementation descriptions do not override the new payload contract.

### Clock audit dispositions

The original migration retained MiniSEED fixed-header truncation to
0.0001-second units. In the normal already-corrected write path this records
the applied amount. That migration did not include a rounding refactor.
The separately authorized v4.6.1 follow-up replaces that truncation with
nearest-tick decimal ties-to-even rounding and validates finite/signed-32-bit
values before file mutation; see [the clock timing matrix](clock-timing-matrix.md).
That follow-up did not regenerate protected products or migration candidates.
Automaid v4.6.2 adds independent per-record read-back of the applied flag and
rounded correction; it changes no clock calculation or GeoCSV value. Canonical
regeneration is complete as of 2026-10-30; earlier candidate acceptance remains
historical and cannot accept the regenerated bytes.

The all-record correction-value bug is resolved by automaid commit `c16bf70`,
which restored record-offset advancement. Current-helper checks passed for
positive and negative corrections across seven records each, with samples and
timestamps unchanged. Do not list it as unresolved or make further utility
changes under the GeoCSV migration's writer-only allowlist.

## 2. Update the allowed producer and consumer together

- In a separately authorized implementation, `automaid/scripts/geocsv.py`
  emits the five agreed identifiers, `Time`, `ClockOffset` and
  `TimeCorrectionApplied`. Preserve the existing clock source variables,
  signs and all correction algorithms. Also emit
  six-digit UTC timestamps (including `#created`) with the selected rounding,
  true empty missing cells and two-space known blank Location. Emit the literal
  Elevation cell `0` on Recording rows, with `float` type and `meters` unit;
  keep it absent elsewhere. Preserve source sensor/version provenance,
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

Publisher 0.4.1 completes the creation/inspection identity preflight: require
the exact schema 1.0.0 first line and literal six-digit UTC `#created`
declaration, including recordless sources; reject duplicate schema declarations
and comments after the column header. Quoted continuation lines remain CSV
data. Full parsing stays consumer-owned. Report schema 4, scientific limits,
manifest/download semantics and hosted deployment status are unchanged.

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

- Exercise the complete specification example and its acceptance/rejection
  cases: all five identifiers; exact timestamp width, UTC, ties-to-even rounding
  and second/date carry; unchanged fallback semantics; numeric corrections;
  literal string `nan`; numeric rejection; missing required fields; and known blank versus unknown
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

The complete [schema contract](../../mermaid-geocsv/docs/format-spec.md) is
written with the selected names, clock signs, file/header grammar, field
requirements/domains, missingness, literal Recording Elevation `0`, six-digit
timestamp rounding, intentional fallback and complete synthetic example.
No naming, sign, row-applicability or required-header decision remains open.
The selected profile is implemented. Review the complete artifact before any
additional functional change. Current source versions are automaid v4.6.2,
consumer 0.27.0 and publisher 0.4.0; schema
`1.0.0` is a separate specification identifier. Any additional scientific/schema
decision, or automaid code edit outside the allowlist, requires discussion.
Do not describe the replacement as ready until design review and the integrated
end-user acceptance path have both passed.
