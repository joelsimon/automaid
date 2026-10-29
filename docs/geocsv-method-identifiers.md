> **AI-written documentation — non-authoritative.** Verify scientific and format details against the linked sources before relying on them.

# GeoCSV migration decisions

Paired design note: [mermaid-geocsv](../../mermaid-geocsv/docs/geocsv-method-identifiers.md). Keep both notes synchronized. Execution plan: [four-repository migration](geocsv-migration-plan.md).

Status: planned replacement format, recorded 2026-10-29; not implemented by
this document. Legacy MERMAID GeoCSV versions are unsupported: no legacy
reader/writer path, aliases or automatic fallback. The reader supports this
controlled MERMAID profile, not arbitrary GeoCSV dialects or schemas.

The planned first header line is exactly:

```text
# MERMAID GeoCSV schema: 1.0.0
```

Emit it once, without a trailing period, alongside the existing GeoCSV metadata.
The reader requires the supported schema declaration and rejects missing or
unsupported revisions clearly. Preserve it in parsed and association provenance.
Schema versioning is independent of producer/package versions and release IDs;
this does not authorize a software 1.0.0 release or require legacy support.

## Agreed names

`MethodIdentifier` has exactly three nonempty colon-separated components:
`Method:Subject:Source`.

```text
Measurement:Position:<GPS-sensor>
Measurement:Pressure-Absolute:<pressure-sensor>
Derivation:Recording:automaid-<version>
Interpolation:Crossing-Descent:automaid-<version>
Interpolation:Crossing-Ascent:automaid-<version>
```

The source identifies the actual sensor or producer version. Parse colons first;
hyphens within source names remain part of that source. Both crossing subjects
may be grouped as `Crossing`. No direction column is added.

Rename `StartTime` to `Time`. `Recording` describes an original waveform segment,
including its corrected first-sample epoch, station position, sample rate and
count. `Position` is geographic station position; `Location` remains the SEED
location code. Rename `TimeDelay` to `ClockOffset` and `TimeCorrection` to
`AppliedClockCorrection`, with their existing values, signs and second units
unchanged. Other column names and units remain unchanged unless separately
specified. Absolute-pressure labeling does not introduce atmospheric correction.

In the writer, **Recording rows only** emit `Elevation = 0 m`, explicitly
representing the adopted sea-level reference rather than the submerged sensor's
physical elevation. Other row kinds keep Elevation absent. The reader preserves
both the explicit zero and absence; it never supplies a default. Association's
approved elevation projections use the Recording source value instead of
inventing zero. Remove the existing hard-coded zero for SAC `stel` and coordinate
elevation; missing required source Elevation raises. This replaces the current
missing Recording Elevation, not the writer's scientific depth calculation.

Repository-wide rule: never substitute defaults for missing, unknown, invalid
or unavailable source metadata in parsing, association, projections or
provenance. Optional values remain missing; required values raise. Preserve
explicit source values without silent repair or filling from another record.

## Clock definitions and source mappings

Reserve clock "drift" for the manual's DRIFT/GPSACK definition:

```text
drift = gps_time - mermaid_time
ClockOffset = mermaid_time - gps_time = -drift
```

`ClockOffset` is the observed-minus-reference time difference at a GPS fix,
in seconds: positive means the instrument clock is ahead. Its current source
is `gps.mseed_time_delay`, which is already `-gps.clockdrift`. The proposed
`ClockDrift` column is withdrawn; no GPS-column sign reversal is selected.

`AppliedClockCorrection` is `event.mseed_time_correction`, the clock correction
already added to the uncorrected original recording start, in seconds:

```text
Time = uncorrected_starttime + AppliedClockCorrection
AppliedClockCorrection = -estimated_ClockOffset_at_recording_time
```

Positive correction moves the recording timestamp later; negative moves it
earlier. The correction uses the drift sign, opposite to the offset evaluated
at that same epoch. GPS-fix offsets and recording corrections refer to different
epochs and need not have equal magnitudes. Never apply the correction again.

The current source computes `clockdrift_correction = ascent.clockdrift * pct`,
then sets `mseed_time_correction = clockdrift_correction` and adds it to the
uncorrected start. These are definitions of existing behavior, not changes to
clock processing. The source mapping was checked at automaid commit
`3d183c1636f92c7fe8619b02adc5892500dbd208`: 24 extracted-source sign/interpolation
assertions passed, and all 3,266 P0006 GPS GeoCSV rows matched the negative of
the raw GPS export's `gps_time-mer_time` value at printed precision.

## Time and scientific meaning

All timestamp strings emitted by the GeoCSV writer, including `Time` and
`#created`, use UTC with exactly six fractional digits and `Z`. Time offsets and
corrections retain six decimal places in seconds. Displayed digits make no
claim about clock accuracy, measurement precision, or uncertainty.

Recording times retain the existing clock-drift correction. GPS, pressure and
crossing times do not acquire that recording correction. Preserve the existing
crossing fallback: the GPS helper can retain an anchor position and its epoch
instead of estimating unstable velocity from closely spaced fixes. This
intentional stability behavior is accepted and documented; it is not a
migration defect to fix.

The mixed layer is the upper approximately 50 m in this model. Descent exits
its base; ascent enters it. The boundary is configured, not measured from a
temperature profile. Faster horizontal drift there is empirical, approximately
an order of magnitude, with no imposed velocity ratio. Legacy code's deeper
"mixed layer" means the deep layer; correction outside the allowed writer is
deferred. See [Nolet et al. (2024), pp. 2369, 2372](https://doi.org/10.1785/0220230377).

Geographic metadata describe the station, not earthquake parameters. Retain the
approximate `WaterPressure / 100` station-depth mapping, existing tolerances,
and optional SAC `stla/stlo/stel/stdp` projections; never add `ev*` projections.

## Missing values

- An empty decoded CSV field is absent; required/applicable fields follow the
  row schema. Missing required values raise. Emit absence as an unquoted empty
  field; quoted empty fields also decode as absent.
- Nonempty strings such as `nan`, `NA`, `NULL` and `NaT` are literal values.
  Nonempty invalid numbers or timestamps raise; numeric values must be finite.
- Whitespace-only values are invalid except exactly two ASCII spaces in
  `Location`: `6,"  ",9` means a known blank SEED location code; `6,,9` is absent.
  The DataFrame and raw association source retain both literal spaces. Never
  trim, normalize, or default this source value. Association can compare that
  defined blank code with ObsPy's `""` without rewriting either identity or the
  source cell; absent Location remains unknown and never becomes a blank code.

References: [CSV field syntax](https://www.rfc-editor.org/rfc/rfc4180.html),
[blank SEED identifiers](https://www.fdsn.org/media/message_attachments/2015/07/20/FDSN-WS-Specifications-1.1.pdf).
Blank Location does not universally mean "primary instrument".

## Scope and remaining decisions

The names, schema declaration, six-digit display, missing-value policy,
Recording-only zero Elevation and intentional fallback above are selected.
Before code changes, complete and review the normative field dictionary,
row-specific requirements, metadata/CSV
grammar and representative examples. Existing columns, apart from the selected
`Time`, `ClockOffset` and `AppliedClockCorrection` renames, are the review
baseline, not a substitute for that specification review. No
direction column or generic GeoCSV compatibility layer is added.

Current work is definitions/documentation only: no automaid code edits,
including writer or version edits. For a separately requested implementation,
the earlier allowlist remains `scripts/geocsv.py` plus the version value in
`scripts/setup.py`; that boundary is not authorization to implement now.
Stop and discuss before any other code edit.

Historical source data, published bytes, release identities and evidence remain
unchanged. Six-digit formatting applies to new output; fetch must preserve any
existing literal release ID and all its fractional digits. No source rewrites,
uploads, pointer changes, deployment, commits, tags or releases are authorized.

Source references: [writer](../scripts/geocsv.py),
[position helper](../scripts/gps.py),
[cycle calculation](../scripts/cycles.py).
Timing and fallback behavior were inspected at automaid revisions
`678ea8f06e445d3da1864236f29a2ebee87ee659` and
`f0c97e98629967de8161f744e9fbbc4474030b9b`; inspect current source before implementation.
