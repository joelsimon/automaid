> **AI-written documentation — non-authoritative.** Verify scientific and format details against the linked source code and the recorded revision before relying on them.

# GeoCSV method identifiers and row times

Paired design note: [mermaid-geocsv](../../mermaid-geocsv/docs/geocsv-method-identifiers.md). Both notes must be maintained together.

## Status and scope

Design discussion recorded on 2026-10-29. The three-component identifier
structure and a common row timestamp are agreed design directions. The field
labels and subject spellings below are proposals awaiting final selection.
Existing producer labels and the `StartTime` column remain implemented behavior.
No format migration or scientific behavior change is implemented by this note.

Geographic values describe the MERMAID receiver (station). The historical
`Algorithm(event)` label refers to recording metadata, not earthquake-source
parameters. A recording row includes the original waveform segment's first
sample time, station position at that time, identity, sampling rate, and sample
count. A later request, download, or clipped window is distinct from that
original segment.

## Proposed identifier terminology

Use `Method:Subject:Source` as the proposed names for the three components of
`MethodIdentifier`:

| Component | Meaning | Reason for this name |
| --- | --- | --- |
| `Method` | How the row's principal values were obtained: `Measurement` or `Interpolation`. | Describes derivation without suggesting a programming-language data type. |
| `Subject` | What the row describes, including the relevant context or epoch. | Covers both a physical quantity such as pressure and a contextual quantity such as position at recording start; `ValueType` could be mistaken for a data type. |
| `Source` | Sensor model or producing software and its version. | Fits both kinds of provenance and avoids the earthquake-related interpretation of `Origin`. |

Identifiers have exactly three nonempty components separated by two colons.
The source component combines software name and version into one token, such
as `automaid-v4.5.14`; it does not contain a colon. That example demonstrates
spelling only and does not claim that the historical version emits new labels.
Actual future output must identify its actual producer version.

Proposed templates:

```text
Measurement:Position:<GPS-sensor>
Measurement:Pressure:<pressure-sensor>
Interpolation:RecordingStartPosition:automaid-<version>
Interpolation:LayerCrossingPosition:automaid-<version>
```

`Position` means geographic station position. The existing `Location` column
has the separate SEED location-code meaning. `RecordingStartPosition` identifies
the spatial purpose of a recording row; its timing is also corrected using an
interpolated clock-drift estimate. The identifier does not imply that waveform
samples or every other column were interpolated.

`LayerCrossingPosition` describes a crossing of the configured layer boundary.
The current calculation uses `mixed_layer_depth_m`; it does not determine the
thermocline depth from a temperature profile in this operation. Selection of
this name, or a more specific boundary name, remains part of the design review.

## Common row timestamp and current producer behavior

Every row emitted by the current MERMAID GeoCSV writer requires `StartTime`.
The planned common name `Time` better describes all four row types: only a
recording row specifically represents a start. The eventual column rename is
a coordinated producer/consumer format change, not an implemented alias.

Time is the row's common reference epoch. It is not necessarily independently
measured: position estimates can depend on the timestamp, and timestamp
provenance differs between row types.

| Row | Meaning of timestamp | Recording clock-drift correction |
| --- | --- | --- |
| GPS position measurement | Timestamp attached to the GPS fix, as supplied by the source parsing path. | Not applied by the recording correction routine. |
| Pressure measurement | Timestamp of the pressure-log entry. | Not applied. |
| Recording start position | Corrected time of the original waveform segment's first sample. | An interpolated drift correction is added before station-position estimation. |
| Layer-crossing position | Intended to be the estimated time of crossing the configured depth boundary; see the fallback finding below. | Crossing time is interpolated from pressure-log timestamps without applying the recording correction. |

For recordings, the current sequence is:

1. Compute the uncorrected first-sample time from recording metadata.
2. Estimate clock drift at that time using the bounding GPS synchronization data.
3. Set `corrected_starttime = uncorrected_starttime + clockdrift_correction`.
4. Request station-position interpolation at the corrected time.
5. Write waveform Stats start time to the GeoCSV recording row.

For layer crossings, estimate the crossing time from pressure values and log
timestamps, then request GPS interpolation at that time. The recording clock
correction routine operates on recording objects and does not correct those
pressure timestamps or crossing times. A timestamp column alone does not imply
identical timing provenance or accuracy across row types.

## Open finding: interpolation fallback can replace the row epoch

The current `gps.linear_interpolation()` helper can replace the requested
interpolation time with an input GPS timestamp:

- A single input fix returns that fix's position and timestamp.
- Retained fixes less than 20 metres apart return an input position and timestamp.

The layer-crossing writer uses the returned position object's `.date`, rather
than the separately calculated crossing time. Consequently, a stored crossing
row can carry a GPS anchor time instead of the intended crossing epoch.
The recording writer takes its timestamp from waveform Stats, so this particular
export path does not replace a recording row's corrected start time.

This is a source-inspection finding. Its incidence in existing data has not
been established. Before implementation, decide how requested row time and
fallback position provenance must be represented. Do not silently relabel an
anchor as an interpolated result or alter historical row timestamps. Resolving
this finding is a scientific behavior decision separate from terminology.

## Coordinated migration and maintenance

Update this note and its linked companion in the same change whenever automaid
changes GeoCSV row construction, timestamp derivation, clock corrections,
position interpolation, or the layer-crossing fallback. Keep the status, table,
findings, and agreed vocabulary synchronized in both repositories. Source-code
inspection references must identify the revision they describe.

The migration must coordinate automaid, mermaid-geocsv, mermaid-publisher, and
bathymetrix-site, including the canonical publisher contract. Compatibility
with historical identifiers and `StartTime`, new-output versioning, fixtures,
and acceptance evidence require an explicit migration plan. Existing source
archives, frozen payloads, manifests, and published release bytes remain
unchanged. This note authorizes no regeneration, publication, or deployment.

## Source inspection

The current-behavior findings above were inspected in automaid commit
`678ea8f06e445d3da1864236f29a2ebee87ee659`. They describe that source revision,
not a newly verified corpus or a guarantee for future producer versions.

- [GeoCSV row construction and required timestamps](../scripts/geocsv.py).
- [Recording clock correction and station position](../scripts/events.py).
- [Layer-crossing calculation and correction scope](../scripts/cycles.py).
- [Position interpolation and timestamp fallback](../scripts/gps.py).
- [Pressure-log timestamp extraction](../scripts/utils.py).
