> **AI-written documentation — non-authoritative.** Verify terminology against the cited source and SEED manual before editing code or comments.

# Clock terminology cleanup TODO

Read-only audit: 2026-10-29, automaid commit
`3d183c1636f92c7fe8619b02adc5892500dbd208`.
These are terminology/documentation tasks, not authorization to change code.
Preserve the existing clock algorithms, signs and values. The selected
[GeoCSV definitions](geocsv-method-identifiers.md) remain authoritative.

## Definitions to preserve

```text
manual drift = GPS time - MERMAID time
ClockOffset = MERMAID time - GPS time = gps.mseed_time_delay
TimeCorrectionApplied = event.mseed_time_correction
corrected_start = uncorrected_start + TimeCorrectionApplied
```

The GPS-fix offset and recording correction refer to different epochs. At the
same epoch, the applied correction has the opposite sign to the estimated
offset. Positive fixed-header correction moves a timestamp later; when already
applied, the MiniSEED applied bit prevents a second application.

## Cleanup checklist

- [ ] In [gps.py](../scripts/gps.py), replace the claim of a generic "miniSEED
  convention of a time delay" with the explicit automaid clock-offset definition.
  `mseed_time_delay` is an internal variable, not a named fixed-data-header field.
  Do not confuse it with filter-stage delay in SEED Blockette 57.
- [ ] In [events.py](../scripts/events.py) `correct_clockdrift()` and the analogous
  [cycles.py](../scripts/cycles.py) wording, distinguish MERMAID clock error
  measured against GPS from drift of the GPS clock itself. Prefer "reads
  behind/ahead" to "slow/fast" where discussing accumulated offset rather than
  oscillator rate. Keep the current sign statements and calculations.
- [ ] In [main.py](../scripts/main.py), narrow the GPS-validation comment that
  mentions correcting diving/surfacing timestamps. The called clock-correction
  routine updates recording objects; pressure-derived crossing times do not
  receive that correction.
- [ ] In [the timing notes](../notes/set_mseed_time_correction.txt), replace the
  ambiguous "ADVANCES" with "moves the timestamp later" for positive fixed-header
  correction. Contrast explicitly with Blockette 57, whose positive correction
  denotes time advance in the opposite direction. The existing equations are
  consistent; the word alone is ambiguous.
- [ ] In the [abandoned header experiment](../tests_and_verifications/mseed_fixed_section_data_header.py),
  correct the explanatory table calling `LONG` unsigned. SEED defines it as a
  signed 32-bit quantity. The production `struct` packing already uses the
  signed type. Retain the experiment's archival status; do not execute it on
  existing data.

## Evidence and scope

The [SEED v2.4 manual](https://www.fdsn.org/pdf/SEEDManual_V2.4.pdf), printed
pages 33, 77 and 109, defines signed LONG, the decimation convention, and the
fixed-header correction units/applied bit. The fixed-header paragraph does not
explicitly state an addition equation; the addition rule was verified against
ObsPy/libmseed source and positive/negative synthetic records with both flag
states. Do not attribute an unstated equation to the manual.

The separate missing-record-offset defect was fixed in automaid commit
`c16bf70`; current-helper synthetic checks confirm correction metadata on every
record. The signed five-second guard remains an inspection item. Existing
MiniSEED correction truncation is accepted; no rounding change is planned.
These are separate from the terminology checklist. The audits changed no code,
source data or existing waveform outputs.
