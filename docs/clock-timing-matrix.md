> **AI-written documentation — non-authoritative.** Verify timing definitions against the linked source code and SEED manual before relying on them.

# Clock timing matrix

Use **`t_mer` = uncorrected MERMAID time** and **`t_gps` = GPS-reference time for the same instant**.

For GPS rows, these describe the synchronization measurement. For Recording
rows, use the original recording's first sample; its `t_gps` is **estimated**
through the existing clock-correction model. The GeoCSV names below are the
locally implemented schema 1.0.0 names; their existing source variables, values
and signs are unchanged.

| Context | Quantity | Definition | MERMAID ahead: `t_mer > t_gps` | MERMAID behind: `t_mer < t_gps` |
| --- | --- | --- | --- | --- |
| GPS | `gps.clockdrift` — manual's drift | `t_gps - t_mer` | Negative | Positive |
| GPS | `gps.mseed_time_delay` | `t_mer - t_gps` | Positive | Negative |
| Position GeoCSV row | `ClockOffset` | `t_mer - t_gps` | Positive | Negative |
| Recording | `event.clockdrift_correction` | Estimated `t_gps - t_mer` | Negative | Positive |
| Recording | `event.mseed_time_correction` | Estimated `t_gps - t_mer` | Negative | Positive |
| Recording GeoCSV row | `TimeCorrectionApplied` | Estimated `t_gps - t_mer` | Negative | Positive |
| Recording | `corrected_starttime` / GeoCSV `Time` | `t_mer + (t_gps - t_mer) = t_gps` | Moves earlier | Moves later |
| MiniSEED fixed header | Field 16: correction, in 0.0001-second ticks | `round_half_even(decimal(t_gps - t_mer) * 10000)` | Negative or zero after rounding | Positive or zero after rounding |

**The GPS and Recording differences refer to different epochs. Their signs
follow the same definitions, but their magnitudes need not match.**

The MiniSEED applied bit is separate from that numeric correction:

| Applied bit | Stored start time, using the first record | Reader behavior |
| --- | --- | --- |
| `0` — not applied | `t_mer` | Adds the stored correction: approximately `t_mer + (t_gps - t_mer) ≈ t_gps`, subject to header quantization |
| `1` — already applied; normal automaid output | `t_gps` | Uses `t_gps` directly; **does not add the correction again** |

A concrete sign check, illustrative rather than measured data:

| Case | `t_mer` | `t_gps` | Manual drift | `ClockOffset` | `TimeCorrectionApplied` |
| --- | --- | --- | --- | --- | --- |
| MERMAID ahead | `12:00:02` | `12:00:00` | `-2 s` | `+2 s` | `-2 s` |
| MERMAID behind | `11:59:58` | `12:00:00` | `+2 s` | `-2 s` | `+2 s` |

## Sources and limits

- [GPS drift/delay](../scripts/gps.py): `mseed_time_delay = -clockdrift`.
- [Recording correction and addition](../scripts/events.py): `mseed_time_correction = clockdrift_correction`, then `corrected_starttime = uncorrected_starttime + clockdrift_correction`.
- [MiniSEED serialization](../scripts/utils.py): nearest-0.0001-second rounding with decimal ties to even, finite/signed-32-bit checks before mutation, applied-bit setting and per-record correction metadata.
- [SEED v2.4 manual](https://www.fdsn.org/pdf/SEEDManual_V2.4.pdf), printed pages 33 and 109: signed LONG representation, correction units and applied-bit semantics. The fixed-header addition rule was separately verified in ObsPy/libmseed source and synthetic records; it is not an explicit equation in that manual paragraph. Blockette 57 on page 77 uses a distinct filter-correction convention.

Source behavior was inspected at automaid commit `9888d0015ca15c2c4c0d14b75b18d9df656b174a`.
The selected names, row requirements and accepted quantization are recorded in
[the migration decisions](geocsv-method-identifiers.md). This matrix documents serialization separately from the unchanged clock model
and already corrected waveform timestamps.

## Local migration verification

Automaid v4.6.0 changes only GeoCSV serialization and its allowed version value.
ClockOffset still reads `gps.mseed_time_delay`; TimeCorrectionApplied still reads
`event.mseed_time_correction`. Recording Time still comes from the already corrected
source Stats. Six-digit UTC serialization uses nearest-microsecond ties to even;
MiniSEED truncation, applied bits and record-offset logic were unchanged by that migration.

## Separately authorized MiniSEED correction follow-up

Automaid v4.6.1 replaces truncation of the fixed-header correction with rounding
to the nearest 0.0001-second tick, using the decimal spelling of the supplied
value and exact ties to even. The maximum quantization error is 50 microseconds.
The rounded tick count must fit the signed 32-bit field; non-finite and
out-of-range corrections raise before file mutation. This changes correction
history serialization only. The clock model, corrected waveform timestamps,
GeoCSV values, applied bits and record-offset advancement are unchanged.
Protected processed products and existing candidate snapshots are not regenerated.


### Conversion audit and verification

The follow-up audit inspected integer conversions and binary header writes in
`scripts/*.py`, plus consumer and publisher timing checks. No second active
float-to-integer timing-header conversion was found. GeoCSV already performs
decimal ties-to-even timestamp rounding. Waveform integer dtypes decode native
samples; source integer parsing preserves integer fields. The rounded measured
sampling frequency in `utils.location` identifies the nominal 40 Hz channel
family and does not quantize the exported sampling rate.

Using conda `pymaid3.10` (Python 3.10.15, ObsPy 1.4.1), two external verification
tests passed: 75 valid-file cases, 21 invalid-file cases and 850 per-record
correction checks across uniform 256-, 512- and 4096-byte big-endian records.
Cases included both signs, exact decimal ties, tick boundaries, signed-32-bit
endpoints, non-finite values and overflow. Only correction bytes and the applied
flag changed; samples, sample counts, sampling rates and start/end times were
preserved. Rejected inputs left every file byte unchanged. Synthetic files
were generated only under `/private/tmp/`; protected products were untouched.

The wider diagnostic also reproduced two pre-existing limits against both the
unchanged helper at `9eddbdddb698cdee63da5d5baa53cb0878fe4a53` and the modified
helper. These are separate from rounding and remain unresolved:

- ObsPy 1.4.1's `set_flags_in_fixed_headers` unpacks headers as big-endian;
  little-endian synthetic input raises in its timestamp decoding.
- Mixed 256/512-byte records are unsafe: the initial record count is inferred
  from file size divided by the first record length. A 16-record synthetic file
  gained 1,068 invalid trailing bytes, although its original record corrections
  were written correctly. Uniform record-length output passed. The normal
  ObsPy writer defaults to big-endian output with a uniform record length.

No flag-setting, record-count or record-offset algorithms were changed in this
follow-up. Support for other header byte orders or mixed record lengths requires
a separately scoped change.

## MiniSEED correction read-back validation

Automaid v4.6.2 adds a separate read-back pass after the correction writer closes
the file. Every record must have the correction-applied activity bit set and its
stored signed integer correction must equal the expected rounded 0.0001-second
tick count, including when that count is zero. A mismatch raises `ValueError`
with the filename, one-based record number, byte offset, expected tick count and
observed flags and tick count. The check independently advances through records
so a missing increment in the writing loop is detected on a subsequent record
whose correction was not written. An exception leaves the modified file in
place; validation does not roll back writes.

Correction rounding, clock modeling, waveform timestamps and GeoCSV output are
unchanged. The existing uniform-record-length and big-endian limitations above
still apply.

Verification used conda `pymaid3.10` and temporary synthetic 256-byte big-endian
records: eight valid correction cases passed 320 raw per-record header checks,
including zero, both signs, decimal ties and signed-32-bit endpoints. Samples,
sample rates and start/end times were preserved. In-memory mutations reproduced
the missing writer increment and an incorrect numeric write; both raised at the
first mismatching record. Suppressing flag writes also raised for a zero
correction. No processed products were modified.
