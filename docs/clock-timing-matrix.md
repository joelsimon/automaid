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
| MiniSEED fixed header | Field 16: correction, in 0.0001-second ticks | `trunc((t_gps - t_mer) / 0.0001)` | Negative or zero after truncation | Positive or zero after truncation |

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
- [MiniSEED serialization](../scripts/utils.py): truncation to 0.0001-second ticks, applied-bit setting and per-record correction metadata.
- [SEED v2.4 manual](https://www.fdsn.org/pdf/SEEDManual_V2.4.pdf), printed pages 33 and 109: signed LONG representation, correction units and applied-bit semantics. The fixed-header addition rule was separately verified in ObsPy/libmseed source and synthetic records; it is not an explicit equation in that manual paragraph. Blockette 57 on page 77 uses a distinct filter-correction convention.

Source behavior was inspected at automaid commit `9888d0015ca15c2c4c0d14b75b18d9df656b174a`.
The selected names, row requirements and accepted quantization are recorded in
[the migration decisions](geocsv-method-identifiers.md). This matrix changes no
code, clock model, timestamps or scientific data.

## Local migration verification

Automaid v4.6.0 changes only GeoCSV serialization and its allowed version value.
ClockOffset still reads `gps.mseed_time_delay`; TimeCorrectionApplied still reads
`event.mseed_time_correction`. Recording Time still comes from the already corrected
source Stats. Six-digit UTC serialization uses nearest-microsecond ties to even;
MiniSEED truncation, applied bits and record-offset logic remain unchanged.
