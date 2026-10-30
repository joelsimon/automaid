# P0006 historical requested-recording sample frequencies

The August 17, 2018 file `06_5B773AE6.MER` contains four historical requested recordings. Comparing their `<EVENT><INFO TRUE_FS=...>` values with `<ENVIRONMENT><TRUE_SAMPLE_FREQ FS_Hz=...>` in the MER files associated with their original acquisition dives gives:

| `<EVENT><INFO DATE=...>` in `06_5B773AE6.MER` | Original acquisition-dive MER file | `<EVENT><INFO TRUE_FS=...>` | Original acquisition-dive `<ENVIRONMENT><TRUE_SAMPLE_FREQ FS_Hz=...>` | Equal? |
| --- | --- | --- | --- | --- |
| `2018-07-06T01:49:29` | `06_5B3F1904.MER` | `40.014095` | `40.014125` | No |
| `2018-07-13T09:48:32` | `06_5B4C1D4F.MER` | `40.014095` | `40.014104` | No |
| `2018-07-28T22:56:20`, `<EVENT><INFO SMP_OFFSET=9618902>` | `06_5B62EA1F.MER` | `40.014095` | `40.014107` | No |
| `2018-07-28T22:56:20`, `<EVENT><INFO SMP_OFFSET=9618905>` | `06_5B62EA1F.MER` | `40.014107` | `40.014107` | Yes |

The original acquisition-dive associations were checked against saved log intervals, and the frequencies were read directly from the original MER files.

Detected recordings in those older MER files have no `<EVENT><INFO TRUE_FS=...>` attribute; their `<EVENT><FORMAT SAMPLING_RATE=...>` value is `20.000000`.

The August 17 transmission file itself contains `<ENVIRONMENT><TRUE_SAMPLE_FREQ FS_Hz=40.014087 />`.

Three historical requested recordings therefore repeat `<EVENT><INFO TRUE_FS=40.014095>` across different acquisition dives without matching those dives' `<ENVIRONMENT><TRUE_SAMPLE_FREQ FS_Hz=...>` values. The reason remains unresolved.

## Local evidence

The four MER files in this directory are byte-identical copies from `/Users/jdsimon/mermaid/server/`. Original files remain unchanged. `source-manifest.json` records absolute source paths, resolved source paths, byte counts and SHA-256 hashes; `SHA256SUMS` permits local copy verification.

`historical-cycle-frequency-check.json` preserves the previously checked acquisition-dive log intervals and comparison results. The saved cycle associations came from `/Users/jdsimon/mermaid/processed/452.020-P-06/452.020-P-06.pickle`; these associations are contextual evidence, separate from the literal frequency values in the copied MER files.
