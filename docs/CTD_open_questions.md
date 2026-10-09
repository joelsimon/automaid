> **AI-written documentation — non-authoritative.** This note records a terminology and interpretation issue for project discussion. Consult the current Argo User’s Manual and the applicable DAC guidance before implementing Argo-like metadata.

# CTD open questions

## What does Argo `JULD` mean for a profile?

Argo profile files store profile-level `JULD` (profile date/time), `LATITUDE`,
`LONGITUDE`, and `JULD_LOCATION` (the time associated with the reported
position). These fields describe a profile; they are not per-pressure-level
timestamps or positions.

A footnote in the Argo User’s Manual says to assume an ascending profile. In
that context, the float begins transmitting when it first reaches the surface,
and the earliest transmission time is recorded as `JULD`. If that first
transmission cannot be used to derive a position, `JULD_LOCATION` is later: it
marks the time associated with the usable location. This can be misread as
`JULD` being the time ascent starts or the time of the first CTD sample. The
footnote describes a surface-transmission convention, not either of those
events.

Do not generalize that footnote to every float or telemetry system without
checking the applicable timing rules. Argo also defines separate cycle timing
variables for ascent start/end, transmission start, first message, and first
location; their order can differ by system. In particular, the Argo cycle
timing guidance says Iridium floats typically obtain a GPS position after
surfacing and before starting satellite communication. Thus, “first
transmission” and “first location” need not be the same event. `JULD_LOCATION`
is the time associated with the selected profile position; it need not equal
`JULD`.

### Interim automaid decision

The temporary `.ctd` companion will use the time and coordinates from the
first usable GPS fix after the CTD profile. These values will be written in
commented headers; the existing three-column `.csv` remains unchanged. The
timestamp is explicitly labeled as the GPS-fix time and is not claimed to be
Argo `JULD`. This project-specific format is not an Argo implementation.

### Questions to resolve later

- Does `first_gps_fix_after_profile()` select the same fix as `gps_after_dive[0]`?
  That is the expected behavior, but verify it against how both selections
  filter and order cycle GPS fixes before relying on the equivalence.
- What event does the timestamp embedded in each SBE41/SBE61 CSV filename
  represent: profile acquisition, first surfaced telemetry, file creation, or
  something else? Resolve this before attempting to map the filename timestamp
  to an Argo profile-time field.
- If there is no usable post-profile GPS fix, should a later format include a
  clearly marked estimate with method and quality information?
- What metadata and position-quality conventions should a future NetCDF
  prototype use?

## References

- [Argo User’s Manual, profile time and position fields](https://oneargo.github.io/argo-format-user-manual/chapter2.html)
- [Argo User’s Manual, DAC–GDAC data management](https://oneargo.github.io/argo-format-user-manual/chapter6.html)
- [Argo cycle timing variables and event sequence](https://argo.ucsd.edu/how-do-floats-work/argo-cycle-timing-variables/)
- [Argo Data Management Meeting 11 report, common position/time method](https://argo.ucsd.edu/wp-content/uploads/sites/361/2020/04/DM11report.pdf)
