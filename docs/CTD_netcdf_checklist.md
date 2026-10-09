> **AI-written documentation — non-authoritative.** This checklist is for a non-submission prototype. Verify details against the selected real GDAC examples, the current Argo User's Manual, and automaid code before implementation.

# Short-term Argo-shaped NetCDF checklist

Goal: create one automaid-derived NetCDF profile that is structurally
comparable to a real Argo profile downloaded from a GDAC. This is a prototype,
not an official Argo profile or an archive submission. Check items off as
completed; record decisions and unresolved fields in the prototype notes.

## Choose and inspect reference files

- [ ] Select a real, current core profile from a GDAC that is as comparable as
  practical to the automaid sensor/profile (sensor type, one profile per file,
  pressure range, and available variables).
- [ ] Save the source GDAC URL, WMO/cycle filename, retrieval date, format
  version, and checksum alongside development notes. Keep the downloaded file
  as an unchanged reference, not an automaid output.
- [ ] Inspect the reference with a NetCDF-aware tool. Record dimensions,
  variables, types, attributes, fill values, units, QC fields, and any raw,
  adjusted, calibration, history, or position fields relevant to the selected
  format version.
- [ ] Select and document the target profile type and format generation. Note
  differences between that target and the chosen reference, especially if the
  reference uses a multi-profile or newer format structure.

## Map automaid observations

- [ ] Identify the exact source profile and sensor model in automaid; retain
  traceability to its source files and processing steps.
- [ ] Map pressure, temperature, and salinity to Argo's corresponding
  parameter variables, and verify units, conversions, scaling, missing values,
  and sample ordering against the source data and sensor documentation.
- [ ] Identify the observation time and the defensible profile-level position.
  Record how the position was obtained and its quality. Do not fabricate a
  position for each depth level.
- [ ] Identify which Argo metadata fields have trustworthy source values and
  which are unavailable. Do not invent a WMO identifier, DAC, PI, deployment,
  calibration, or QC result.
- [ ] Define prototype-only values and missing-field handling. Make the
  filename and metadata clearly identify the file as a test product that must
  not be mistaken for an official Argo profile.

## Create and inspect the prototype

- [ ] Write one NetCDF file using the selected reference and format version as
  structural guides; preserve the reference file unchanged.
- [ ] Reopen the output and inspect dimensions, variable names and shapes,
  units, time/position, fill values, and file attributes.
- [ ] Compare the output systematically against the reference. Document
  intentional differences, missing required Argo information, and any
  prototype-only metadata.
- [ ] Confirm that values round-trip correctly and agree with the originating
  automaid profile within the expected precision. Check pressure/temperature/
  salinity alignment and the profile's ordering.
- [ ] Run the OneArgo format checker only if its current documentation says it
  supports the selected profile type and format version. Record the checker
  version and result; a pass is not DAC acceptance or scientific QC. If the
  checker does not support profile files, record that limitation and have the
  prototype reviewed manually against the selected reference and manual.
- [ ] Add the prototype's path, reference provenance, generation method,
  comparison notes, known limitations, and review outcome to this checklist or
  a linked development note.

## Completion criteria

- [ ] A single prototype NetCDF file can be generated reproducibly from a
  documented automaid input profile.
- [ ] A reviewer can compare it with a documented real GDAC file and trace its
  observational values back to automaid sources.
- [ ] Prototype-only fields and gaps are explicit, and the file cannot
  reasonably be mistaken for a validated Argo submission.

Completing this checklist does **not** complete the long-term Argo pathway.
See [Argo pathway goals](CTD_argo_goals.md) for DAC engagement, eligibility,
quality control, and hosting responsibilities.
