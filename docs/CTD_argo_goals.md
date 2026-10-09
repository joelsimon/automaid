> **AI-written documentation — non-authoritative.** This page records project goals, not an Argo approval or format specification. Confirm technical requirements with the current Argo User's Manual and a recognized Argo DAC before relying on them.

# Argo pathway goals

## Short term: example-shaped NetCDF prototype

Produce one NetCDF file from an automaid CTD profile that resembles the
structure and conventions of comparable, real Argo profile files retrieved
from a GDAC. Compare against actual files, rather than building from the
semicolon-separated CSV or from an informal field list alone.

This is a **pseudo-compliant prototype for inspection and workflow
development**. It is not an official Argo product, is not validated as
submission-ready, and must not be sent to a GDAC as an Argo profile. Use a
filename and metadata that make its prototype status unmistakable; do not
borrow a real float's WMO identifier or imply that automaid is an Argo DAC.
Record missing or invented values explicitly instead of presenting them as
observations.

The short-term scope is one representative profile. Match the chosen example's
applicable format generation, dimensions, variable layout, attributes, and
QC-field shape closely enough to compare them meaningfully. Document which
values came from MERMAID data and which are placeholders. The exact target
profile type and format version remain to be selected after examining real
GDAC examples and checking with Argo documentation.

See the [short-term checklist](CTD_netcdf_checklist.md) for completion steps.

## Long term: establish an official Argo data pathway

Determine whether MERMAID platforms and their CTD observations qualify for an
Argo contribution, then establish the institutional and operational
arrangements required to maintain that contribution. This work precedes any
claim that an automaid output is archive-ready.

The long-term goal includes:

- Agreeing on eligibility, an Argo mission and platform identity, and the
  appropriate core-profile or other data pathway with the Argo program and a
  recognized DAC.
- Naming a PI and responsible contacts for the platform's full data lifecycle.
- Establishing telemetry and delivery, required metadata and technical files,
  real-time QC, delayed-mode QC, long-term curation, and a funded responsible
  group for each approved parameter.
- Implementing the DAC-agreed production format, provenance/history,
  correction process, and routine submission to both GDACs.
- Validating the end-to-end process with the DAC before operational release.

Argo's guidelines require an agreed data pathway, quality-control
responsibilities, and long-term curation arrangements before a float becomes
part of the Argo system. A file that resembles a GDAC example does not satisfy
those program requirements.

## Primary references

- [Argo float guidelines](https://argo.ucsd.edu/about/what-makes-a-float-part-of-argo/table-of-guidelines-for-argo-floats/)
- [Argo User's Manual: DAC–GDAC data management](https://oneargo.github.io/argo-format-user-manual/chapter6.html)
- [Argo User's Manual](https://oneargo.github.io/argo-format-user-manual/)
- [ADMT tools for DACs](https://www.argodatamgt.org/DACs.html)
