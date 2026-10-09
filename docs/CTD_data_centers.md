> **AI-written documentation — non-authoritative.** Verify this summary against the official Argo Data Management Team documentation and the linked services before relying on it.

# Where Argo files go

## Authoritative public archive

For Argo data, the counterpart to EarthScope's seismic data archive is the
pair of **Argo Global Data Assembly Centres (GDACs)**:

- [US GDAC (USGODAE)](https://www.usgodae.org/argo/argo.html), operated in the
  United States by the Fleet Numerical Meteorology and Oceanography Center
  (FNMOC) in Monterey, California.
- [Coriolis GDAC](https://data-argo.ifremer.fr/), operated by Ifremer/Coriolis
  in Brest, France.

Both GDACs distribute the Argo collection managed under the Argo Data
Management Team (ADMT), including profile, trajectory, metadata, and technical
files. They are the standard source when users need the official Argo NetCDF
files. The two centres synchronize their holdings.

Primary overview: [Argo data from GDACs](https://argo.ucsd.edu/data/data-from-gdacs/).

## How the archive is organized

GDAC servers offer equivalent data through HTTP/HTTPS and FTP (and support
other access methods). The principal directories are:

- `dac/`: files grouped by the Data Assembly Centre responsible for the float.
- `geo/`: files grouped geographically by ocean basin.
- `latest_data/`: recently updated files.
- `aux/`: auxiliary, experimental sensor data; these files are distributed
  but are not curated as core Argo data.

For example, a Coriolis DAC profile path may look like:

```text
pub/dac/coriolis/4902602/profiles/R4902602_001.nc
```

Here `coriolis` is the DAC, `4902602` is the float's WMO identifier,
`profiles/` contains cycle profile files, `R` identifies real-time mode, and
`001` is the cycle number. The actual path and file available depend on the
float's DAC and GDAC contents; this is an illustrative example from Argo
documentation.

Top-level index files list available files and useful metadata, allowing users
to search by such properties as time, location, and DAC. Monthly snapshots of
the GDAC collection are also archived with DOIs for citation and reproducible
research. See the [GDAC access guide](https://argo.ucsd.edu/data/data-from-gdacs/)
for current links and retrieval options.

## OceanOPS and other portals

[OceanOPS](https://argo.ucsd.edu/oceanops/) monitors the Argo observing array
and provides float-network information, deployment and sensor metadata, status,
and maps. It is useful for finding and understanding floats; it is **not one of
the two GDACs and is not the canonical archive for Argo profile files**. After
identifying a float in OceanOPS, obtain its standard profile files from a GDAC.

Other services, including Argovis, ERDDAP, and selection tools, can provide
search, visualization, APIs, or converted downloads. These are useful access
routes, but the GDACs remain the reference for the official Argo NetCDF file
collection.

## Primary references

- [Argo data access: GDACs](https://argo.ucsd.edu/data/data-from-gdacs/)
- [Argo data sources](https://argo.ucsd.edu/data/)
- [Argo Data Management Team documentation](https://www.argodatamgt.org/Documentation)
- [Argo User's Manual, online edition](https://oneargo.github.io/argo-format-user-manual/)
- [Argo User's Manual DOI](https://doi.org/10.13155/29825)

The User's Manual is linked here rather than copied into this repository; use
the ADMT documentation page or DOI landing page for the maintained manual and
related reference tables.
