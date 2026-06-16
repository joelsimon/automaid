automaid

Python tools for processing MERMAID instrument data.

automaid converts raw MERMAID transmissions into processed scientific products, including:

* data classification and organization
* clock-drift correction
* float-position interpolation
* seismic SAC and miniSEED generation
* event and dive visualizations
* KML export products

Authors

* @joelsimon
* @sebastienbx
* @oseanfro

Installation

Create a virtual environment and install automaid in editable mode:

python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"

Native dependencies

The legacy processing workflow requires two platform-specific wavelet inversion executables.

Build them locally and copy them into scripts/bin/:

cd scripts/src/V103
make
cp icdf24_v103_test ../../bin/
cd ../V103EC
make
cp icdf24_v103ec_test ../../bin/
cd ../../..

Required executable paths:

scripts/bin/icdf24_v103_test
scripts/bin/icdf24_v103ec_test

Usage

Run the package CLI:

automaid process \
    --server /path/to/server \
    --processed /path/to/processed \
    --database /path/to/database

Equivalent short flags:

automaid process -s /path/to/server -p /path/to/processed -d /path/to/database

Before processing begins, the CLI verifies:

* input paths exist
* output directories exist or can be created
* required native executables are present and executable

Legacy entry point

The original script entry point remains available:

python scripts/main.py \
    --server /path/to/server \
    --processed /path/to/processed \
    --database /path/to/database

If the MERMAID environment variable points to a directory containing:

server/
processed/
database/

the legacy workflow can use those locations as defaults.

### Limitations
Some legacy processing options remain configured within the legacy workflow and have not yet been exposed through the CLI.

### Citation

Simon, J. D., Bonnieux, S., Rocca, F., Simons, F. J., & The EarthScope-Oceans Consortium (2024).

earthscopeoceans/automaid: v4.0.0.

Zenodo.

https://doi.org/10.5281/zenodo.5057096
