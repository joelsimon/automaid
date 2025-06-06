# -*- coding: utf-8 -*-
#
# Part of automaid -- a Python package to process MERMAID files
# pymaid environment (Python v3.10)
#
# Developer: Frédéric rocca <FRO>
# Contact:  frederic.rocca@osean.fr
# Last modified by FRO: 09-Sep-2024
# Last tested: Python 3.10.13, 22.04.3-Ubuntu
#
# Developer: Joel D. Simon (JDS)
# Contact:  jdsimon@alumni.princeton.edu
# Last modified: 13-May-2025
# Last tested: Python Python 3.10.15, Darwin Kernel Version 23.6.0

import os
import re
import glob
import pytz
import shutil
import argparse
import datetime
import functools

import kml
import gps
import setup
import cycles
import utils
import events
import vitals
import geocsv
import preprocess
import sbe41
import sbe61

## !!! TEMP VARS !! ##
csv_file = True
## !!! TEMP VARS !! ##

# Log a creation date for metadata files in ISO 8601, milliseconds precision,
# with "Z" suffix for UTC: "YYYY-MM-DDTHH:MM:SS.sssZ"
creation_datestr = datetime.datetime.now(pytz.UTC).isoformat()[:23] + "Z"

# Get current version number.
version = setup.get_version()

# Set depth of mixed layer (in meters) for drift interpolation
mixed_layer_depth_m = 50

# Minimum number unique GPS fixes before/after each dive required for
# interpolation
min_gps_fix = 2

# Maximum-allowed time in seconds before(after) diving(surfacing) to record the
# minimum-required number of unique GPS fixes
max_gps_time = 5400

# Toggle preliminary (rapid) location estimates on and off
preliminary_location_ok = False

# Set automaid scripts' directory path
main_path = os.path.abspath(__file__)
scripts_path = os.path.dirname(main_path)

# Set default processed,server paths and database path
def_mermaid_path = os.environ["MERMAID"]
def_server_path = os.path.join(def_mermaid_path, "server")
def_processed_path = os.path.join(def_mermaid_path, "processed")
def_database_path = os.path.join(def_mermaid_path, "database")

# Parse (optional) command line inputs to override default paths
parser = argparse.ArgumentParser()
# problem: metavar=''   prints: "-s , --server"
# problem: metavar='\b' prints: "-s, --server", but misaligns the help statement...
parser.add_argument('-s',
                    '--server',
                    default=def_server_path,
                    dest='server',
                    #metavar='\b',
                    help="server directory (default: {:s})".format(def_server_path))
parser.add_argument('-p',
                    '--processed',
                    default=def_processed_path,
                    dest='processed',
                    #metavar='',
                    help="processed directory (default: {:s})".format(def_processed_path))
parser.add_argument('-d',
                    '--database',
                    default=def_database_path,
                    dest='database',
                    #metavar='',
                    help="database directory (default: {:s})".format(def_database_path))
args = parser.parse_args()
server_path = os.path.abspath(args.server)
processed_path = os.path.abspath(args.processed)
database_path = os.path.abspath(args.database)

# Set an inclusive time range of analysis for a specific float
# (by default, deployment to present...adjust here or there)
filterDate = utils.deploy2present()

# Boolean set to true in order to delete every processed data and redo everything
redo = False

# Figures commented by default (take heaps of memory)
write_png = False
write_html = False
write_sac = True
write_mseed = True
write_mhpsd = True

# Use WebGL implementation of graph to
# increase speed, improve interactivity, and the ability to plot even more data
optimized_html = False

# Integrate the plotly library into every html file
# If false user must have internet connection to access graph
# but the file size will be reduced considerably
local_html = True

# Dictionary to write last-cycle vital data to output files
lastcycle = {}


def sort_events(eventA,eventB) :
    if eventA.corrected_starttime and eventB.corrected_starttime :
        return eventA.corrected_starttime - eventB.corrected_starttime
    elif eventA.corrected_starttime and not eventB.corrected_starttime :
        return eventA.corrected_starttime - eventB.uncorrected_starttime
    elif not eventA.corrected_starttime and eventB.corrected_starttime :
        return eventA.uncorrected_starttime - eventB.corrected_starttime
    else :
        return eventA.uncorrected_starttime - eventB.uncorrected_starttime


def sort_mfloats(a,b):
    buoy_nbA = a.split("-")[-1]
    buoy_nbB = b.split("-")[-1]
    nbA = int(buoy_nbA,10)
    nbB = int(buoy_nbB,10)
    return nbA - nbB

def main():
    # Set working directory in "scripts"
    os.chdir(scripts_path)

    # Create processed directory if it doesn't exist
    if not os.path.exists(processed_path):
        os.mkdir(processed_path)

    # Search MERMAID floats
    vitfile_path = os.path.join(server_path, "[0-9]*.*-*-*[0-9].vit")
    mfloats = [p.split("/")[-1][:-4] for p in glob.glob(vitfile_path)]

    # Create database directory if it doesn't exist
    if not os.path.exists(database_path):
        os.mkdir(database_path)

    # Update Database
    preprocess.database_update(database_path)

    # Sort *.vit path
    mfloats_sorted = sorted(mfloats, key=functools.cmp_to_key(sort_mfloats))

    # For each MERMAID float
    for mfloat in mfloats_sorted:
        print("Processing {:s} .LOG & .MER files...".format(mfloat))

        # Set the path for the float
        mfloat_path = os.path.join(processed_path, mfloat, "")

        # Get float number
        mfloat_nb = re.findall("(\d+)$", mfloat)[0]

        # Delete the directory if the redo flag is true
        if redo and os.path.exists(mfloat_path):
            shutil.rmtree(mfloat_path)

        # Create directory for the float
        if not os.path.exists(mfloat_path):
            os.mkdir(mfloat_path)

        # Remove existing files in the processed directory (the script may have been previously
        # executed, copied the files, then failed)
        for f in glob.glob(mfloat_path + "*.*"):
            os.remove(f)

        # Copy appropriate files in the directory and remove files outside of the time range
        files_to_copy = []

        # All files begin with buoy_nb followed by underscore
        # Add underscore avoids errors between similar buoy numbers (Ex: 01_* and 0101_*)
        files_to_copy += glob.glob(os.path.join(server_path, mfloat_nb +  "_*"))

        # Add .cmd, .out, and .vit files
        files_to_copy += glob.glob(os.path.join(server_path, mfloat + "*"))

        # Copy files
        for f in files_to_copy:
            shutil.copy(f, mfloat_path)

        # Concatenate all files for this float
        preprocess.concatenate_files(mfloat_path);

        # Decrypt all files for this float
        preprocess.decrypt_all(mfloat_path);

        # Determine the time range of analysis (generally; birth to death of a MERMAID)
        if mfloat in filterDate.keys():
            begin = filterDate[mfloat][0]
            end = filterDate[mfloat][1]
        else:
            begin = datetime.datetime(1000, 1, 1)
            end = datetime.datetime(3000, 1, 1)

        # Convert in cycle files
        preprocess.convert_in_cycle(mfloat_path,begin,end);

        # Really: collect all the .MER files (next we correlate their environments to .LOG files)
        print(" ...compiling a list of events from {:s} .MER files (GPS & seismic data)..." \
              .format(mfloat))
        mevents = events.Events(mfloat_path)
        # Build list of all S41 profiles recorded
        ms41s = sbe41.Profiles(mfloat_path)
        # Build list of all S61 profiles recorded
        ms61s = sbe61.Profiles(mfloat_path)

        # Collect all the .CYCLE files
        print(" ...matching those events to {:s} .LOG ('dive') files (GPS & dive metadata)..." \
              .format(mfloat))
        cycle_logs = cycles.get_cycles(mfloat_path, mevents, ms41s, ms61s)

        kml.generate(mfloat_path, mfloat, cycle_logs)


        # Verify dive logs are sorted as expected
        if cycle_logs!= sorted(cycle_logs, key=lambda x: x.start_date):
            raise ValueError('`cycle_logs` improperly sorted')

        for i, cycle_log in enumerate(cycle_logs):
            # Create the directory
            if not os.path.exists(cycle_log.processed_path):
                os.mkdir(cycle_log.processed_path)

            # Reformat and write .LOG in individual dive directory
            cycle_log.write_datetime_cycle()

            # Write .MER environment in individual directories
            cycle_log.write_mermaid_environment_files()

            # Write .S41 environment in individual directories
            cycle_log.write_s41_environment_file();

            # Write .S61 environment in individual directories
            cycle_log.write_s61_environment_file();

            # Generate dive plot
            cycle_log.write_cycle_html(csv_file,optimize=optimized_html,include_plotly=local_html)
            # <-- timestamps not corrected for clockdrift

            # The GPS list is None outside of requested begin/end dates, within
            # which it defaults to an empty list if it is truly empty
            if cycle_log.gps_list is None:
                continue

            # Validate that the GPS may be used to correct various MERMAID
            # timestamps, including diving/surfacing and event starttimes
            cycle_log.validate_gps(min_gps_fix, max_gps_time)

            # Apply clock corrections to the events associated with this
            # completed dive
            cycle_log.correct_clockdrifts()

            # Set output (.sac, .mseed) file names of the events associated with
            # this cycle using the adjusted and corrected event dates
            cycle_log.set_processed_file_names()

            # Interpolate station locations at various points in the dive
            cycle_log.compute_station_locations(mixed_layer_depth_m, preliminary_location_ok)

            # Format station-location metadata for ObsPy and attach to complete dive object
            cycle_log.set_events_obspy_trace_stats()

            # Write profiles html
            cycle_log.write_profile_html(optimize=optimized_html,include_plotly=local_html)

            # Write profiles data on CSV
            if csv_file :
                cycle_log.write_profile_csv();

            # Write requested output files
            if write_png:
                cycle_log.write_events_png()

            if write_html:
                cycle_log.write_events_html(optimize=optimized_html,include_plotly=local_html)

            if write_sac:
                cycle_log.write_events_sac()

            if write_mseed:
                cycle_log.write_events_mseed()

            if write_mhpsd:
                cycle_log.write_events_mhpsd(creation_datestr)

        # Verify events sublists are sorted as expected
        events_list = [event for cycle in cycle_logs for event in cycle.events]
        # Sort event lists by corrected starttime is exist => use uncorrected_starttime elsewhere
        if events_list != sorted(events_list, key=functools.cmp_to_key(sort_events)):
            raise ValueError('`cycle_logs[*].events` improperly sorted')
        # Plot vital data
        kml.generate(mfloat_path, mfloat, cycle_logs)
        vitals.plot_battery_voltage(mfloat_path, mfloat + ".vit", begin, end)
        vitals.plot_internal_pressure(mfloat_path, mfloat + ".vit", begin, end)
        vitals.plot_pressure_offset(mfloat_path, mfloat + ".vit", begin, end)
        if len(cycle_logs) > 1:
            vitals.plot_corrected_pressure_offset(mfloat_path, cycle_logs, begin, end)

        # NB, at this point, the total event lists associated with `dive_logs`
        # and `cycle_logs` may differ because the former collects all events
        # and the latter winnows that list to only include unique events (via
        # `dives.set_processed_file_names`, which removes redundant events from
        # individual `cycle_logs.events` lists); ergo, one may use the
        # existence of `event.station_loc` to determine what events in
        # `dive_logs` were actually retained in `cycle_logs` (see e.g.,
        # `events.write_traces_txt`)

        # Write csv and txt files containing all GPS fixes from .LOG and .MER
        gps.write_gps(cycle_logs, creation_datestr, processed_path, mfloat_path)

        # Write text file detailing event-station location interpolation parameters
        gps.write_gps_interpolation_txt(cycle_logs,creation_datestr, processed_path, mfloat_path)

        # Write text file detailing which SINGLE .LOG and .MER files define
        # (possibly incomplete) dives
        cycles.write_logs_txt(cycle_logs, creation_datestr,  processed_path, mfloat_path)

        # Write text file detailing .CYCLE files (init,complete dives, last dive)
        cycles.write_cycles_txt(cycle_logs, creation_datestr,  processed_path, mfloat_path,mfloat)

        # Write a text file relating all SAC and mSEED to their associated .LOG
        # and .MER files
        events.write_traces_txt(cycle_logs, creation_datestr, processed_path, mfloat_path)

        # Write a text file with our best-guess at the location of MERMAID at
        # the time of recording
        events.write_loc_txt(cycle_logs, creation_datestr, processed_path, mfloat_path)

        # Write mseed2sac and automaid metadata csv and text files
        events.write_obspy_trace_stats(cycle_logs, creation_datestr, processed_path, mfloat_path)

        # Write GeoCSV files
        geocsv_meta = geocsv.GeoCSV(cycle_logs, creation_datestr, mixed_layer_depth_m)
        geocsv_meta.write(os.path.join(processed_path, mfloat_path, 'geo.csv'))

        # Clean directories
        files_to_delete = list()
        files_to_delete += glob.glob(mfloat_path + "/" + mfloat_nb + "_*.MER")
        files_to_delete += glob.glob(mfloat_path + "/" + mfloat_nb + "_*.S41")
        files_to_delete += glob.glob(mfloat_path + "/" + mfloat_nb + "_*.S61")
        files_to_delete += glob.glob(mfloat_path + "/" + mfloat_nb + "_*.LOG")
        files_to_delete += glob.glob(mfloat_path + "/" + mfloat_nb + "_*.BIN")
        files_to_delete += glob.glob(mfloat_path + "/" + "*.CYCLE")
        for f in files_to_delete:
            os.remove(f)
        # Save the last complete dive of this float to later write output list
        # of external pressure measurements for the entire array
        if cycle_logs:
            lastcycle[mfloat] = cycle_logs[-1]

        # Remove lingering incomplete "IcCycle" folders, if completed
        mfloat_files = os.listdir(mfloat_path)
        incomplete_cycles = list(filter(lambda x: 'IcCycle' in x, mfloat_files))
        for incomplete_cycle in incomplete_cycles:
            complete_cycle = incomplete_cycle.replace('IcCycle', '')
            if os.path.exists(os.path.join(mfloat_path, complete_cycle)):
                shutil.rmtree(os.path.join(mfloat_path, incomplete_cycle))

    # Done looping through all dives for each float
    #______________________________________________________________________________________#

    # Print a text file of corrected external pressures measured on the final
    # dive, and warn if any are approaching the limit of 300 mbar (at which
    # point adjustment is required)
    vitals.write_corrected_pressure_offset(lastcycle, processed_path)

if __name__ == "__main__":
    main()
