# Part of automaid -- a Python package to process MERMAID files
# pymaid environment (Python v2.7)
#
# Original author: Sebastien Bonnieux
# Current maintainer: Joel D. Simon (JDS)
# Contact: jdsimon@alumni.princeton.edu | joeldsimon@gmail.com
# Last modified by JDS: 17-Dec-2020, Python 2.7.15, Darwin-18.7.0-x86_64-i386-64bit

import os
import glob
import re
import subprocess
from obspy import UTCDateTime
from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.core.trace import Stats
import plotly.graph_objs as graph
import plotly.offline as plotly
import matplotlib.pyplot as plt
import utils
import gps
import sys
import numpy as np
import setup

# Get current version number.
version = setup.get_version()

class Events:
    '''The Events (plural) class references a SINGLE .MER file, and all events that
     live within it, which may be associated with the environments of multiple
     other .MER files

     Multiple events (event binary blocks) may exist in a single Events.mer_name

    '''

    def __init__(self, base_path=None, mer_name=None):
        self.mer_name = mer_name
        self.base_path = base_path
        self.events = list()
        self.__version__ = version

        # If just a base path to (e.g., a server directory) is passed, load all
        # .MER files contained there; otherwise read a single input file
        if self.mer_name is None:
            mer_files = glob.glob(os.path.join(self.base_path, "*.MER"))
        else:
            mer_files = glob.glob(os.path.join(self.base_path, self.mer_name))

        for mer_file in mer_files:
            # This .MER file name
            mer_binary_name = mer_file.split("/")[-1]

            # The </EVENT> binary blocks contained in this .MER file
            with open(mer_file, "r") as f:
                content = f.read()
            events = content.split("</PARAMETERS>")[-1].split("<EVENT>")[1:]

            for event in events:
                # The header of this specific </EVENT> block (NOT the </ENVIRONMENT> of
                # the same .MER file, which may be unrelated (different time))
                mer_binary_header = event.split("<DATA>\x0A\x0D")[0]

                # The actual binary data contained in this </EVENT> block (the seismogram)
                mer_binary_binary = event.split("<DATA>\x0A\x0D")[1].split("\x0A\x0D\x09</DATA>")[0]

                # N.B:
                # "\x0A" is "\n": True
                # "\x0D" is "\r": True
                # "\x09" is "\t": True
                # https://docs.python.org/2/reference/lexical_analysis.html#string-and-bytes-literals
                # I don't know why Seb choose to represent the separators as
                # hex, I believe a valid split would be "...split('\n\r\t</DATA>')..."

                # The double split above is not foolproof; if the final data
                # block in the .MER file ends without </DATA> (i.e., the file
                # was not completely transmitted), the object 'binary' will just
                # return everything to the end of the file -- verify that the we
                # actually have the expected number of bytes (apparently len()
                # returns the byte-length of a string, though I am not super
                # happy with this solution because I would prefer to know the
                # specific encoding used for event binary...)
                actual_binary_length = len(mer_binary_binary)
                bytes_per_sample = int(re.search('BYTES_PER_SAMPLE=(\d+)', mer_binary_header).group(1))
                num_samples = int(re.search('LENGTH=(\d+)', mer_binary_header).group(1))
                expected_binary_length = bytes_per_sample * num_samples

                if actual_binary_length == expected_binary_length:
                    self.events.append(Event(mer_binary_name, mer_binary_header, mer_binary_binary))

            # Sort by date the list of events contained in this .MER file
            self.events.sort(key=lambda x: x.date)


    def __repr__(self):
        return "Events('{}', '{}')".format(self.base_path, self.mer_name)


    def get_events_between(self, begin, end):
        catched_events = list()
        for event in self.events:
            if begin < event.date < end:
                catched_events.append(event)
        return catched_events


class Event:
    '''The Event (singular) class references TWO .MER files, which may be the same,
    through Event.mer_binary_name (.MER file containing the </EVENT> binary
    data), and Event.mer_environment_name (.MER file containing the
    </ENVIRONMENT> metadata (e.g., GPS, clock drift, sampling freq. etc.)
    associated with that event)

    Only a SINGLE event (event binary block) is referenced by
    Event.mer_binary_name and Event.mer_environment_name

    '''

    def __init__(self, mer_binary_name=None, mer_binary_header=None, mer_binary_binary=None):
        self.mer_binary_name = mer_binary_name
        self.mer_binary_header = mer_binary_header
        self.mer_binary_binary = mer_binary_binary
        self.__version__ = version

        # Defaults
        self.mer_environment_name = None
        self.mer_environment = None

        self.data = None
        self.measured_fs = None
        self.decimated_fs = None
        self.trig = None
        self.depth = None
        self.temperature = None
        self.criterion = None
        self.snr = None
        self.scales = None

        self.date = None
        self.station_loc = None
        self.clockdrift_correction = None
        self.obspy_trace_stats = None

        self.is_requested = False

        self.scales = re.findall(" STAGES=(-?\d+)", self.mer_binary_header)[0]
        catch_trig = re.findall(" TRIG=(\d+)", self.mer_binary_header)
        if len(catch_trig) > 0:
            # Event detected with STA/LTA algorithm
            self.trig = int(catch_trig[0])
            date = re.findall(" DATE=(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{6})", mer_binary_header, re.DOTALL)
            self.date = UTCDateTime.strptime(date[0], "%Y-%m-%dT%H:%M:%S.%f")
            self.depth = int(re.findall(" PRESSURE=(-?\d+)", self.mer_binary_header)[0])
            self.temperature = int(re.findall(" TEMPERATURE=(-?\d+)", self.mer_binary_header)[0])
            self.criterion = float(re.findall(" CRITERION=(\d+\.\d+)", self.mer_binary_header)[0])
            self.snr = float(re.findall(" SNR=(\d+\.\d+)", self.mer_binary_header)[0])
        else:
            # Event requested by user
            self.is_requested = True
            date = re.findall(" DATE=(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})", mer_binary_header, re.DOTALL)
            self.date = UTCDateTime.strptime(date[0], "%Y-%m-%dT%H:%M:%S")

    def __repr__(self):
        # Hacked repr dunder because I can't print binary...
        if self.mer_binary_binary:
            bin_str = '<int32 binary>'
        else:
            bin_str = self.mer_binary_binary

        return "Event('{}', '{}', {})".format(self.mer_binary_name, self.mer_binary_header, bin_str)

    def set_environment(self, mer_environment_name, mer_environment):
        self.mer_environment_name = mer_environment_name
        self.mer_environment = mer_environment

    def find_measured_sampling_frequency(self):
        # Get the frequency recorded in the .MER environment header
        fs_catch = re.findall("TRUE_SAMPLE_FREQ FS_Hz=(\d+\.\d+)", self.mer_environment)
        self.measured_fs = float(fs_catch[0])
        #self.measured_fs = 40

        # Divide frequency by number of scales
        int_scl = int(self.scales)
        if int_scl >= 0:
            self.decimated_fs = self.measured_fs / (2. ** (6 - int_scl))
        else:
            # This is raw data sampled at 40Hz
            self.decimated_fs = self.measured_fs

    def correct_date(self):
        # Calculate the date of the first sample
        if self.is_requested:
            # For a requested event
            rec_file_date = re.findall("FNAME=(\d{4}-\d{2}-\d{2}T\d{2}_\d{2}_\d{2})", self.mer_binary_header)
            rec_file_date = UTCDateTime.strptime(rec_file_date[0], "%Y-%m-%dT%H_%M_%S")

            rec_file_ms = re.findall("FNAME=\d{4}-\d{2}-\d{2}T\d{2}_\d{2}_\d{2}\.?(\d{6}?)", self.mer_binary_header)
            if len(rec_file_ms) > 0:
                rec_file_date += float("0." + rec_file_ms[0])

            sample_offset = re.findall("SMP_OFFSET=(\d+)", self.mer_binary_header)
            sample_offset = float(sample_offset[0])
            self.date = rec_file_date + sample_offset/self.measured_fs
        else:
            # For a detected event
            # The recorded date is the STA/LTA trigger date, subtract the time before the trigger.
            self.date = self.date - float(self.trig) / self.decimated_fs

    def correct_clockdrift(self, gps_descent, gps_ascent):
        # Correct the clock drift of the Mermaid board with GPS measurement
        pct = (self.date - gps_descent.date) / (gps_ascent.date - gps_descent.date)
        self.clockdrift_correction = gps_ascent.clockdrift * pct
        # Apply correction
        self.date = self.date + self.clockdrift_correction

    def compute_station_location(self, drift_begin_gps, drift_end_gps):
        '''Fills attribute self.station_loc, the interpolated location of MERMAID when
        it recorded an event

        '''
        self.station_loc = gps.linear_interpolation([drift_begin_gps, drift_end_gps], self.date)

    def invert_transform(self, bin_path=os.path.join(os.environ["AUTOMAID"], "scripts", "bin")):
        # If scales == -1 this is a raw signal, just convert binary data to np array of int32
        if self.scales == "-1":
            self.data = np.frombuffer(self.mer_binary_binary, np.int32)
            return

        # Get additional information on flavor of invert wavelet transform
        normalized = re.findall(" NORMALIZED=(\d+)", self.mer_environment)[0]
        edge_correction = re.findall(" EDGES_CORRECTION=(\d+)", self.mer_environment)[0]

        # Change to binary directory because these scripts can fail with full paths
        start_dir = os.getcwd();
        os.chdir(bin_path)

        # The following scripts READ wavelet coefficients (what MERMAID
        # generally sends) from a file named "wtcoeffs" and WRITE the inverted
        # data to a file name, e.g., "wtcoeffs.icdf24_5"
        wtcoeffs_data_file_name = "wtcoeffs"
        inverted_data_file_name = "wtcoeffs.icdf24_" + self.scales

        # Delete any previously-inverted data just to be absolutely sure we are
        # working with this events' data only (an interruption before the second
        # call to delete these files could result in their persistence)
        if os.path.exists(wtcoeffs_data_file_name):
            os.remove(wtcoeffs_data_file_name)

        if os.path.exists(inverted_data_file_name):
            os.remove(inverted_data_file_name)

        # Write cdf24 data
        with open(wtcoeffs_data_file_name, 'w') as f:
            f.write(self.mer_binary_binary)

        # The inverse wavelet transform C code (icdf24_v103(ec)_test) is called
        # below in a subprocess and its output is verified; determine if edge
        # correction needs to be accounted for
        icdf24_arg_list = list()
        if edge_correction == "1":
            icdf24_arg_list.append("icdf24_v103ec_test")

        else:
            icdf24_arg_list.append("icdf24_v103_test")


        # Extend the icdf24_v103(ec)_test argument list with other data values
        icdf24_arg_list.extend([self.scales, normalized, wtcoeffs_data_file_name])

        # Perform inverse wavelet transform
        stdout = subprocess.check_output(icdf24_arg_list)

        # Ensure the inverse wavelet transform worked as expected, meaning that
        # it generated an output file of int32 data
        if not os.path.exists(inverted_data_file_name):
            cmd  = ' '.join(map(str, icdf24_arg_list)) # prints python list as comma-separated string
            err_mess = "\nFailed: inverse wavelet transformation\n"
            err_mess += "In directory: {:s}\n".format(bin_path)
            err_mess += "Attempted command: {:s}\n".format(cmd)
            err_mess += "Using: event around {:s} in {:s}\n\n".format(self.date, self.mer_binary_name)
            err_mess += "Command printout:\n'{:s}'".format(stdout)

            # This output message is more helpful than the program crashing on
            # the next line
            sys.exit(err_mess)

        # Read the inverted data
        self.data = np.fromfile(inverted_data_file_name, np.int32)

        # Delete the files of coefficient and inverted data, otherwise a latter
        # .MER with an incomplete binary event block can come along and use the
        # same data
        os.remove(wtcoeffs_data_file_name)
        os.remove(inverted_data_file_name)

        # Return to start directory.
        os.chdir(start_dir)

    def get_export_file_name(self):
        export_file_name = UTCDateTime.strftime(UTCDateTime(self.date), "%Y%m%dT%H%M%S") + "." + self.mer_binary_name
        if not self.trig:
            export_file_name = export_file_name + ".REQ"
        else:
            export_file_name = export_file_name + ".DET"
        if self.scales == "-1":
            export_file_name = export_file_name + ".RAW"
        else:
            export_file_name = export_file_name + ".WLT" + self.scales
        return export_file_name

    def __get_figure_title(self):
        title = "" + self.date.isoformat() \
                + "     Fs = " + str(self.decimated_fs) + "Hz\n" \
                + "     Depth: " + str(self.depth) + " m" \
                + "     Temperature: " + str(self.temperature) + " degC" \
                + "     Criterion = " + str(self.criterion) \
                + "     SNR = " + str(self.snr)
        return title

    def plotly(self, export_path, force_redo=False):
        # Check if file exist
        export_path_html = export_path + self.get_export_file_name() + ".html"
        if not force_redo and os.path.exists(export_path_html):
            return

        # Add acoustic values to the graph
        data_line = graph.Scatter(x=utils.get_date_array(self.date, len(self.data), 1./self.decimated_fs),
                                  y=self.data,
                                  name="counts",
                                  line=dict(color='blue',
                                            width=2),
                                  mode='lines')

        data = [data_line]

        layout = graph.Layout(title=self.__get_figure_title(),
                              xaxis=dict(title='Coordinated Universal Time (UTC)', titlefont=dict(size=18)),
                              yaxis=dict(title='Counts', titlefont=dict(size=18)),
                              hovermode='closest'
                              )

        plotly.plot({'data': data, 'layout': layout},
                    filename=export_path_html,
                    auto_open=False)

    def plot_png(self, export_path, force_redo=False):
        # Check if file exist
        export_path_png = export_path + self.get_export_file_name() + ".png"
        if not force_redo and os.path.exists(export_path_png):
            return

        pascals = [utils.counts2pascal(d) for d in self.data]

        # Plot frequency image
        plt.figure(figsize=(9, 4))
        plt.title(self.__get_figure_title(), fontsize=12)
        plt.plot(utils.get_time_array(len(pascals), 1./self.decimated_fs),
                 pascals,
                 color='b')
        plt.xlabel("Time (s)", fontsize=12)
        plt.ylabel("Pascal", fontsize=12)
        plt.tight_layout()
        plt.grid()
        plt.savefig(export_path_png)
        plt.clf()
        plt.close()

    def attach_obspy_trace_stats(self, kstnm, kinst, force_without_loc=False):
        '''Attaches attribute: obspy_trace_stats, an obspy.core.trace.Stats instance.

        obspy_trace_stats holds metadata common to both miniSEED and SAC formats.
        obspy_trace_stats.sac holds extra metadata only found in the SAC format.

        Floats are NOT converted to np.float32() in either case.

        NB: the SAC header value shown to the world (e.g., "sac.delta"), and the private SAC header
        written to disk (e.g., "sac._hf[0]"), differ in type.  The relevant float header values that
        actually get written to disk with sac.write are stored in the private "._hf" attribute,
        which is not generated with initialization of the raw Stats() container. Therefore, if
        printing those values to, e.g. a text file, ensure the relevant F (float) fields are cast to
        np.float32 first.

        For example:
        >> from obspy.core.trace import Trace
        >> from obspy.io.sac.sactrace import SACTrace
        >> trace = Trace()
        >> sac = SACTrace.from_obspy_trace(trace)  <-- this gets called by sac.write (within stream.write)
        >> sac.delta = 1/20
        >> isinstance(sac.delta, float)            <-- True: this is the public attr shown to the world
        >> isinstance(sac.delta, np.float32)       <-- False
        >> isinstance(sac._hf[0], float)           <-- False
        >> isinstance(sac._hf[0], np.float32)      <-- True: this is the private attr written to disk

        For more detail see: http://www.adc1.iris.edu/files/sac-manual/manual/file_format.html

        Update function "write_metadata" if the fields in this method are changed.

        '''

        # Fill metadata common to SAC and miniSEED formats
        stats = Stats()
        stats.network = utils.network()
        stats.station = kstnm
        stats.location = "00"
        stats.channel = utils.band_code(self.decimated_fs) + "DH"  # SEED manual Appendix A
        stats.starttime = self.date
        stats.sampling_rate = self.decimated_fs
        stats.npts = len(self.data)

        # Extra metadata only written to SAC files (stel, cmpaz, and cmpinc are included here for
        # writing to mseed2sac_metadata.csv, but they are left unfilled)
        keys = ['stla',
                'stlo',
                'stel',
                'stdp',
                'scale',
                'cmpaz',
                'cmpinc',
                'user0',
                'user1',
                'user2',
                'user3',
                'kinst',
                'kuser0',
                'kuser1']
        def_float = -12345.

        # Default SAC header (we may not will not fill all of these keys)
        stats.sac = dict.fromkeys(keys, def_float)

        # Fill station-location header fields.
        if not force_without_loc:
            stats.sac["stla"] = self.station_loc.latitude;
            stats.sac["stlo"] = self.station_loc.longitude;

        # Elevation is 0 (our reference is truly sea level)
        stats.sac["stel"] = 0

        # Add scaling factor to convert digital counts to Pa
        stats.sac["scale"] = utils.sac_scale()

        # Add dip (CMPINC; "component incidence") in SAC dip convention, using as guide:
        # https://github.com/iris-edu/mseed2sac/blob/master/doc/mseed2sac.md
        #
        # SAC dip convention: "degrees down from vertical up/outward",
        # i.e., BHN, BHE = 90, BHZ = 0
        #
        # SEED dip convection: "degrees down from horizontal"
        # i.e., BHN, BHE = 0, BHZ = -90
        stats.sac["cmpinc"] = 0 # SAC dip

        # Add azimuth: horizontal projection of component vector measured clockwise from north
        # It is 0 for vertical components. Theoretically, BHN, BHZ = 90, BHE = 90
        stats.sac["cmpaz"] = 0

        # NB: I checked how IRIS serves up hydrophone data (in MATLAB):
        # >> s = irisFetch.Stations('channel', '*', '*', '*', '?DH')
        #
        # For all 3233 channels from 2147 stations that were returned:
        # dip = -90, 0, or 90
        # azimuth = 0 or 360
        #
        # For dip = -90, I assume that is the SEED dip convention
        # For dip = +90, I do not know; I thought perhaps it might be some(thing like a?)
        # right-hand-rule convention, but not all +90 dips are associated with 360 azimuth

        # REQ events do not record their depth at the time of acquisition, and because the onboard
        # detection algorithm was not triggered there are no trigger parameters to report
        if not self.is_requested:
            stats.sac["stdp"] = self.depth # meters (from external pressure sensor; down is positive)
            stats.sac["user0"] = self.snr
            stats.sac["user1"] = self.criterion
            stats.sac["user2"] = self.trig # sample index

        # Clock drift is computed for both DET and REQ, unless prevented by GPS error (computation
        # not determined by DET or REQ status)
        stats.sac["user3"] = self.clockdrift_correction if self.clockdrift_correction else def_float # seconds

        # Generic instrument (e.g., '452.020')
        stats.sac['kinst'] = kinst

        # automaid version number
        stats.sac["kuser0"] = self.__version__

        # String describing detection/request status, and number of wavelet scales transmitted
        # (e.g., 'DET.WLT5')
        reqdet_scales = self.get_export_file_name().split('.')[-2:]
        stats.sac['kuser1'] = '.'.join(reqdet_scales)

        self.obspy_trace_stats = stats

    def to_mseed(self, export_path, kstnm, kinst, force_without_loc=False, force_redo=False):
        # NB, writes mseed2sac writes, e.g., "MH.P0025..BDH.D.2018.259.211355.SAC", where "D" is the
        # quality indicator, "D -- The state of quality control of the data is indeterminate" (SEED
        # v2.4 manual pg. 108)

        # Check if the station location has been calculated
        if self.station_loc is None and not force_without_loc:
            #print self.get_export_file_name() + ": Skip mseed generation, wait the next ascent to compute location"
            return

        # Format the metadata into miniSEED and SAC header formats
        if not self.obspy_trace_stats:
            self.attach_obspy_trace_stats(kstnm, kinst, force_without_loc)

        # Check if file exist
        export_path_msd = export_path + self.get_export_file_name() + ".mseed"
        if not force_redo and os.path.exists(export_path_msd):
            return

        # Get stream object
        stream = self.get_stream(export_path, kstnm, kinst, force_without_loc)

        # Save stream object
        stream.write(export_path_msd, format='MSEED')

    def to_sac(self, export_path, kstnm, kinst, force_without_loc=False, force_redo=False):
        # Check if the station location has been calculated
        if self.station_loc is None and not force_without_loc:
            #print self.get_export_file_name() + ": Skip sac generation, wait the next ascent to compute location"
            return

        # Format the metadata into miniSEED and SAC header formats
        if not self.obspy_trace_stats:
            self.attach_obspy_trace_stats(kstnm, kinst, force_without_loc)

        # Check if file exist
        export_path_sac = export_path + self.get_export_file_name() + ".sac"
        if not force_redo and os.path.exists(export_path_sac):
            return

        # Get stream object
        stream = self.get_stream(export_path, kstnm, kinst, force_without_loc)

        # Save stream object
        stream.write(export_path_sac, format='SAC')

    def get_stream(self, export_path, kstnm, kinst, force_without_loc=False):
        # Check if an interpolated station location exists
        if self.station_loc is None and not force_without_loc:
            return


        # Save data into a Stream object
        trace = Trace()
        trace.stats = self.obspy_trace_stats
        trace.data = self.data

        stream = Stream(traces=[trace])

        return stream


def write_traces_txt(mdives, processed_path, mfloat_path):
    event_dive_tup = ((event, dive) for dive in mdives for event in dive.events if event.station_loc)

    traces_file = os.path.join(processed_path, mfloat_path, "traces.txt")
    fmt_spec = '{:>40s}    {:>15s}    {:>15s}    {:>15s}    {:>15s}    {:>15s}    {:>15s}    {:>15s}\n'

    version_line = "automaid {} ({})\n\n".format(setup.get_version(), setup.get_url())
    header_line = "                               file_name            bin_mer      prev_dive_log  prev_dive_env_mer      this_dive_log  this_dive_env_mer      next_dive_log  next_dive_env_mer\n".format()

    with open(traces_file, "w+") as f:
        f.write(version_line)
        f.write(header_line)

        for e, d in sorted(event_dive_tup, key=lambda x: x[0].date):
            f.write(fmt_spec.format(e.get_export_file_name(),
                                    e.mer_binary_name,
                                    d.prev_dive_log_name,
                                    d.prev_dive_mer_environment_name,
                                    d.log_name,
                                    d.mer_environment_name,
                                    d.next_dive_log_name,
                                    d.next_dive_mer_environment_name))


def write_loc_txt(mdives, processed_path, mfloat_path):
    '''Writes interpolated station locations at the time of event recording for all events for each
    individual float

    '''

    event_list = [event for dive in mdives for event in dive.events if event.station_loc]

    loc_file = os.path.join(processed_path, mfloat_path, "loc.txt")
    fmt_spec = "{:>40s}    {:>10.6f}    {:>11.6f}    {:>6.0f}\n"

    version_line = "automaid {} ({})\n\n".format(setup.get_version(), setup.get_url())
    header_line = "                               file_name   interp_STLA    interp_STLO      STDP\n"

    with open(loc_file, "w+") as f:
        f.write(version_line)
        f.write(header_line)

        for e in sorted(event_list, key=lambda x: x.date):
            f.write(fmt_spec.format(e.get_export_file_name(),
                                    np.float32(e.obspy_trace_stats.sac["stla"]),
                                    np.float32(e.obspy_trace_stats.sac["stlo"]),
                                    np.float32(e.obspy_trace_stats.sac["stdp"])))

def write_metadata(mdives, processed_path, mfloat_path):
    '''Write mseed2sac metadata and automaid metadata files.

    Update this function if the fields in method "attach_obspy_trace_stats" are changed.

    In total four files are written:

    mseed2sac_metadata.csv (actually used by mseed2sac)
    mseed2sac_metadata.txt (same info; more human-readable)

    automaid_metadata.csv (ALL and ONLY SAC info defined in automaid)
    automaid_metadata.txt (same info; more human-readable)

    msee2sac_metadata.csv/txt:

        Usage: mseed2sac -m mseed2sac_metadata.csv *mseed

        From: https://github.com/iris-edu/mseed2sac/blob/master/doc/mseed2sac.md

        (01) Network (KNETWK)
        (02) Station (KSTNM)
        (03) Location (KHOLE)
        (04) Channel (KCMPNM)
        (05) Latitude (STLA)
        (06) Longitude (STLO)
        (07) Elevation (STEL), in meters [not currently used by SAC]
        (08) Depth (STDP), in meters [not currently used by SAC]
        (09) Component Azimuth (CMPAZ), degrees clockwise from north
        (10) Component Incident Angle (CMPINC), degrees from vertical
        (11) Instrument Name (KINST), up to 8 characters
        (12) Scale Factor (SCALE)
        (13) Scale Frequency, unused
        (14) Scale Units, unused
        (15) Sampling rate, unused
        (16) Start time, used for matching
        (17) End time, used for matching


    automaid_metadata.csv/txt:

        Prints ALL and ONLY the non-default SAC headers filled by automaid:

        (01) file name (from automaid; not a SAC header field)
        (02) KNETWK
        (03) KSTNM
        (04) KCMPNM
        (05) STLA
        (06) STLO
        (07) STEL
        (08) STDP
        (09) CMPAZ
        (10) CMPINC
        (11) KINST
        (12) SCALE
        (13) USER0 (SNR)
        (14) USER1 (criterion)
        (15) USER2 (trig)
        (16) USER3 (clockdrift correction)
        (17) KUSER0 (automaid version)
        (18) KUSER1 (REQ or DET and scales)
        (19) samplerate (not a SAC header field)
        (20) start (not a SAC header field)
        (21) end (not a SAC header field)


    '''

    ## NB, concerning filename abbreviations:
    ## m2s_* == mseed2sac
    ## atm_* == automaid*

    # Version line is the same for both
    version_line = "#automaid {} ({})\n\n".format(setup.get_version(), setup.get_url())

    # Generate header lines for all four files: generate .csv by replacing
    # spaces with commas in text format
    m2s_header_line_txt = "#net    sta   loc   chan           lat            lon      elev     depth   azimuth    SACdip  instrument     scale  scalefreq scaleunits samplerate                  start                    end\n"
    m2s_header_line_csv = ','.join(m2s_header_line_txt.split())  + '\n'

    atm_header_line_txt = "                               filename KNETWK    KSTNM KCMPNM          STLA           STLO STEL      STDP CMPAZ CMPINC      KINST     SCALE            USER0            USER1     USER2            USER3      KUSER0      KUSER1 samplerate                  start                    end\n"
    atm_header_line_csv = '#' + ','.join(atm_header_line_txt.split()) + '\n'
    atm_header_line_txt = '#' + atm_header_line_txt # add pount after comma substitution

    # Field specifiers for mseed2sac_metadata.csv and mseed2sac_metadata.txt
    m2s_fmt = ['{:>2s}',    # Network (KNETWK)
               '{:>5s}',    # Station (KSTNM)
               '{:>2s}',    # Location (KHOLE)
               '{:>3s}',    # Channel (KCMPNM)
               '{:>10.6f}', # Latitude (STLA)
               '{:>11.6f}', # Longitude (STLO)
               '{:>6.0f}',  # Elevation (STEL), in meters [not currently used by SAC]
               '{:>6.0f}',  # Depth (STDP), in meters [not currently used by SAC]
               '{:>6.0f}',  # Component Azimuth (CMPAZ), degrees clockwise from north
               '{:>6.0f}',  # Component Incident Angle (CMPINC), degrees from vertical
               '{:>8s}',    # Instrument Name (KINST), up to 8 characters
               '{:6.0f}',   # Scale Factor (SCALE)
               '{:>7.1f}',  # Scale Frequency, unused
               '{:>7s}',    # Scale Units, unused
               '{:>7.0f}',  # Sampling rate, unused
               '{:>19s}',   # Start time, used for matching
               '{:>19s}\n'] # End time, used for matching

    # Add four spaces between each field to format the text file
    m2s_fmt_txt  = '    '.join(m2s_fmt)

    # Add comma between each field and remove field width (non-decimal) to format the csv
    m2s_fmt_csv  = ','.join(m2s_fmt)
    m2s_fmt_csv  = re.sub(':>\d*', ':', m2s_fmt_csv)

    # Field specifiers for automaid_metadata.csv and automaid_metadata.txt format
    atm_fmt = ['{:>40s}',   # file name (from automaid; not a SAC header field)
               '{:>3s}',    # KNETWK
               '{:>5s}',    # KSTNM
               '{:>3s}',    # KCMPNM
               '{:>10.6F}', # STLA
               '{:>11.6f}', # STLO
               '{:>1.0F}',  # STEL
               '{:>6.0f}',  # STDP
               '{:>2.0F}',  # CMPAZ
               '{:>3.0f}',  # CMPINC
               '{:s}',      # KINST
               '{:>.0f}',   # SCALE
               '{:>13.6f}', # USER0 (detection SNR)
               '{:>13.6f}', # USER1 (detecion criterion)
               '{:>6.0f}',  # USER2 (detecion trigger sample index)
               '{:>13.6f}', # USER3 (clockdrift correction)
               '{:>8s}',    # KUSER0 (automaid version)
               '{:>8s}',    # KUSER1 (REQ or DET and scales)
               '{:>7.0f}',  # samplerate (not a SAC header field)
               '{:>19s}',   # start (not a SAC header field)
               '{:>19s}\n'] # end (not a SAC header field)

    # Add four spaces between each field to format the text file
    atm_fmt_txt  = '    '.join(atm_fmt)

    # Add comma between each field and remove field width (non-decimal) to format the csv
    atm_fmt_csv  = ','.join(atm_fmt)
    atm_fmt_csv  = re.sub(':>\d*', ':', atm_fmt_csv)

    # The base path (the folder) is the same for all four files
    base_path = os.path.join(processed_path, mfloat_path)
    m2s_path =  os.path.join(base_path, 'mseed2sac_metadata')
    atm_path =  os.path.join(base_path, 'automaid_metadata')

    # These are mseed2sac_metadata values that do not differ(yet?) between MERMAIDs
    scalefreq = np.float32(1.)
    scaleunits = 'Pa'

    # Open all four files
    with open(m2s_path+".txt", "w+") as m2s_f_txt, \
         open(m2s_path+".csv", "w+") as m2s_f_csv, \
         open(atm_path+'.txt', "w+") as atm_f_txt, \
         open(atm_path+'.csv', "w+") as atm_f_csv:

        ## Write version line and header line to all four files

        m2s_f_csv.write(version_line)
        m2s_f_csv.write(m2s_header_line_csv)

        m2s_f_txt.write(version_line)
        m2s_f_txt.write(m2s_header_line_txt)

        atm_f_csv.write(version_line)
        atm_f_csv.write(atm_header_line_csv)

        atm_f_txt.write(version_line)
        atm_f_txt.write(atm_header_line_txt)

        # Loop over all events for which a station location was computed
        event_list = [event for dive in mdives for event in dive.events if event.station_loc]
        for e in sorted(event_list, key=lambda x: x.date):

            ## Collect metadata and convert to np.float32()

            # For mseed2sac_metadata.csv/txt
            net = e.obspy_trace_stats["network"]
            sta = e.obspy_trace_stats["station"]
            loc = e.obspy_trace_stats["location"]
            chan = e.obspy_trace_stats["channel"]
            lat = np.float32(e.obspy_trace_stats.sac["stla"])
            lon = np.float32(e.obspy_trace_stats.sac["stlo"])
            elev = np.float32(e.obspy_trace_stats.sac["stel"])
            depth = np.float32(e.obspy_trace_stats.sac["stdp"])
            azimuth = np.float32(e.obspy_trace_stats.sac["cmpaz"])
            SACdip = np.float32(e.obspy_trace_stats.sac["cmpinc"])
            instrument = e.obspy_trace_stats.sac["kinst"]
            scale = np.float32(e.obspy_trace_stats.sac["scale"])
            # scalefreq (defined above)
            # scaleunits (defined above)
            samplerate = np.float32(e.obspy_trace_stats["sampling_rate"])
            start = str(e.obspy_trace_stats["starttime"])[:19]
            end = str(e.obspy_trace_stats["endtime"])[:19]

            # Additional fields defined by automaid that are not in mseed2sac_metadata*
            filename = e.get_export_file_name()
            # KNETWK = net  (LHS are SAC names; RHS are their mseed2sac equivalents)
            # KSTNM = sta
            # KCMPNM = chan
            # STLA = lat
            # STLO = lon
            # ELEV = elev
            # STDP = depth
            # CMPAZ = azimuth
            # CMPINC = SACdip
            # KINST = instrument
            # SCALE = scale
            USER0 = np.float32(e.obspy_trace_stats.sac["user0"])
            USER1 = np.float32(e.obspy_trace_stats.sac["user1"])
            USER2 = np.float32(e.obspy_trace_stats.sac["user2"])
            USER3 = np.float32(e.obspy_trace_stats.sac["user3"])
            KUSER0 = e.obspy_trace_stats.sac["kuser0"]
            KUSER1 = e.obspy_trace_stats.sac["kuser1"]
            # samplerate (these last three already defined)
            # start
            # end

            ## Group into correct order

            # mseed2sac_metadata.csv fields
            m2s_meta = [net,
                        sta,
                        loc,
                        chan,
                        lat,
                        lon,
                        elev,
                        depth,
                        azimuth,
                        SACdip,
                        instrument,
                        scale,
                        scalefreq,
                        scaleunits,
                        samplerate,
                        start,
                        end]

            # automaid_metadata.csv fields, with SAC names commented
            atm_meta = [filename,
                        net,
                        sta,
                        chan,
                        lat,
                        lon,
                        elev,
                        depth,
                        azimuth,
                        SACdip,
                        instrument,
                        scale,
                        USER0,
                        USER1,
                        USER2,
                        USER3,
                        KUSER0,
                        KUSER1,
                        samplerate,
                        start,
                        end]

            ## Write to all files.

            m2s_f_txt.write(m2s_fmt_txt.format(*m2s_meta))
            m2s_f_csv.write(m2s_fmt_csv.format(*m2s_meta))

            atm_f_txt.write(atm_fmt_txt.format(*atm_meta))
            atm_f_csv.write(atm_fmt_csv.format(*atm_meta))
