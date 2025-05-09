# -*- Coding: utf-8 -*-
#
# Part of automaid -- a Python package to process MERMAID files
# pymaid environment (Python v3.10)
#
# Developer: Joel D. Simon <JDS>
# Developer: Frédéric rocca <FRO>
# Contact:  frederic.rocca@osean.fr
# Last modified by JDS: 21-Apr-2025
# Last modified by FRO: 09-Sep-2024
# Last tested: Python 3.10.13, 22.04.3-Ubuntu

import os
import re

from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth

import setup

# Get current version number.
version = setup.get_version()

class GPS:
    def __init__(self, date=None, latitude=None,
                 longitude=None,hdop=None,vdop=None,clockdrift=None, clockfreq=None,
                 source=None, rawstr_dict=None):
        self.__version__ = version
        self.date = date
        self.latitude = latitude
        self.longitude = longitude
        self.hdop = hdop
        self.vdop = vdop
        self.clockdrift = clockdrift
        self.clockfreq = clockfreq
        self.rawstr_dict = rawstr_dict
        self.source = source
        # The miniSEED convention of a time delay is of opposite sign of the
        # `clockdrift` (or the `clockdrift_correction`) of this program
        #
        # (+) clockdrift = MER time EARLY w.r.t GPS = (-) MER (record) time delay
        # (-) clockdrift = MER time LATE w.r.t GPS = (+) MER (record) time delay
        if self.clockdrift is not None:
            self.mseed_time_delay = -self.clockdrift
        else:
            self.mseed_time_delay = None
    def __len__(self):
        # To check if a single GPS instance passed into something that expects a list
        return 1


class GPS_interp(GPS):
    def __init__(self, date=None, latitude=None, longitude=None,
                 clockdrift=None, clockfreq=None, hdop=None, vdop=None,
                 interp_dict=None):

        super(GPS_interp, self).__init__(date=date, latitude=latitude,
                                         longitude=longitude,
                                         clockdrift=clockdrift,
                                         clockfreq=clockfreq, hdop=hdop,
                                         vdop=vdop)

        self.interp_dict = interp_dict

def linear_interpolation(gps_list, date):
    '''linear_interpolation(gps_list, date)

    Attempts to interpolate for the location of MERMAID at the requested date given an input list of
    GPS instances.  If a single GPS instance is give as input, interpolation cannot be performed and
    the input is returned (the interpolation is fixed to the input).

    Only two GPS instances from the input GPS list are retained for interpolation, and, if possible
    this method retains only the first two GPS instances that are separated by more than 10 minutes.
    If the two retained locations are within 20 m of one another interpolation is not performed and
    the interpolated location and date are fixed to an input (i.e., the requested interpolation date
    is not returned in the output GPS instance).

    The interpolation dictionary, ".interp_dict", attached to each instance attempts to explain the
    outcome of this method.

    '''

    # Defaults for the interpolation dictionary to be attached to this GPS instance
    i = None
    j = None

    # "input" -- difference between first and last locations retained (which may not even be actual
    # GPS fixes in the cases when the inputs are already an interpolated points)
    input_drift_dist_m = None
    input_drift_time = None
    input_drift_vel_ms = None
    input_lat_drift_dist_deg = None
    input_lat_drift_vel_degs = None
    input_lon_drift_dist_deg = None
    input_lon_drift_vel_degs = None

    # "interp" -- difference between nearest/reference GPS point and interpolation point
    interp_drift_dist_m = None
    interp_drift_time = None
    interp_drift_vel_ms = None
    interp_lat_drift_dist_deg = None
    interp_lat_drift_vel_degs = None
    interp_lon_drift_dist_deg = None
    interp_lon_drift_vel_degs = None

    # Return prematurely if the GPS list is not a list
    if len(gps_list) == 1:
        interp_lat = gps_list[0].latitude
        interp_lon = gps_list[0].longitude
        description = "interpolation not attempted (GPS list of length 1)"

        if date == gps_list[0].date:
            # Set time and drift distances to 0; leave velocities undefined
            input_drift_dist_m = 0.0
            input_drift_time = 0.0
            input_lat_drift_dist_deg = 0.0
            input_lon_drift_dist_deg = 0.0

            interp_drift_dist_m = 0.0
            interp_drift_time = 0.0
            interp_lat_drift_dist_deg = 0.0
            interp_lon_drift_dist_deg = 0.0
            description += "; interpolation not required (interpolation date is gps_list.date)"

        else:
            date = gps_list[0].date # overwrite: we did not interpolate at the requested date
            description += "; location and date fixed to input gps_list"

        interp_dict = locals()

        return GPS_interp(date=date,
                          latitude=interp_lat,
                          longitude=interp_lon,
                          clockdrift=None,
                          clockfreq=None,
                          hdop=None,
                          vdop=None,
                          interp_dict=interp_dict)

    # Return prematurely if the requested date is included in the GPS list
    for gps in gps_list:
        if date == gps.date:
            interp_lat = gps.latitude
            interp_lon = gps.longitude

            # Set time and drift distances to 0; leave velocities undefined
            input_drift_dist_m = 0.0
            input_drift_time = 0.0
            input_lat_drift_dist_deg = 0.0
            input_lon_drift_dist_deg = 0.0

            interp_drift_dist_m = 0.0
            interp_drift_time = 0.0
            interp_lat_drift_dist_deg = 0.0
            interp_lon_drift_dist_deg = 0.0
            description = "interpolation not required (interpolation date in gps_list)"

            interp_dict = locals()

            return GPS_interp(date=date,
                              latitude=interp_lat,
                              longitude=interp_lon,
                              clockdrift=None,
                              clockfreq=None,
                              hdop=None,
                              vdop=None,
                              interp_dict=interp_dict)

    # Otherwise, try to interpolate...

    # Ensure input list is sorted
    gps_list.sort(key=lambda x: x.date)

    # Identify the reference GPS points (gps_list[i]):
    # If date is before any gps fix compute drift from the two first gps fix
    if date < gps_list[0].date:
        # In this case: gps_list[i] is the FIRST GPS fix AFTER the interpolation date
        i = 0
        j = 1
        # Try to get a minimum time between two gps fix of 10 minutes
        while abs(gps_list[j].date - gps_list[i].date) < 10*60 and j < len(gps_list)-1:
            j += 1
        # Try to get a minimum distance between two gps fix of 20 meters
        while gps2dist_azimuth(gps_list[j].latitude, gps_list[j].longitude,
                               gps_list[i].latitude, gps_list[i].longitude)[0] < 20 and j < len(gps_list)-1:
            j += 1


    # If date is after any gps fix compute drift from the two last gps fix
    elif date > gps_list[-1].date:
        # In this case gps_list[i] is the LAST GPS fix BEFORE the interpolation date
        i = -1
        j = -2
        # Try to get a minimum time between two gps fix of 10 minutes
        while abs(gps_list[j].date - gps_list[i].date) < 10 * 60 and abs(j) < len(gps_list):
            j -= 1
        # Try to get a minimum distance between two gps fix of 20 meters
        while gps2dist_azimuth(gps_list[j].latitude, gps_list[j].longitude,
                               gps_list[i].latitude, gps_list[i].longitude)[0] < 20 and abs(j) < len(gps_list):
            j -= 1

    else:
        # If date is between two gps fix find the appropriate gps fix
        # In this case gps_list[i] is the LAST GPS fix BEFORE the interpolation date
        # In this case gps_list[j] is the FIRST GPS fix AFTER the interpolation date
        i = 0
        j = 1
        while not gps_list[i].date < date < gps_list[j].date and j < len(gps_list)-1:
            i += 1
            j += 1

    description= "interpolation attempted using multiple GPS points"
    # If the distance between the two GPS points retained is less than 20 m, don't interpolate just
    # return the one nearest in time to the requested date (don't simply fix a known location to the
    # requested interpolation date because it may happen that the input locations may be very near
    # to one another, but their dates may be very far from the requested date -- imagine if you gave
    # this algorithm two points, which differ by 1 meter and 1 second, and requested an interpolated
    # location date one after the last input date)
    if gps2dist_azimuth(gps_list[j].latitude, gps_list[j].longitude, gps_list[i].latitude, gps_list[i].longitude)[0] < 20:
        interp_lat = gps_list[i].latitude
        interp_lon = gps_list[i].longitude
        date = gps_list[i].date  # overwrite: we did not interpolate at the requested date
        description += "; retained points too close (spatially) for interpolation; location and date fixed to one of input gps_list"

    else:
        input_drift_dist_m = gps2dist_azimuth(gps_list[j].latitude, gps_list[j].longitude, \
                                              gps_list[i].latitude, gps_list[i].longitude)[0]
        input_drift_time = gps_list[j].date - gps_list[i].date
        input_drift_vel_ms = input_drift_dist_m / input_drift_time

        input_lat_drift_dist_deg = gps_list[j].latitude - gps_list[i].latitude
        input_lat_drift_vel_degs = input_lat_drift_dist_deg / input_drift_time

        # Do longitude arithmetic on [0:360] instead of [-180:180] system to
        # avoid wrapping issues when float crosses antimeridian
        # Ex.
        #    lon[i] = -174 (converted: 186)
        #    lon[j] = +178
        # Drifted west -8 degrees, not +352 if we do lon[j]-lon[i]
        # At conclusion of this algorithm, subtract 360 (if >180) to convert
        # back to [-180:+180]
        lonj = gps_list[j].longitude if gps_list[j].longitude >= 0 else gps_list[j].longitude + 360
        loni = gps_list[i].longitude if gps_list[i].longitude >= 0 else gps_list[i].longitude + 360

        # This is a bit of a cheat because it assumes an longitude lines are equally spaced on the
        # sphere, which they are not
        input_lon_drift_dist_deg = lonj - loni
        input_lon_drift_vel_degs = input_lon_drift_dist_deg / input_drift_time

        interp_drift_time = date - gps_list[i].date
        interp_lat_drift_vel_degs = input_lat_drift_vel_degs # they must equal
        interp_lat_drift_dist_deg = interp_lat_drift_vel_degs * interp_drift_time
        interp_lat = gps_list[i].latitude + interp_lat_drift_dist_deg

        interp_lon_drift_vel_degs = input_lon_drift_vel_degs # they must equal
        interp_lon_drift_dist_deg = interp_lon_drift_vel_degs * interp_drift_time
        interp_lon = loni + interp_lon_drift_dist_deg

        # Previously longitudes converted to [0:360] from [-180:180]
        if loni > 180:
            loni -= 360

        if lonj > 180:
            lonj -= 360

        if interp_lon > 180:
            interp_lon -= 360

        # This is also a bit of flub -- the interpolated drift distance computed here is using our
        # (ever so slightly) incorrect longitude, so when projected on a sphere we get a slightly
        # different distance than in our equal-box lat/lon projection; as such, the interpolated
        # drift velocity, which in reality must equal the drift velocity computed from the input,
        # will be slightly different
        # print("************************")
        # print("first {} source {}".format(gps_list[i].rawstr_dict,gps_list[i].source))
        # print("second {} source {}".format(gps_list[j].rawstr_dict,gps_list[j].source))
        # print("latitude 1 {} / latitude 2 {}".format(gps_list[i].latitude,gps_list[j].latitude))
        # print("longitude 1 {} / longitude 2 {}".format(gps_list[i].longitude,gps_list[j].longitude))
        # print("interp_lat_drift_dist_deg {} ".format(interp_lat_drift_dist_deg))
        # print("interp_drift_time {} ".format(interp_drift_time))
        # print("interp_lat_drift_vel_degs {} ".format(interp_lat_drift_vel_degs))
        #
        # print("input_lat_drift_vel_degs {} ".format(input_lat_drift_vel_degs))
        # print("input_lat_drift_dist_deg {} ".format(input_lat_drift_dist_deg))
        # print("input_drift_time {} ".format(input_drift_time))
        # print("input_drift_vel_ms {} ".format(input_drift_vel_ms))
        # print("input_drift_dist_m {}".format(input_drift_dist_m))
        # print("************************")
        interp_drift_dist_m = gps2dist_azimuth(interp_lat, interp_lon, gps_list[i].latitude, gps_list[i].longitude)[0]
        interp_drift_vel_ms = interp_drift_dist_m / interp_drift_time
        description += "; executed successfully"

    interp_dict = locals()

    # Nicety: >>> from pprint import pprint
    #         >>> pprint(interp_dict)
    return GPS_interp(date=date,
                      latitude=interp_lat,
                      longitude=interp_lon,
                      clockdrift=None,
                      clockfreq=None,
                      hdop=None,
                      vdop=None,
                      interp_dict=interp_dict)

def valid_clockfreq(GPS_object):
    '''Returns True if the clock frequency associated with a single GPS object is
    valid, between 3000000 and 4000000 Hz (MERMAID Manual Ref : 452.000.852
    Version 00, pg. 16), implying that MERMAID's onboard clock was properly
    synchronized.

    Note that, especially, the last(first) GPS fix taken before(after) the dive
    cannot display an anomalous clock frequency because that means the MERMAID
    clock had a synchronization error. Most usually, JDS has found, the
    clockdrift just keeps growing and it therefore seemingly would be acceptable
    to use the previous(next) valid GPS point to compute clockdrifts.

    For example, in this case it looks like the GPS points taken on March 18
    both failed to synchronize MERMAID's clock thus the clockdrift keep growing
    until the next GPS on the 24th.  This would imply it would be acceptable to
    correct for the clockdrift on the 18th using the clockdrifts on the 24th.

    Date                 Date Source      Loc Source       Clockdrift Clockfreq
    2019-03-01T14:59:55, 25_5C846B54.MER, 25_5C794585.LOG: +0.000062 ; 3686332
    2019-03-01T15:06:34, 25_5C846B54.MER, 25_5C794585.LOG: +0.000000 ; 3686332
    2019-03-10T02:14:28, 25_5C8F954C.MER, 25_5C846B61.LOG: -0.659393 ; 3686330
    2019-03-18T12:58:18, 25_5C96D514.MER, 25_5C8F9558.LOG: -1.045867 ; -1
    2019-03-18T13:04:25, 25_5C96D514.MER, 25_5C8F9558.LOG: -1.046020 ; -1
    2019-03-24T00:53:14, 25_5C96D514.MER, 25_5C8F9558.LOG: -1.709381 ; 3686333
    2019-03-24T00:55:16, 25_5CA1F464.MER, 25_5C96D552.LOG: -0.000030 ; 3686333
    2019-03-24T00:58:58, 25_5CA1F464.MER, 25_5C96D552.LOG: +0.000000 ; 3686333

    However, one can see here that the first synchronization on April 26th had
    an error, and then the following one 90 s later was normal. So we would
    think we could use the second GPS fix after surfacing to compute
    clockdrifts...but we cannot because it (appears in this case) that the clock
    partially reset, from -2.23 to -0.63 seconds.  So we have lost information
    on timing here.

    Date                 Date Source      Loc Source       Clockdrift Clockfreq
    2019-04-18T08:08:34, 25_5CC34C10.MER, 25_5CB83061.LOG: +0.000000 ; 3686333
    2019-04-18T08:18:51, 25_5CC34C10.MER, 25_5CB83061.LOG: +0.000000 ; 3686333
    2019-04-18T08:21:51, 25_5CC34C10.MER, 25_5CB83061.LOG: +0.000000 ; 3686332
    2019-04-26T18:20:18, 25_5CC34C10.MER, 25_5CB83061.LOG: -2.234710 ; 0
    2019-04-26T18:21:48, 25_5CCE6E90.MER, 25_5CC34C1D.LOG: -0.634887 ; 3686333
    2019-04-26T18:26:07, 25_5CC34C1D.LOG, 25_5CC34C1D.LOG: +0.000031 ; 3686332
    2019-04-26T18:28:25, 25_5CCE6E90.MER, 25_5CC34C1D.LOG: -0.000030 ; 3686332

    '''

    if 3000000 <= GPS_object.clockfreq <= 4000000:
        return True

    else:
        return False

def get_gps_from_mer_environment(mer_environment_name, mer_environment):

    '''

        Collect GPS fixes from MER environments within an inclusive datetime range

    '''
    gps_out = []

    # Mermaid environment can be empty
    if mer_environment is None:
        return gps_out

    # Get gps information in the mermaid environment
    gps_mer_list = mer_environment.split("</ENVIRONMENT>")[0].split("<GPSINFO")[1:]
    for gps_mer in gps_mer_list:
        rawstr_dict = {'fixdate': None, 'latitude': None, 'longitude': None, 'clockdrift': None}
        # .MER times are given simply as, e.g., "2020-10-20T02:36:55"
        fixdate = re.findall(" DATE=(\d+-\d+-\d+T\d+:\d+:\d+)", gps_mer)
        if len(fixdate) > 0:
            fixdate = fixdate[0]
            rawstr_dict['fixdate'] = fixdate
            fixdate = UTCDateTime(fixdate)
        else:
            fixdate = None

        # .MER latitudes are given as, e.g., "-2233.9800" (degrees decimal minutes) where the first 3
        # chars are the degrees (= S22deg33.9800mn) in .LOG parlance, with extra precision here
        latitude = re.findall(" LAT=([+,-])(\d{2})(\d+\.\d+)", gps_mer)
        if len(latitude) > 0:
            rawstr_dict['latitude'] = re.search("LAT=(.*) LON", gps_mer).group(1)
            latitude = latitude[0]
            if latitude[0] == "+":
                sign = 1
            elif latitude[0] == "-":
                sign = -1
            latitude = sign*(float(latitude[1]) + float(latitude[2])/60.)
        else:
            latitude = None

        # .MER longitudes are given as, e.g., "-14122.6800" (degrees decimal minutes) where the first
        # 4 chars are the degrees (= W141deg22.6800mn) in .LOG parlance, with an extra precision here
        longitude = re.findall(" LON=([+,-])(\d{3})(\d+\.\d+)", gps_mer)
        if len(longitude) > 0:
            rawstr_dict['longitude'] = re.search("LON=(.*) />", gps_mer).group(1)
            longitude = longitude[0]
            if longitude[0] == "+":
                sign = 1
            elif longitude[0] == "-":
                sign = -1
            longitude = sign*(float(longitude[1]) + float(longitude[2])/60.)
        else:
            longitude = None

        # .MER clockdrifts are given as, e.g.,
        # "<DRIFT YEAR=48 MONTH=7 DAY=4 HOUR=12 MIN=41 SEC=20 USEC=-563354 />"
        # which describe the drift using the sign convention of "drift = gps_time - mermaid_time"
        # (manual Ref: 452.000.852, pg. 32), NB: not all (any?) fields must exist (this is a
        # variable-length string); very often only USEC=*" will exist
        clockdrift = re.findall("<DRIFT( [^>]+) />", gps_mer)
        if len(clockdrift) > 0:
            rawstr_dict['clockdrift'] = re.search("<DRIFT (.*) />", gps_mer).group(1)
            clockdrift = clockdrift[0]
            _df = 0
            catch = re.findall(" USEC=(-?\d+)", clockdrift)
            if catch:
                _df += 10 ** (-6) * float(catch[0])
            catch = re.findall(" SEC=(-?\d+)", clockdrift)
            if catch:
                _df += float(catch[0])
            catch = re.findall(" MIN=(-?\d+)", clockdrift)
            if catch:
                _df += 60 * float(catch[0])
            catch = re.findall(" HOUR=(-?\d+)", clockdrift)
            if catch:
                _df += 60 * 60 * float(catch[0])
            catch = re.findall(" DAY=(-?\d+)", clockdrift)
            if catch:
                _df += 24 * 60 * 60 * float(catch[0])
            catch = re.findall(" MONTH=(-?\d+)", clockdrift)
            if catch:
                # An approximation of 30 days per month is sufficient this is just to see if there is something
                # wrong with the drift
                _df += 30 * 24 * 60 * 60 * float(catch[0])
            catch = re.findall(" YEAR=(-?\d+)", clockdrift)
            if catch:
                _df += 365 * 24 * 60 * 60 * float(catch[0])
            clockdrift = _df
        else:
            clockdrift = None

        clockfreq = re.findall("<CLOCK Hz=(-?\d+)", gps_mer)
        if len(clockfreq) > 0:
            clockfreq = clockfreq[0]
            clockfreq = int(clockfreq)
        else:
            clockfreq = None

        # Check if there is an error of clock synchronization
        # if clockfreq <= 0:
            # err_msg = "WARNING: Error with clock synchronization in file \"" + mer_environment_name + "\"" \
            #        + " at " + fixdate.isoformat() + ", clockfreq = " + str(clockfreq) + "Hz"
            # print err_msg

        # Add date to the list
        if fixdate is not None and latitude is not None and longitude is not None and clockdrift \
           is not None and clockfreq is not None:
                gps_out.append(GPS(date=fixdate,
                                             latitude=latitude,
                                             longitude=longitude,
                                             clockdrift=clockdrift,
                                             clockfreq=clockfreq,
                                             source=mer_environment_name,
                                             rawstr_dict=rawstr_dict))
        else:
            raise ValueError

    gps_out = sorted(gps_out, key=lambda x: x.date)
    return gps_out


def get_gps_from_log_content(log_name, log_content):
    '''

        Collect GPS fixes from LOG files within an inclusive datetime range

    '''
    gps_out = []
    gps_log_list = log_content.split("GPS fix...")[1:]
    for gps_log in gps_log_list:
        rawstr_dict = {'fixdate': None, 'latitude': None, 'longitude': None, 'clockdrift': None}
        # .LOG GPS times are given as integer UNIX Epoch times prepending the latitude longitude line
        # .LOG latitudes are given as, e.g., "S22deg33.978mn" (degrees and decimal minutes)
        latitude = re.findall("(\d+):\[SURF *, *\d+\]([S,N])(\d+)deg(\d+.\d+)mn", gps_log)
        if len(latitude) > 0:
            rawstr_dict['latitude'] = re.search("[S,N][0-9]+deg[0-9]+\.[0-9]+mn", gps_log).group(0)
            latitude = latitude[0]
            fixdate = latitude[0]
            rawstr_dict['fixdate'] = fixdate
            fixdate = UTCDateTime(int(fixdate))
            if latitude[1] == "N":
                sign = 1
            elif latitude[1] == "S":
                sign = -1
            latitude = sign*(float(latitude[2]) + float(latitude[3])/60.)
        else:
            fixdate = None
            latitude = None

        # .LOG latitudes are given as, e.g., "W141deg22.679mn" (degrees and decimal minutes)
        longitude = re.findall("([E,W])(\d+)deg(\d+.\d+)mn", gps_log)
        if len(longitude) > 0:
            rawstr_dict['longitude'] = re.search("[E,W][0-9]+deg[0-9]+\.[0-9]+mn", gps_log).group(0)
            longitude = longitude[0]
            if longitude[0] == "E":
                sign = 1
            elif longitude[0] == "W":
                sign = -1
            longitude = sign*(float(longitude[1]) + float(longitude[2])/60.)
        else:
            longitude = None

        hdop = re.findall("hdop (\d+.\d+)", gps_log)
        if len(hdop) > 0:
            hdop = hdop[0]
            hdop = float(hdop)
        else:
            hdop = None

        vdop = re.findall("vdop (\d+.\d+)", gps_log)
        if len(vdop) > 0:
            vdop = vdop[0]
            vdop = float(vdop)
        else:
            vdop = None

        clockdrift = re.findall("GPSACK:(.\d+),(.\d+),(.\d+),(.\d+),(.\d+),(.\d+),(.\d+)?;", gps_log)
        if len(clockdrift) > 0:
            clockdrift = clockdrift[0]
            rawstr_dict['clockdrift'] = clockdrift
            # YEAR + MONTH + DAY + HOUR + MIN + SEC + USEC
            clockdrift = 365 * 24 * 60 * 60 * float(clockdrift[0]) \
                + 30 * 24 * 60 * 60 * float(clockdrift[1]) \
                + 24 * 60 * 60 * float(clockdrift[2]) \
                + 60 * 60 * float(clockdrift[3]) \
                + 60 * float(clockdrift[4]) \
                + float(clockdrift[5]) \
                + 10 ** (-6) * float(clockdrift[6])
        else:
            clockdrift = None

        clockfreq = re.findall("GPSOFF:(-?\d+);", gps_log)
        if len(clockfreq) > 0:
            clockfreq = clockfreq[0]
            clockfreq = int(clockfreq)
        else:
            clockfreq = None

        if fixdate is not None and latitude is not None and longitude is not None:
            gps_out.append(GPS(date=fixdate,
                                             latitude=latitude,
                                             longitude=longitude,
                                             clockdrift=clockdrift,
                                             clockfreq=clockfreq,
                                             hdop=hdop,
                                             vdop=vdop,
                                             source=log_name,
                                             rawstr_dict=rawstr_dict))

    gps_out = sorted(gps_out, key=lambda x: x.date)
    return gps_out


def merge_gps_list(gps_from_log,gps_from_mer_environment) :
    '''

    Merge gps from log and gps from mermaid

    All are unique position but objects are not same

    '''
    gps_list = gps_from_log + gps_from_mer_environment
    gps_list = sorted(gps_list, key=lambda x: x.date)
    return gps_list


def write_gps(cycles, creation_datestr, processed_path, mfloat_path):
    '''

    Write complete (raw, full, all, nonunique) GPS data from .LOG and .MER.
    Differs from GeoCSV, which writes unique (merged .MER time and .LOG
    position) GPS fixes.

    '''

    gps_genexp = (gps for cycle in cycles for gps in cycle.gps_list)

    # Version and creation-date lines are the same for both csv and txt files
    version_line = "#automaid {} ({})\n".format(setup.get_version(), setup.get_url())
    created_line = "#created {}\n".format(creation_datestr)

    # Specify field headers of both csv and txt files
    header_line_txt = "           gps_time       gps_lat        gps_lon  gps_hdop  gps_vdop    gps_time-mer_time mer_clockfreq               source       raw_gps_lat        raw_gps_lon\n"
    header_line_csv = '#' + ','.join(header_line_txt.split()) + '\n'
    header_line_txt = '#' + header_line_txt # add pound sign after comma substitution

    # Specify generic format of both csv and txt files
    fmt = ['{:>19s}',
           '{:>10.6f}',
           '{:>11.6f}',
           '{:>6.3f}',
           '{:>6.3f}',
           '{:>17.6f}',
           '{:>10.0f}',
           '{:>17s}',
           '{:>14s}',
           '{:>15s}\n']

    # Add four spaces between each field for the txt file
    fmt_txt  = '    '.join(fmt)

    # Add comma between each field and remove field width (non-decimal) to format the csv
    fmt_csv  = ','.join(fmt)
    fmt_csv  = re.sub(':>\d*', ':', fmt_csv)

    # Specify file paths
    base_path = os.path.join(processed_path, mfloat_path)
    csv_file =  os.path.join(base_path, 'gps.csv')
    txt_file =  os.path.join(base_path, 'gps.txt')

    with open(csv_file, "w+") as f_csv, open(txt_file, "w+") as f_txt:
        # Write the version and header lines to both the csv and txt file
        f_csv.write(version_line)
        f_csv.write(created_line)
        f_csv.write(header_line_csv)

        f_txt.write(version_line)
        f_txt.write(created_line)
        f_txt.write(header_line_txt)

        for g in sorted(gps_genexp, key=lambda x: x.date):
            if g.hdop is None:
                g.hdop = float("NaN")
            if g.vdop is None:
                g.vdop = float("NaN")
            if g.clockdrift is None:
                g.clockdrift = float("NaN")
            if g.clockfreq is None:
                g.clockfreq = float("NaN")

            # Parse and format the raw strings.
            raw_lat = g.rawstr_dict['latitude']
            raw_lon = g.rawstr_dict['longitude']

            # Collect list of GPS data
            gps_data = [str(g.date)[:19] + 'Z',
                        g.latitude,
                        g.longitude,
                        g.hdop,
                        g.vdop,
                        g.clockdrift,
                        g.clockfreq,
                        g.source,
                        raw_lat,
                        raw_lon]

            # Write data to .csv and .txt formats
            f_csv.write(fmt_csv.format(*gps_data))
            f_txt.write(fmt_txt.format(*gps_data))


def write_gps_interpolation_txt(cycles, creation_datestr, processed_path, mfloat_path):
    '''Writes MERMAID GPS interpolation file, detailing GPS and interpolation parameters for the three
    main regimes of each dive: descent and drift in the surface layer, drift in the mixed layer, and
    ascent and drift in the surface layer.


    '''

    # NB, the comments here assume a normal dive where all GPS fixes are obtained and MERMAID dives
    # deeper than the mix layer depth; see especially dives.compute_station_locations and
    # gps.linear_interpolation to understand of edge-cases where perhaps some GPS fixes are missing
    # and/or MERMAID didn't dive into the mixed layer.  In all cases, GPS interpolation is still
    # broken into three regimes: descent drift, "deep" drift, and ascent drift.  Descent drift uses
    # the surface-drift velocity before the dive to interpolate forward in time for the location
    # where MERMAID dove into the mixed layer (left the surface layer); ascent drift uses the
    # surface-drift velocity after the dive to interpolate backward in time for the location where
    # MERMAID ascended into the surface layer (left the mixed layer); "deep" drift uses the velocity
    # of drift between those two points to estimate where MERMAID was when it recorded events while
    # drifting in the mixed layer.

    # "input" to gps.linear_interpolation are those GPS instances that we give the algorithm
    def parse_input_params(leg):
        input_params = [leg['input_drift_time']               if leg['input_drift_time'] else float("Nan"),
                        leg['input_drift_time'] / 60.0        if leg['input_drift_time'] else float("Nan"),
                        leg['input_drift_dist_m']             if leg['input_drift_dist_m'] else float("Nan"),
                        leg['input_drift_dist_m'] / 1000      if leg['input_drift_dist_m'] else float("Nan"),
                        leg['input_drift_vel_ms']             if leg['input_drift_vel_ms'] else float("Nan"),
                        leg['input_drift_vel_ms'] * 3.6       if leg['input_drift_vel_ms'] else float("Nan"), # km/hr
                        leg['input_drift_vel_ms'] * 3.6 * 24  if leg['input_drift_vel_ms'] else float("Nan")] # km/day

        input_params = map(abs, input_params)
        input_params = list(input_params)
        input_fmt_spec = '{:>6.0f}        {:>7.1f}        {:>6.0f}        {:>4.1f}        {:>5.2f}        {:>7.2f}        {:>7.2f}'

        return (input_params, input_fmt_spec)

    # "interp" from gps.linear_interpolation are those GPS instances the algorithm computes given
    # the input
    def parse_interp_params(leg):
        interp_params = [leg['interp_drift_time']              if leg['interp_drift_time'] else float("Nan"),
                        leg['interp_drift_time'] / 60.0        if leg['interp_drift_time'] else float("Nan"),
                        leg['interp_drift_dist_m']             if leg['interp_drift_dist_m'] else float("Nan"),
                        leg['interp_drift_dist_m'] / 1000      if leg['interp_drift_dist_m'] else float("Nan"),
                        leg['interp_drift_vel_ms']             if leg['interp_drift_vel_ms'] else float("Nan"),
                        leg['interp_drift_vel_ms'] * 3.6       if leg['interp_drift_vel_ms'] else float("Nan"), # km/hr
                        leg['interp_drift_vel_ms'] * 3.6 * 24  if leg['interp_drift_vel_ms'] else float("Nan")] # km/day

        interp_params = map(abs, interp_params)
        interp_params = list(interp_params)
        interp_fmt_spec = '{:>6.0f}        {:>7.1f}        {:>6.0f}        {:>4.1f}        {:>5.2f}        {:>7.2f}        {:>7.2f}'

        return (interp_params, interp_fmt_spec)

    # Generate (unique) list of cycle with events whose interpolated locations we are able to compute
    cycle_set = set(cycle for cycle in cycles for event in cycle.events if event.station_loc)

    # Print GPS interpolation information for every dive that includes an event all three dive regimes
    gps_interp_file = os.path.join(processed_path, mfloat_path, "gps_interpolation.txt")
    version_line = "#automaid {} ({})\n".format(setup.get_version(), setup.get_url())
    created_line = "#created {}\n".format(creation_datestr)

    with open(gps_interp_file, "w+") as f:
        f.write(version_line)
        f.write(created_line)

        for cycle in sorted(cycle_set, key=lambda x: x.start_date):

            # Compute the percentage of the total interpolate distance for the three regimes:
            # (1) surface-layer drift during the descent
            #
            # (2) mixed_layer drift
            #     .station.loc['interp_dist_m'] differs for each event  (drift to event in mixed layer)
            #     .station.loc['input_dist_m'] same for all events (total mixed-layer drift)
            #
            # (3) surface-layer drift during the ascent

            leg_descent = cycle.descent_last_loc_before_event
            leg_ascent = cycle.ascent_first_loc_after_event
            if leg_descent is None or leg_ascent is None:
                continue

            interp_dist_descent = leg_descent.interp_dict['interp_drift_dist_m']
            input_dist_mixed =  cycle.events[0].station_loc.interp_dict['input_drift_dist_m']
            interp_dist_ascent = leg_ascent.interp_dict['interp_drift_dist_m']

            if all([interp_dist_descent, input_dist_mixed, interp_dist_ascent]):
                bad_interp = False
                total_interp_dist = sum([interp_dist_descent, input_dist_mixed, interp_dist_ascent])
                interp_perc_descent = (interp_dist_descent / total_interp_dist) * 100
                input_perc_mixed = (input_dist_mixed / total_interp_dist) * 100
                interp_perc_ascent = (interp_dist_ascent / total_interp_dist) * 100

            else:
                bad_interp = True
                interp_perc_descent = float("nan")
                input_perc_mixed = float("nan")
                interp_perc_ascent = float("nan")

            # Write headers to each cycle block
            f.write("CYCLE ID: {:>4d} ".format(cycle.cycle_nb))
            f.write("DATES: {:>19s} --> {:19s}\n\n".format(str(cycle.start_date)[:19] + 'Z', str(cycle.end_date)[:19] + 'Z'))
            f.write("DRIFT_REGIME               TIME_S       TIME_MIN        DIST_M     DIST_KM      VEL_M/S      VEL_KM/HR     VEL_KM/DAY      DIST_%                                   SAC_MSEED_TRACE\n")

            # Parse the GPS ('input') components of surface drift before cycle: these are actual GPS points
            gps_surface_descent, gps_fmt_spec = parse_input_params(leg_descent.interp_dict)

            gps_fmt_spec = "gps_surface                " + gps_fmt_spec + "\n"
            f.write(gps_fmt_spec.format(*gps_surface_descent))

            # Parse the interpolated components of surface drift before cycle: between last GPS point
            # and crossing into mixed layer
            interp_surface_descent, interp_fmt_spec = parse_interp_params(leg_descent.interp_dict)
            interp_surface_descent.append(interp_perc_descent)

            interp_fmt_spec = "interp_surface             " + interp_fmt_spec + "        {:>4.1f}\n"
            f.write(interp_fmt_spec.format(*interp_surface_descent))

            # For every event recorded during the cycle: parse the interpolated components of the
            # mixed-layer drift from leaving the surface layer (passing into the "deep" or
            # mixed-layer drift regime) and recording an event
            for event in cycle.events:
                # if event.station_loc_is_preliminary:
                #     continue

                interp_drift_to_event_mixed_layer, interp_fmt_spec = parse_interp_params(event.station_loc.interp_dict)
                interp_drift_to_event_mixed_layer.append(event.processed_file_name)

                interp_fmt_spec = " interp_mixed(to_event)    " + interp_fmt_spec + "                    {:>42s}\n"
                f.write(interp_fmt_spec.format(*interp_drift_to_event_mixed_layer))

            # The total interpolated drift in the mixed layer -- that drift that occurs between the
            # last point of the ascent and the first point of the ascent -- is the same for every
            # event; just use the last event instance
            total_drift_mixed_layer, interp_fmt_spec = parse_input_params(cycle.events[0].station_loc.interp_dict)
            total_drift_mixed_layer.append(input_perc_mixed)

            interp_fmt_spec = "interp_mixed(total)        " + interp_fmt_spec + "        {:>4.1f}\n"
            f.write(interp_fmt_spec.format(*total_drift_mixed_layer))

            # Parse the interpolated components of surface drift after cycle: crossing out of mixed
            # layer and recording first GPS point
            interp_surface_ascent, interp_fmt_spec = parse_interp_params(leg_ascent.interp_dict)
            interp_surface_ascent.append(interp_perc_ascent)

            interp_fmt_spec = "interp_surface             " + interp_fmt_spec + "        {:>4.1f}\n"
            f.write(interp_fmt_spec.format(*interp_surface_ascent))

            # Parse the GPS ('input') components of surface drift after cycle: these are actual GPS points
            gps_surface_ascent, gps_fmt_spec = parse_input_params(leg_ascent.interp_dict)

            gps_fmt_spec = "gps_surface                " + gps_fmt_spec + "\n"
            f.write(gps_fmt_spec.format(*gps_surface_ascent))

            # If the interpolation failed, print some helpful statements at end of block
            if bad_interp:
                f.write('\n')
                if leg_descent.interp_dict['input_drift_dist_m'] is None:
                    f.write("*Interpolation issue before dive (surface-layer drift): {:s}\n" \
                            .format(leg_descent.interp_dict['description']))

                if cycle.events[0].station_loc.interp_dict['input_drift_dist_m'] is None:
                    f.write("*Interpolation issue during dive (mixed-layer drift): {:s}\n" \
                            .format(cycle.events[0].station_loc.interp_dict['description']))

                if leg_ascent.interp_dict['input_drift_dist_m'] is None:
                    f.write("*Interpolation issue after dive (surface-layer drift): {:s}\n" \
                            .format(leg_ascent.interp_dict['description']))

            f.write('\n__________END__________\n\n')
