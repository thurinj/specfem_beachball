#!/usr/bin/env python

import pygmt
import obspy
import numpy as np
from obspy.taup import TauPyModel
import os
from pygmt.helpers import build_arg_string, use_alias
import argparse

parser = argparse.ArgumentParser(description='Get the path for CMTSOLUTION and STATIONS files.')
parser.add_argument('CMTSOLUTION', type=str, help='path for CMTSOLUTION')
parser.add_argument('STATIONS', type=str, help='path for STATIONS file')
args = parser.parse_args()

event_file = args.CMTSOLUTION
stations_file = args.STATIONS


class Station():
    # initialize the class
    def __init__(self, id, network, latitude, longitude):
        self.id = id
        self.network = network
        self.latitude = latitude
        self.longitude = longitude

    # define a string representation of the object
    def __str__(self):
        return "{} - {} - {} - {}".format(self.id, self.network, self.latitude, self.longitude)


class CMTSolution():
    """
    Class to hold the information from a CMTSOLUTION file.
    """

    def __init__(self, cmt_file):
        """
        Initialize event from a given CMTSOLUTION file.
        """
        with open(cmt_file, 'r') as f:
            lines = f.readlines()
            self.event_id = lines[0].split()[-1]
            self.time_shift = float(lines[2].split()[-1])
            self.half_duration = float(lines[3].split()[-1])
            self.latitude = float(lines[4].split()[-1])
            self.longitude = float(lines[5].split()[-1])
            self.depth = float(lines[6].split()[-1])
            self.m_rr = float(lines[7].split()[1])
            self.m_tt = float(lines[8].split()[1])
            self.m_pp = float(lines[9].split()[1])
            self.m_rt = float(lines[10].split()[1])
            self.m_rp = float(lines[11].split()[1])
            self.m_tp = float(lines[12].split()[1])
            # Read in the full moment tensor
            moment_tensor = np.zeros((3, 3))
            moment_tensor_vector = np.zeros(6)
            moment_tensor[0][0] = self.m_rr
            moment_tensor[0][1] = self.m_rp
            moment_tensor[1][0] = self.m_rp
            moment_tensor[1][1] = self.m_tt
            moment_tensor[0][2] = self.m_rt
            moment_tensor[2][0] = self.m_rt
            moment_tensor[1][2] = self.m_tp
            moment_tensor[2][1] = self.m_tp
            moment_tensor[2][2] = self.m_pp
            self.tensor = moment_tensor
            moment_tensor_vector[0] = self.m_rr
            moment_tensor_vector[1] = self.m_tt
            moment_tensor_vector[2] = self.m_pp
            moment_tensor_vector[3] = self.m_rt
            moment_tensor_vector[4] = self.m_rp
            moment_tensor_vector[5] = self.m_tp
            self.tensor_vector = moment_tensor_vector


def stations_list(filename):
    stations = []
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            stations.append(line.split())
    stations_coordinates = []
    for station in stations:
        stations_coordinates.append(Station(station[0], station[1], station[2], station[3]))
    return stations_coordinates


def beachball_pygmt(filename, stations, mt, plot_all=False):
    mt = event.tensor_vector
    focal_mechanism = np.append(np.append([0, 0, 10], mt), [25, 0, 0])

    # Initialize the pygmt plot with the beachball plot
    fig = pygmt.Figure()
    fig.meca(region=[-1.5, 1.5, -1.5, 1.5], projection='m0/0/5c', scale='3c',
             convention="mt", G='grey50', spec=focal_mechanism, N=False)
    # Launch the GMT polar function from pygmt wrapper.
    _pygmt_polar(stations, symbol='t0.40c', comp_fill='black', mt_outline=None)

    # fig.show(dpi=300, method="external")
    fig.savefig(filename, show=True)

# Define aliases for the pygmt function. Please refer to GMT 6.2.0 `polar` function documentation for a complete overview of all the available options and option details.
@use_alias(
    D='offset',
    J='projection',
    M='size',
    S='symbol',
    E='ext_fill',
    G='comp_fill',
    F='background',
    Qe='ext_outline',
    Qg='comp_outline',
    Qf='mt_outline',
    T='station_labels'
)
def _pygmt_polar(trace_list, **kwargs):
    """ Wrapper around GMT polar function. ]
    Color arguments must be in {red, green, blue, black, white} and the symbol in {a,c,d,h,i,p,s,t,x} (see GMT polar function for reference).
    """

    # Define some default values to format the plot.
    defaultKwargs = {
        'D': '0/0',
        'J': 'm0/0/5c',
        'M': '12.92c',
        'T': '+f0.18c'
    }
    kwargs = {**defaultKwargs, **kwargs}

    colorcodes = {
        "red": "255/0/0",
        "green": "0/255/0",
        "blue": "0/0/255",
        "white": "255, 255, 255",
        "black": "0/0/0"
    }
    for key in kwargs:
        try:
            kwargs[key] = colorcodes[kwargs[key]]
        except:
            pass

    tmp_filename = 'polar_temp.txt'
    with open(tmp_filename, 'w') as f:
        for sta in trace_list:
            pol = '+'
            f.write('{} {} {} {}'.format(sta.network+'.'+sta.id, sta.azimuth, sta.takeoff_angle, pol))
            f.write('\n')

    arg_string = " ".join([tmp_filename, build_arg_string(kwargs)])
    with pygmt.clib.Session() as session:
        session.call_module('polar', arg_string)

    os.remove(tmp_filename)

################################################################################
#                           Code starts here                                   #
################################################################################


station_coords = stations_list(stations_file)
event = CMTSolution(event_file)

taupmodel = TauPyModel(model="ak135")
for sta in station_coords:
    sta_lat = float(sta.latitude)
    sta_lon = float(sta.longitude)
    dist_in_m, az, baz = obspy.geodetics.base.gps2dist_azimuth(
        event.latitude, event.longitude, sta_lat, sta_lon)
    sta.back_azimuth = baz
    sta.azimuth = az
    sta.dist = dist_in_m * 1e-3
    sta.dist_in_degree = obspy.geodetics.base.kilometers2degrees(sta.dist)
    arrivals = taupmodel.get_travel_times(
        source_depth_in_km=event.depth, distance_in_degree=sta.dist_in_degree, phase_list=["P", "p"])
    first_arrival = arrivals[0]
    sta.takeoff_angle = first_arrival.takeoff_angle
    sta.incidence_angle = first_arrival.incident_angle
    sta.arrival_time = first_arrival.time

beachball_pygmt('big_beachball_test.pdf', station_coords, event)
