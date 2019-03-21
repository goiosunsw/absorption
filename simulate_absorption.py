#!/usr/bin/env python
"""
    simulate_absorption.py
    ----------------------

    Generates simulated absorption measurement signals
    and plots absorption calculated using tf_from_rec.py

    Plots a transfer function between the two tracks or 
    the attenuation coefficient based on those two tracks.

    Author: Andre Almeida <a.almeida@unsw.edu.au> 
    Creation date: 2019/03/18
    Version: 0.0.1
"""

import os
import sys
import argparse

import numpy as np
import matplotlib.pyplot as pl

from tf_from_rec import loss_factor

def generate_specs(mic_dist=0.1, sr=44100,
                   sample_dist=0.3, ls_dist=0.3,
                   reflection_coef = 1.0,
                   speed_of_sound=345.0,
                   loop_length=4096, 
                   fmin=100, fmax=4000):

    f_orig = np.linspace(0, sr, loop_length)
    f = f_orig[(f_orig>fmin) & (f_orig<fmax)]
    k = 2*np.pi*f/speed_of_sound

    p_out_ls = np.ones(len(f),dtype='complex')
    total_dist = mic_dist+ls_dist+sample_dist
    p_out_smpl = p_out_ls * np.exp(1j*k*total_dist)
    p_in_smpl = p_out_smpl * reflection_coef
    
    p_out_1 = p_out_ls * np.exp(1j*k*(ls_dist+mic_dist))
    p_out_2 = p_out_ls * np.exp(1j*k*(mic_dist))
    p_in_1 = p_in_smpl * np.exp(1j*k*sample_dist)
    p_in_2 = p_in_1 * np.exp(1j*k*mic_dist)

    p_1 = p_in_1 + p_out_1
    p_2 = p_in_2 + p_out_2

    return f, p_1, p_2


def parse_arguments():
    ap1 = argparse.ArgumentParser()
    ap1.add_argument("-l", "--loop-length",  nargs='?', default=4096, type=int,
                     help="number of samsples per loop")
    ap1.add_argument("-n", "--num-loops",  nargs='?', type=int,
                     help="number of loops to use in recording")
    ap1.add_argument("-f", "--freq-lims",  nargs=2, default=[0, None],
                     help="minimum and maximum frequency of excitation")
    ap1.add_argument("-e", "--err",  action='store_true', default=False,
                     help="plot spectral error")
    ap1.add_argument("-a", "--attenuation", default=None, nargs='*',
                     help="""calculate attenuation instead of transfer function. Arguments
                     are microphone distance and speed of sound""")
    #ap2 = argparse.ArgumentParser()
    ap1.add_argument("inputs",  nargs='*',
                     help="sequence of file names and channels")
    ns, rest = ap1.parse_known_args()

    # ns = ap2.parse_args(rest, ns)

    filedict = {}
    f = ''
    chans = []
    for arg in ns.inputs:
        try:
            ch = int(arg)
            chans.append(ch)
        except ValueError:
            try:
                filedict[f].extend(chans)
            except KeyError:
                if len(f) > 0:
                    filedict[f] = chans
            f = arg
            chans = []

    try:
        filedict[f].extend(chans)
    except KeyError:
        if len(f) > 0:
            filedict[f] = chans

    filedict[f] = chans

    print(filedict)
    ns.filedict = filedict
    return ns


if __name__ == '__main__':
    args = parse_arguments()

    f, p1, p2 = generate_specs()

    l = loss_factor(f,p1,p2)

    pl.figure()
    pl.plot(f,l)
    pl.show()
