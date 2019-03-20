#!/usr/bin/env python
"""
    tf_from_rec.py
    --------------

    Reads two or more tracks from a wavefile containing loops
    of broadband noise. 

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


def file_reader(filename, channels):
    file_base, file_ext = os.path.splitext(filename)
    data = []
    if file_ext.lower() == '.aup':
        import audacity
        aud = audacity.Aup(filename)
        for chan in channels:
            w = aud.get_channel_data(chan)
            data.append(w)
        sr = aud.rate
    else:
        sys.stderr.write("Format not recognized: %s" % file_ext)
        return
    return(sr, data)


def avgspec(w, nwind, idx_start=0, idx_end=None):
    nhop = nwind
    if idx_end is None:
        idx_end = len(w)

    specs = []
    for ii in range(0, idx_end-nwind, nhop):
        specs.append(np.fft.fft(w[ii:ii+nwind]))
    specs = np.array(specs)
    meanspec = np.mean(specs, axis=0)
    devspecs = specs-np.tile(meanspec, (specs.shape[0], 1))
    stdspec = np.sqrt(np.mean(np.abs(devspecs**2), axis=0))
    return meanspec, stdspec


def calc_specs(filedict,
               nwind=4096,
               nloops=None,
               flims=None):

    nl = nloops
    discard_loops_start = 2
    discard_loops_end = 2
    idx_start = nwind * discard_loops_start

    specs = []
    ferrs = []

    if flims is None:
        flims = [0, None]
    if flims[1] is None:
        flims[1] = np.inf

    for ii, (fn, chans) in enumerate(filedict.items()):
        sr, data = file_reader(fn, chans)
        fvec = np.arange(nwind)/nwind*sr
        idx = (fvec > flims[0]) & (fvec < flims[1])
        #print (fn,chans)
        for chno, w in zip(chans, data):
            print(fn, chno, len(w))
            if nloops is None:
                samples_from_end = nwind*discard_loops_end
                nl = (len(w)-idx_start-samples_from_end)//nwind
            idx_end = idx_start + nl*nwind
            av_spec, std_spec = avgspec(w, nwind=nwind,
                                        idx_start=idx_start,
                                        idx_end=idx_end)
            specs.append(av_spec[idx])
            ferrs.append(std_spec[idx])

    return fvec[idx], specs, ferrs


def transfer_function(filedict,
                      nwind=4096,
                      nloops=None,
                      flims=None,
                      plot_error=False):

    fvec, specs, ferrs = calc_specs(filedict, nwind=nwind,
                                    nloops=nloops, flims=flims)
    rspec = specs[0]
    rerr = ferrs[0]

    fig, ax = pl.subplots(2, sharex=True)
    idx = (fvec > flims[0]) & (fvec < flims[1])

    for spec, espec in zip(specs[1:], ferrs[1:]):
        tf = spec/rspec
        lns = ax[0].semilogy(fvec[idx], np.abs(tf[idx]))
        if plot_error:
            color = lns[0].get_color()
            amin = (np.abs(spec)-np.abs(espec))/(np.abs(rspec)+np.abs(rerr))
            amax = (np.abs(spec)+np.abs(espec))/(np.abs(rspec)-np.abs(rerr))
            ax[0].fill_between(fvec[idx], amin, amax, alpha=0.3, color=color)
        ax[1].plot(fvec[idx], np.angle(tf[idx]))
    ax[0].set_ylabel("$|TF|$")
    ax[1].set_ylabel("$\\angle(TF)$")
    ax[1].set_xlabel("Frequency (Hz)")
    pl.show()


def loss_factor(fvec, spec1, spec2, dist=0.1, c=345.0):
    k = 2*np.pi * fvec / c
    kL = k*dist
    sinkl = np.sin(kL)
    ejkl = np.exp(1j*kL)
    #p,u=io.get_pressure_flow()
    #pout = p+u*io.parameters.z0
    #pin = p-u*io.parameters.z0
    pout = (spec2 - spec1*np.conjugate(ejkl))/sinkl/1j
    pin = (spec2 - spec1*ejkl)/(-sinkl)/1j

    return 1-np.abs(pin/pout)**2


def attenuation(filedict,
                nwind=4096,
                nloops=None,
                flims=None,
                mic_dist=0.1,
                speed_of_sound=345.0):

    fvec, specs, ferrs = calc_specs(filedict, nwind=nwind,
                                    nloops=nloops, flims=flims)
    fig, ax = pl.subplots(1, sharex=True)
    for spec2, spec1 in zip(specs[0::2], specs[1::2]):
        lf = loss_factor(fvec, spec1, spec2, dist=mic_dist, c=speed_of_sound)
        ax.plot(fvec, np.abs(lf))
    pl.ylabel("loss factor")
    pl.xlabel("Frequency (Hz)")
    fig.canvas.toolbar.push_current()  # save the 'un zoomed' view to stack
    ax.set_ylim([0, 1])
    fig.canvas.toolbar.push_current()  # save 'zoomed' view to stack

    pl.show()


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

    if args.attenuation is None:
        transfer_function(args.filedict,
                          nwind=args.loop_length,
                          nloops=args.num_loops,
                          flims=[float(xx) for xx in args.freq_lims],
                          plot_error=args.err)
    else:
        att_args = args.attenuation
        try:
            mic_dist = float(att_args[0])
        except IndexError:
            mic_dist = 0.1
            sys.stderr.write("Mic distance set to %f\n" % mic_dist)

        try:
            speed_of_sound = float(att_args[1])
        except IndexError:
            speed_of_sound = 345.0
            sys.stderr.write("Speed of sound set to %f\n" % speed_of_sound)
        attenuation(args.filedict,
                    nwind=args.loop_length,
                    nloops=args.num_loops,
                    flims=[float(xx) for xx in args.freq_lims],
                    mic_dist=mic_dist,
                    speed_of_sound=speed_of_sound)


