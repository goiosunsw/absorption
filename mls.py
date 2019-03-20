#!/usr/bin/env python
"""Module providing the functionality to generate Maximal Length Sequences. 
Based on wiki's description and polynomial representation.
http://en.wikipedia.org/wiki/Maximum_length_sequence

Main code from
https://github.com/mubeta06/python/blob/master/signal_processing/sp/mls.py

"""

import numpy


bittaps = {2: [1],
           3: [2],
           4: [3],
           5: [2],
           6: [5],
           7: [6],
           8: [4, 5, 6],
           9: [5],
           10: [7],
           11: [9]}


class Error(Exception):

    """This is a user-defined exception for errors raised by this module."""

    pass


def lfsr(taps, buf):
    """Function implements a linear feedback shift register
    taps:   List of Polynomial exponents for non-zero terms other than 1 and n
    buf:    List of buffer initialisation values as 1's and 0's or booleans
    """
    nbits = len(buf)
    sr = numpy.array(buf, dtype='bool')
    out = []
    for i in range((2**nbits)-1):
        feedback = sr[0]
        out.append(feedback)
        for t in taps:
            feedback ^= sr[t]
        sr = numpy.roll(sr, -1)
        sr[-1] = feedback
    return out


def mls(n, seed=None):
    """Generate a Maximal Length Sequence 2^n - 1 bits long
    """
    if n in bittaps:
        taps = bittaps[n]
    else:
        raise Error('taps for %s bits unknown' % str(n))
    if seed is None:
        seed = list(numpy.ones(n))
    elif len(seed) != n:
        raise Error('length of seed must equal n')
    return numpy.where(lfsr(taps, seed), 1, -1)


def main(nbits):
    """Main Program"""
    import pylab
    nbits = int(nbits)
    m = mls(nbits)
    fm = pylab.fft(m)
    fig, ax = pylab.subplots(2, sharex=True)
    ax[0].plot(20 * pylab.log10(numpy.abs(fm)))
    ax[1].plot(numpy.angle(fm))


if __name__ == '__main__':
    import sys
    main(sys.argv[1])
