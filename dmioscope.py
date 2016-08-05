#!/usr/bin/python
# dmshist.py - IO frequency by disk location
#
# Usage: dmshist.py <dm_device>
#
from __future__ import print_function

from subprocess import Popen, PIPE, STDOUT
import shlex
import sys
import os
import errno
import math
import time

_dm_report_fields = "read_count,write_count"
_dm_report_cmd = "dmstats report --noheadings -o"

READ_COUNT = 0
WRITE_COUNT = 1

_counters = [
    "READ_COUNT",
    "WRITE_COUNT"
]

_test=False
_verbose=False

def log_info(str):
    print(str)

def log_verbose(str):
    if not _verbose:
        return
    print(str)

def log_error(str):
    print(str, file=sys.stderr)


# Heckbert's Axis Labelling Algorithm ("nice numbers"), from:
# "Graphics Gems", Paul S. Heckbert, 1988.

def frange(start, end, step):
    """ Yield a range of floats, beginning at 'start' and
        incrementing by 'step' until 'end'.
    """
    while start <= end:
        yield start
        start += step

def find_nice(x, round_val):
    """ Find a "nice" number that is approximately equal to 'x'.
        Round the number if 'round_val'=1, or take ceiling otherwise.
    """
    expv = math.floor(math.log10(x))
    f = x / math.pow(10, expv)
    if (round_val):
        if (f < 1.5):
            nf = 1
        elif (f < 3):
            nf = 2
        elif (f < 7):
            nf = 5
        else:
            nf = 10
    else:
        if (f <= 1):
            nf = 1
        elif (f <= 2):
            nf = 2
        elif (f <= 5):
            nf = 5
        else:
            nf = 10
    return nf * math.pow(10, expv)

def label_value_axis(minval, maxval, nticks):
    span = find_nice(maxval - minval, 0)
    log_verbose("find_nice(%d, 1) = %d" % (maxval - minval, span))
    d = find_nice(span / (nticks - 1), 1)
    log_verbose("find_nice(%d, 1) = %d" % (span / (nticks - 1), d))
    graphmin = math.floor(minval / d) * d
    graphmax = math.ceil(maxval / d) * d
    #nfrac = max(-math.floor(math.log10(d)), 0)
    labels = []
    for x in frange(graphmin, graphmax + 0.5 * d, d):
        labels.append(x)
    return labels

def _sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

class Bin(object):
    count = 0
    start = 0
    width = 0

    def __init__(self, start, width):
        self.start = start
        self.width = width

    def split(self):
        """ Split the count of this bin uniformly into two sub-bins """
        return self.count / 2


class IOHistogram(object):
    counter = 0 # READ_COUNT
    nr_bins = 0
    dev_size = 0
    bins = None

    def __init__(self, bounds, counter):
        """ Initialise an IOHistogram using the sector boundaries listed
            in bins. Boundary values are given as an upper bound on the
            current bin, with an implicit lower bound of 0 on the first.
            The final bin should span all remaining space on the device.

            All offsets are in 512b sectors.

            The new histogram is bound to the counter field specified by
            counter.

            For e.g.:

            ioh = IOHistogram([2048, 4096, 6144, 8192], READ_COUNT);

            Will create a histogram with four bins from 0..1M, 1M..2M,
            2M..3M, and 3M..4M, tracking the READ_COUNT counter.
        """
        start = 0
        nr_bins = 0
        self.bins = []
        for val in bounds:
            _bin = Bin(start, val - start)
            self.bins.append(_bin)
            log_info("Initialised bin @ %d width %d"
                     % (start, (val - start)))
            start = val
            nr_bins += 1
        self.dev_size = start
        self.nr_bins = nr_bins
        log_info("Initialised %s histogram with %d bins."
                 % (_counters[counter], nr_bins))

    def _parse_data(self, data, update):
        """ Populate or update the histogram using the string counter
            values in data. If update is non-zero counter values will
            not be reset and new values are summed to yield a
            cumulative distribution.
        """
        _bin = 0
        # data contains one row per histogram bin
        for line in data.splitlines():
            counters = line.split(":")

            if _bin > self.nr_bins:
                log_error("Unexpected row in histogram data")
                return False

            if not update:
                self.bins[_bin].count = 0

            self.bins[_bin].count += int(counters[self.counter])
            _bin += 1

    def populate(self, data):
        """ Populate the histogram using the string counter values in
            data. Any existing counter values are discarded.
        """
        return self._parse_data(data, 0)

    def update(self, data):
        """ Update the histogram using the string counter values in
            data. Counter values are not reset and new values are
            summed to yield a cumulative distribution.
        """
        return self._parse_data(data, 1)

    def max_count(self):
        """ Return the maximum count value contained in any bin.
        """
        max_count = 0
        for _bin in self.bins:
            if _bin.count > max_count:
                max_count = _bin.count
        return max_count

    def sum(self):
        """ Return the sum of all count values for all bins.
        """
        bin_sum = 0
        for _bin in self.bins:
            bin_sum += _bin.count
        return bin_sum

    def io_distribution(self, percent):
        """ Return the proportion of disk reached by the specified
            percentage of IO requests.
        """
        size = 0
        total = 0
        count = float(self.sum())
        thresh = (count * percent) / 100.0
        sorted_bins = sorted(self.bins, key=lambda b: b.count, reverse=True)

        for _bin in sorted_bins:
            total += _bin.count
            size += _bin.width;
            if (total > thresh):
                break

        return ((100.0 * size) / self.dev_size)

    def print_histogram(self, columns=80):
        """ Print an ASCII representation of a histogram and its values,
            suitable for display on a terminal of at least 'colums'
            width.
        """
        max_count = self.max_count()
        hist_sum = self.sum()
        x_label_width = 8
        vbar = "|"
        axis = "+"
        hbar = "-"

        row_width = min([_bin.width for _bin in self.bins])

        if not hist_sum:
            return None

        prefix_width = x_label_width + len(vbar) + 2
        chars = columns - prefix_width

        axis_labels = " " * (prefix_width - 1)
        labels = label_value_axis(0, max_count, 5)
        label_width = chars / (len(labels) - 1)
        for label in labels:
            label_str = str(label)
            label_len = len(label_str)
            axis_labels += label_str
            axis_labels += (label_width - label_len) * " "

        counts_per_char = int(math.ceil(labels[-1] / float(chars)))

        print(axis_labels)            
        header = (1 + x_label_width) * " " + axis + chars * hbar
        print(header)
        for _bin in self.bins:
            scale = _bin.width / row_width
            label = _sizeof_fmt(_bin.start << 9)
            for step in xrange(scale):
                row = label + (1 + (x_label_width - len(label))) * " " + vbar
                row += (_bin.count / (scale * counts_per_char)) * "#"
                print(row)
                label = "" # only label 1st row

        print("95%% of IO reaches %f%% of disk." % self.io_distribution(95.0))
        print("90%% of IO reaches %f%% of disk." % self.io_distribution(90.0))
        print("75%% of IO reaches %f%% of disk." % self.io_distribution(75.0))
        print("50%% of IO reaches %f%% of disk." % self.io_distribution(50.0))
        print("")

def _get_cmd_output(cmd,stderr=False):
    args = shlex.split(cmd)
    print(cmd)
    try:
        p = Popen(args, shell=False, stdout=PIPE,
                  stderr=STDOUT if stderr else PIPE,
                  bufsize=-1, close_fds=True)
        stdout, stderr = p.communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return None
        else:
            raise e

    if p.returncode != 0:
        log_verbose(stdout)
        return None

    return stdout


def main(args):
    _dev = args[1]
    interval = 3

    out = _get_cmd_output("dmstats delete --allregions %s" % _dev)

    # Uniform 4GiB bins through a 32GiB device.
    bounds = [
        16777216, 25165824, 33554432,
        41943040, 50331648, 58720256, 67108864
    ]

    ioh = IOHistogram(bounds, READ_COUNT);
    ioh_cum = IOHistogram(bounds, READ_COUNT);

    start = 0
    for bound in bounds:
        out = _get_cmd_output("dmstats create --start %d --length %d %s"
                              % (start, bound - start, _dev))
        start = bound

    while True:
        time.sleep(interval)

        out = _get_cmd_output(_dm_report_cmd + _dm_report_fields)
        print(out)
        ioh.populate(out)
        ioh_cum.update(out)
        ioh.print_histogram(columns=80)
        ioh_cum.print_histogram(columns=80)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        log_error("%s: missing file argument." % sys.argv[0])
        sys.exit(1)
    main(sys.argv)

# vim: set et ts=4 sw=4 :
