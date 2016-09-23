#!/usr/bin/python
# dmioscope.py - IO frequency by disk location
#
# Copyright (C) 2016 Red Hat, Inc. All rights reserved.
#
# Bryn M. Reeves <bmr@redhat.com>
#
# This copyrighted material is made available to anyone wishing to use,
# modify, copy, or redistribute it subject to the terms and conditions
# of the GNU General Public License v.2.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#
# Usage:
#
#  dmioscope.py [-h] [-a] [-b NR_BINS] [-c] [-n] [-p] [-r ROWS] [-s]
#                      [-t THRESH] [-v] [-w WIDTH]
#                      [interval] [count] dev [dev ...]
#
#  IOScope arguments.
#
#  positional arguments:
#    interval              Duration of the report interval in seconds.
#    count                 Number of intervals to report.
#    dev                   Device(s) to monitor.
#
#  optional arguments:
#    -h, --help            show this help message and exit
#    -a, --adaptive        Adapt the number of bins to observed IO volume.
#    -b NR_BINS, --bins NR_BINS
#                          Divide devices into nr equally sized bins.
#    -c, --current         Show the current interval plot.
#    -n, --no-adaptive     Do not adapt the number of bins according to
#                          observed IO volume.
#    -p, --percent         Show distribution values as percentages.
#    -r ROWS, --rows ROWS  Specify the maxumum number of rows to use.
#    -s, --summary         Show the accumulated summary plot.
#    -t THRESH, --threshold THRESH
#                          Threshold at which to split a bin when using
#                          --adaptive
#    -v, --verbose         Enable verbose debugging messages.
#    -w WIDTH, --width WIDTH
#                          Specify the maximum terminal width to use.
#
from __future__ import print_function

from subprocess import Popen, PIPE, STDOUT
import argparse
import shlex
import sys
import os
import errno
import math
import time
import fcntl
import struct
import termios

_dm_report_fields = "region_id,read_count,write_count"
_dm_report_cmd = "dmstats report --noheadings -o"

READS_COUNT = 0
WRITES_COUNT = 1

_counters = [
    "READS_COUNT",
    "WRITES_COUNT"
]

CLEAR_SCREEN = "\033c"

_test = False
_verbose = False


class DmstatsException(Exception):
    """ Generic class representing errors communicating with dmstats.
    """
    pass


def log_info(str):
    print(str)


def log_verbose(str):
    if not _verbose:
        return
    print(str)


def log_error(str):
    print(str, file=sys.stderr)


def _terminal_size():
    """ Return the current terminal size as a ('w', 'h') tuple.
    """
    h, w, hp, wp = struct.unpack('HHHH',
                                 fcntl.ioctl(0, termios.TIOCGWINSZ,
                                             struct.pack('HHHH', 0, 0, 0, 0)))
    return w, h


def _device_sectors(dev_path):
    """ Return device size in 512b sectors.
    """
    req = 0x80081272  # BLKGETSIZE64, result is bytes as uint64_t.
    buf = ' ' * 8
    fmt = 'L'

    if os.path.sep in dev_path and not dev_path.startswith(os.path.sep):
        # vg/lv name
        dev_path = os.path.join("/dev", dev_path)
    else:
        # dm name
        os.path.join("/dev/mapper", dev_path)

    with open(dev_path) as dev:
        buf = fcntl.ioctl(dev.fileno(), req, buf)
    size_bytes = struct.unpack(fmt, buf)[0]
    return size_bytes >> 9

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
    if not maxval:
        maxval = nticks
    span = find_nice(maxval - minval, 0)
    d = find_nice(span / (nticks - 1), 1)
    graphmin = math.floor(minval / d) * d
    graphmax = math.ceil(maxval / d) * d
    # nfrac = max(-math.floor(math.log10(d)), 0)
    labels = []
    for x in frange(graphmin, graphmax + 0.5 * d, d):
        labels.append(x)
    return labels


def _sizeof_fmt(num, suffix='B'):
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


class Bin(object):
    count = 0
    start = 0
    width = 0

    def __init__(self, start, width, count=0):
        self.start = start
        self.width = width
        self.count = count

    def bound(self):
        """ Return the upper bound on this bin.
        """
        return self.start + self.width

    def split(self):
        """ Split the count of this bin uniformly in two.
        """
        self.count /= 2
        self.width /= 2
        return Bin(self.start + self.width, self.width, count=self.count)


_min_size_fraction = 40

# IOHistogram class constants
RENDER_COUNTS = 1
RENDER_TOTALS = 2


class IOHistogram(object):

    # IOHistogram instance variables
    counter = 0  # READS_COUNT
    nr_bins = 0
    dev_size = 0
    device = None
    bounds = None  # current bounds set
    regions = None  # current region_ids
    bins = None
    totals = None
    region_map = dict()

    def _init_bins(self):
        self.bins = []
        self.totals = []
        bounds = self.bounds

        start = 0
        nr_bins = 0
        for val in bounds:
            _bin = Bin(start, val - start)
            self.bins.append(_bin)
            _tot_bin = Bin(start, val - start)
            self.totals.append(_tot_bin)
            log_verbose("Initialised bin @ %d width %d"
                        % (start, (val - start)))
            start = val
            nr_bins += 1
        self.dev_size = start
        self.nr_bins = nr_bins
        return nr_bins

    def _make_bounds(self, nr_regions):
        """ Make a list of evenly spaced bounds values spanning an
            entire device.
        """
        if not self.dev_size:
            log_error("Found zero-sized device: %s" % device)
            return None

        step = self.dev_size / nr_regions

        bounds = [x for x in xrange(step, self.dev_size + 1, step)]
        log_verbose(bounds)
        return bounds

    def __init__(self, device, counter, bounds, adapt=True):
        """ Initialise an IOHistogram using the sector boundaries listed
            in bins. Boundary values are given as an upper bound on the
            current bin, with an implicit lower bound of 0 on the first.
            The final bin should span all remaining space on the device.

            All offsets are in 512b sectors.

            The new histogram is bound to the counter field specified by
            counter.

            For e.g.:

            ioh = IOHistogram([2048, 4096, 6144, 8192], READS_COUNT);

            Will create a histogram with four bins from 0..1M, 1M..2M,
            2M..3M, and 3M..4M, tracking the READS_COUNT counter.
        """
        self.device = device
        self.bounds = bounds
        log_verbose(self.bounds)

        self.adapt = adapt

        self.dev_size = _device_sectors(self.device)
        self.min_size = self.dev_size / _min_size_fraction

        nr_bins = self._init_bins()
        self.counter = counter
        log_info("Initialised %s histogram with %d bins, min_size=%d."
                 % (_counters[counter], nr_bins, self.min_size))

    def __init__(self, device, counter, initial_bins=1, adapt=True):
        """ Initialise an IOHistogram with inital_regions bins evenly
            spaced across the specified device.

            All offsets are in 512b sectors.

            The new histogram is bound to the counter field specified by
            counter.

            For e.g.:

            ioh = IOHistogram([2048, 4096, 6144, 8192], READS_COUNT);

            Will create a histogram with four bins from 0..1M, 1M..2M,
            2M..3M, and 3M..4M, tracking the READS_COUNT counter.
        """
        self.device = device

        self.dev_size = _device_sectors(self.device)
        self.min_size = self.dev_size / _min_size_fraction

        self.bounds = self._make_bounds(initial_bins)
        log_verbose(self.bounds)

        self.adapt = adapt

        nr_bins = self._init_bins()

        if self.adapt:
            hist_type = "adaptive"
        else:
            hist_type = "fixed"

        log_info("Initialised %s %s histogram with %d bins, min_size=%d."
                 % (_counters[counter], hist_type, nr_bins, self.min_size))

    def update(self, data, test=False):
        """ Populate or update the histogram using the string counter
            values in data. The current value of the histogram is set
            to the counter values in data and this is added to the
            historical totals for the histogram.
        """
        if test:
            return self._test_update()

        log_verbose("Updating %d histogram bins" % self.nr_bins)
        # data contains one row per histogram bin
        for line in data.splitlines():
            fields = line.split(":")
            counters = fields[1:]
            region = int(fields[0])
            log_verbose("%d: %s" % (region, counters[self.counter]))

            _bin = self.region_map[region]
            value = int(counters[self.counter])
            self.bins[_bin].count = value
            self.totals[_bin].count += value

    def _test_update(self):
        data = ""
        for i in xrange(self.nr_bins):
            data += "%d:%d\n" % (i, i)
        self.update(data)

    def min_width(self):
        return min([_bin.width for _bin in self.bins])

    def max_count(self):
        """ Return the maximum count value contained in any bin.
        """
        max_count = 0
        for _bin in self.bins:
            if _bin.count > max_count:
                max_count = _bin.count
        return max_count

    def max_freq(self, total=False):
        """ Return the maximum frequency value contained in any bin.
        """
        max_freq = 0.0
        min_width = self.min_width()

        if total:
            bins = self.totals
        else:
            bins = self.bins

        for _bin in bins:
            scale = _bin.width / min_width
            freq = float(_bin.count) / scale
            if freq > max_freq:
                max_freq = freq
        return max_freq

    def sum(self, total=False):
        """ Return the sum of all count values for all bins.
        """
        bin_sum = 0
        if total:
            bins = self.totals
        else:
            bins = self.bins

        for _bin in bins:
            bin_sum += _bin.count
        return bin_sum

    def io_distribution(self, percent, total=False):
        """ Return the proportion of disk reached by the specified
            percentage of IO requests.
        """
        size = 0
        count_sum = 0
        count = float(self.sum(total=total))
        thresh = (count * percent) / 100.0

        if total:
            bins = self.totals
        else:
            bins = self.bins

        sorted_bins = sorted(bins, key=lambda b: b.count, reverse=True)

        for _bin in sorted_bins:
            count_sum += _bin.count
            size += _bin.width
            if (count_sum > thresh):
                break

        return ((100.0 * size) / self.dev_size)

    def print_histogram(self, columns=80, render=RENDER_COUNTS):
        """ Print an ASCII representation of a histogram and its values,
            suitable for display on a terminal of at least 'colums'
            width.
        """

        if render & RENDER_TOTALS:
            max_freq = self.max_freq(total=True)
            hist_sum = self.sum(total=True)
        else:
            max_freq = self.max_freq(total=False)
            hist_sum = self.sum(total=False)

        x_label_width = 8
        vbar = "|"
        axis = "+"
        hbar = "-"

        row_width = self.min_width()

        prefix_width = x_label_width + len(vbar) + 2
        chars = columns - prefix_width

        labels = label_value_axis(0, max_freq, 5)

        axis_labels = " " * (prefix_width - len(str(labels[0])))
        label_width = int(round(float(chars) / (len(labels) - 1)))
        for label in labels:
            label_str = str(label)
            label_len = len(label_str)
            axis_labels += label_str
            axis_labels += (label_width - label_len) * " "

        counts_per_char = float(chars) / labels[-1]

        # regions is a list of either a 1-tuples, or a 2-tuples containing the
        # count and total bins for this region.
        if render == RENDER_COUNTS:
            regions = zip(self.bins)
        elif render == RENDER_TOTALS:
            regions = zip(self.totals)
        else:
            regions = zip(self.bins, self.totals)

        print(axis_labels)
        header = (1 + x_label_width) * " " + axis + chars * hbar
        print(header)

        for region in regions:
            _bin = region[0]
            if len(region) > 1:
                # stacked count and totals render
                _tot = region[1]
            else:
                # single count or totals render
                _tot = None

            scale = _bin.width / row_width
            label = _sizeof_fmt(_bin.start << 9)
            for step in xrange(scale):
                row = label + (1 + (x_label_width - len(label))) * " " + vbar
                row += int((_bin.count * counts_per_char) / scale) * "#"
                if _tot:
                    # plot accumulated total over count
                    count_diff = _tot.count - _bin.count
                    row += int(((count_diff) * counts_per_char) / scale) * "@"
                print(row)
                label = ""  # only label 1st row

        # FIXME: command line option
        # points = [50.0, 66.6, 75.0, 90.0, 95.0, 99.0]
        points = [90.0]
        for point in points:
            print("%.2f%% of IO reaches %.2f%% of disk."
                  % (point, self.io_distribution(point, total=False)))
        print("")

    def update_region_map(self):
        index = xrange(len(self.regions))
        self.region_map = dict(zip(self.regions, index))

    def update_bin_regions(self, merge=False):
        if not self.adapt:
            return

        log_verbose("Updating bin regions (nr_bins=%d)" % self.nr_bins)
        # zip bins, totals, and bounds into a tuple for splitting.
        inbins = zip(self.bins, self.totals, self.bounds, self.regions)
        outbins = []
        min_split = self.min_size * 2

        # copy from inbins to outbins, splitting and merging as needed.
        for inbin in inbins:
            (_bin, _tot, bound, region) = inbin
            if (_bin.count > _threshold and _bin.width >= min_split):
                log_verbose("  Splitting region @ %d" % _bin.start)
                newbin = _bin.split()
                newtot = _tot.split()

                newbound = newbin.bound()
                bound = _bin.bound()

                self._remove_bin_region(region)

                oldregion = self._create_bin_region(_bin.start, _bin.width)
                newregion = self._create_bin_region(newbin.start, newbin.width)

                outbins.append((_bin, _tot, bound, oldregion))
                outbins.append((newbin, newtot, newbound, newregion))

                self.nr_bins += 1
            else:
                # merge neighbouring bins?
                mergeable = (len(outbins) and
                             _bin.count <= _merge_threshold and
                             outbins[-1][0].count <= _merge_threshold)
                if (merge and mergeable):
                        lastbin = outbins.pop()
                        log_verbose("  Merging regions @ %d, %d" %
                                    (lastbin[0].start, _bin.start))

                        # current bin
                        lastbin[0].count += _bin.count
                        lastbin[0].width += _bin.width

                        # totals bin
                        lastbin[1].count += _tot.count
                        lastbin[1].width += _tot.width

                        # update bounds
                        newbound = lastbin[0].bound()

                        # update dmstats regions
                        # FIXME: this is simple but leads to many create /
                        # delete operations when large numbers of consecutive
                        # regions are merged.
                        self._remove_bin_region(region)
                        self._remove_bin_region(lastbin[3])

                        newregion = self._create_bin_region(lastbin[0].start,
                                                            lastbin[0].width)

                        outbins.append(
                            (lastbin[0], lastbin[1], newbound, newregion)
                        )

                        self.nr_bins -= 1

                else:
                    log_verbose("  Keeping region @ %d" % _bin.start)
                    outbins.append(inbin)

        (_bins, _tots, bounds, regions) = zip(*outbins)

        self.bins = list(_bins)
        self.totals = list(_tots)
        self.bounds = list(bounds)
        self.regions = list(regions)
        self.update_region_map()

        log_verbose("Updated bins: %d bins, %d totals, "
                    "%d bounds, %d regions." %
                    (len(self.bins), len(self.totals),
                     len(self.bounds), len(self.regions)))

        out = _get_cmd_output("dmstats list %s" % self.device)
        log_verbose(out)

    def _create_bin_region(self, start, length):
        out = _get_cmd_output("dmstats create --start %d --length %d %s"
                              % (start, length, self.device))
        if not out:
            return -1

        log_verbose("Created region_id %d @ %d length=%d)" %
                    (int(out.strip().split()[-1]), start, length))

        return int(out.strip().split()[-1])

    def create_bin_regions(self):
        start = 0
        regions = []
        for bound in self.bounds:
            regions.append(self._create_bin_region(start, bound - start))
            start = bound
        self.regions = regions
        self.update_region_map()

    def _remove_bin_region(self, region_id):
        dev_region = (region_id, self.device)
        out = _get_cmd_output("dmstats delete --regionid %d %s" % dev_region)
        log_verbose("Removed region_id %d from %s" % dev_region)

    def remove_bin_regions(self):
        for region in self.regions:
            self._remove_bin_region(region)

_log_commands = True


def _get_cmd_output(cmd, stderr=False):
    args = shlex.split(cmd)

    if _log_commands:
        log_verbose("_get_cmd_output('%s')" % cmd)

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

# threshold at which to split a bin in two.
_threshold = 5000
_merge_threshold = 0


def _parse_options(args):
    parser = argparse.ArgumentParser(description="IOScope arguments.")

    # positional parameters
    parser.add_argument("interval", metavar="interval", nargs="?", default=5,
                        help="Duration of the report interval in seconds.",
                        type=float)

    parser.add_argument("count", metavar="count", nargs="?", default=1,
                        help="Number of intervals to report.", type=int)

    parser.add_argument("devices", metavar="dev", nargs="+",
                        help="Device(s) to monitor.")

    # argument groups

    adaptp = parser.add_mutually_exclusive_group(required=False)

    # switches and options
    adaptp.add_argument("-a", "--adaptive", action="store_true",
                        dest="adaptive", default=True,
                        help="Adapt the number of bins to observed IO volume.")

    parser.add_argument("-b", "--bins", action="store", type=int,
                        dest="bins", metavar="NR_BINS", default=1,
                        help="Divide devices into nr equally sized bins.")

    parser.add_argument("--clear", action="store_true",
                        dest="clear", default=False, help="Clear the screen "
                        "before each update.")

    parser.add_argument("-c", "--current", action="store_true",
                        dest="current", default=True,
                        help="Show the current interval plot.")

    parser.add_argument("-m", "--merge", action="store_true",
                        dest="merge", default=False, help="Allow merging of "
                        "bins with low IO volume in adaptive mode.")

    parser.add_argument("-M", "--merge-threshold", action="store", type=int,
                        dest="merge_thresh", default=0, help="Threshold at "
                              "which to merge adjacent regions with low IO.",
                        metavar="MERGE_THRESH")

    adaptp.add_argument("-n", "--no-adaptive", action="store_false",
                        dest="adaptive", help="Do not adapt the number of bins"
                        " according to observed IO volume.")

    parser.add_argument("-p", "--percent", action="store_true",
                        dest="percent", default=False,
                        help="Show distribution values as percentages.")

    parser.add_argument("-r", "--rows", action="store", type=int,
                        dest="rows", default=None, help="Specify the maxumum "
                        "number of rows to use.")

    parser.add_argument("-s", "--summary", action="store_true",
                        dest="summary", default=False,
                        help="Show the accumulated summary plot.")

    parser.add_argument("-t", "--threshold", action="store", type=int,
                        dest="thresh", default=5000, help="Threshold at "
                        "which to split a bin when using --adaptive")

    parser.add_argument("-v", "--verbose", action="store_true",
                        dest="verbose", default=False, help="Enable verbose "
                        "debugging messages.")

    parser.add_argument("-w", "--width", action="store", type=int,
                        dest="width", default=None, help="Specify the maximum "
                        "terminal width to use.")

    return parser.parse_args()

_devices = []
_histograms = []


def _remove_all_regions():
    for ioh in _histograms:
        ioh.remove_bin_regions()


def main(argv):
    global _devices, _histograms, _merge_threshold, _threshold, _verbose

    args = _parse_options(argv)

    interval = args.interval
    count = args.count
    clear = args.clear

    _devices = args.devices
    _threshold = args.thresh
    _verbose = args.verbose

    adapt = args.adaptive
    merge = args.merge
    _merge_threshold = args.merge_thresh

    nr_bins = args.bins

    if not args.width or not args.rows:
        (w, h) = _terminal_size()
        width = w if not args.width else args.width
        rows = h if not args.rows else args.rows

    log_verbose("Starting %dx%d IOScope histogram" % (width, rows))
    log_verbose("adapt=%s, merge=%s, thresh=%d, merge_thresh=%d, "
                "nr_bins=%d" % (adapt, merge, _threshold,
                                _merge_threshold, nr_bins))

    width -= 8  # maximum length of final label

    if not count:
        count = -1

    for dev in _devices:
        out = _get_cmd_output("dmstats delete --allregions %s" % dev)
        ioh = IOHistogram(dev, READS_COUNT, initial_bins=nr_bins, adapt=adapt)
        ioh.create_bin_regions()
        _histograms.append(ioh)

    while count:
        for dev in _devices:
            # Sleep at the top of the interval to accumulate some initial
            # data to display.
            time.sleep(interval)
            if clear:
                print(CLEAR_SCREEN)
            cmdstr = _dm_report_cmd + _dm_report_fields + " %s"
            out = _get_cmd_output(cmdstr % dev)

            ioh.update(out, test=_test)

            if args.summary:
                print("%s: accumulated IO distribution" % dev)
                ioh.print_histogram(columns=width, render=RENDER_TOTALS)
            if args.current:
                print("%s: current IO distribution" % dev)
                ioh.print_histogram(columns=width)

            ioh.update_bin_regions(merge=merge)

            count -= 1

    _remove_all_regions()

if __name__ == '__main__':
    try:
        status = main(sys.argv[1:])
    except KeyboardInterrupt:
        status = 0
    except DmstatsException:
        status = 1
    finally:
        _remove_all_regions()

    sys.exit(status)

# vim: set et ts=4 sw=4 :
