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
from builtins import range
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


def log_info(msg):
    print(msg)


def log_verbose(msg):
    if not _verbose:
        return
    print(msg)


def log_warn(msg):
    print("WARNING: " + msg, file=sys.stderr)


def log_error(msg):
    print(msg, file=sys.stderr)


def _terminal_size():
    """ Return the current terminal size as a ('w', 'h') tuple.
    """
    h, w, hp, wp = struct.unpack('HHHH',
                                 fcntl.ioctl(0, termios.TIOCGWINSZ,
                                             struct.pack('HHHH', 0, 0, 0, 0)))
    return w, h


def _get_cmd_output(cmd):
    """ Call `cmd` via `Popen` and return the status and combined `stdout`
        and `stderr` as a 2-tuple, e.g.:

        (0, "vg00/lvol0: Created new region with 1 area(s) as region ID 5\n")

    """
    args = shlex.split(cmd)

    if _log_commands:
        log_verbose("_get_cmd_output('%s')" % cmd)

    try:
        p = Popen(args, shell=False, stdout=PIPE,
                  stderr=STDOUT, bufsize=-1, close_fds=True)
        # stderr will always be None
        (stdout, stderr) = p.communicate()
    except OSError as e:
        log_error("Could not call '%s': %s" % (args[0], e))
        raise DmstatsException

    if _log_commands or p.returncode != 0:
        log_verbose(stdout.strip())

    return (p.returncode, stdout)


class DmstatsException(Exception):
    """ Generic class representing errors communicating with dmstats.
    """
    pass

# device-mapper executables
DMSETUP = "dmsetup"
DMSTATS = "dmstats"

# dmstats command verbs
DMS_CREATE = "create"
DMS_DELETE = "delete"
DMS_REPORT = "report"
DMS_LIST = "list"

# dmstats command verb argument templates
DMS_CREATE_ARGS = "--start %d --length %d %s"
DMS_DELETE_ARGS = "--regionid %d %s"
DMS_REPORT_ARGS = "--noheadings -o %s %s"
DMS_LIST_ARGS = ""


class DmStats(object):
    """ Class to encapsulate interaction with device-mapper statistics.
    """

    # bound device
    device = None
    uuid = None
    major = -1
    minor = -1

    # command state
    verb = None
    args = None
    status = 0

    def __init__(self, device, uuid=None, major=None, minor=None):
        """ Initialise a `DmStats` object and bind it to the specified
            DM device name, UUID, or major/minor pair. An error is logged
            and a DmstatsException raised if the device does not exist.
        """
        info_cmdstr = DMSETUP + " info -c --noheadings -o name "
        if device:
            result = _get_cmd_output(info_cmdstr + device)
        if uuid:
            self.uuid = uuid
            device = "with UUID " + uuid
            result = _get_cmd_output(info_cmdstr + "-u " + uuid)
        if major and minor:
            self.major = major
            self.minor = minor
            device = "with device no. %d:%d" % (major, minor)
            maj_min_opt = "-j %d -m %d" % (major, minor)
            result = _get_cmd_output(info_cmdstr + maj_min_opt)
        if result[0]:
            log_error("Error getting device information.")
            log_error("device %s does not exist" % device)
            raise DmstatsException

        # canonical dm name
        self.device = result[1].strip()

    def _issue_command(self):
        """ Issue the command configured in `verb` and `args`, and return
            the combined output as a string. The exit status is stored in
            `status`. Called by command methods.
        """
        cmdstr = "%s %s %s" % (DMSTATS, self.verb, self.args)
        result = _get_cmd_output(cmdstr)
        self.status = result[0]
        return result[1]

    def _status(self, message):
        """ Check the status of the last command and raise a DmstatsException
            if it is non-zero. The `message` argument specified a message to
            be logged in the event of error.
        """
        if self.status:
            log_error(message)
            raise DmstatsException

    def create(self, start, length):
        """ Call `dmstats` to create a new region on the bound device,
            beginning at `start`, and `length` sectors in size.

            Returns the region_id of the newly created region on success,
            or raises DmstatsException on error.
        """
        self.verb = DMS_CREATE
        self.args = DMS_CREATE_ARGS % (start, length, self.device)

        out = self._issue_command()

        self._status("Could not create region on device %s" % self.device)

        region_id = int(out.strip().split()[-1])

        log_verbose("Created region_id %d @ %d length=%d)" %
                    (region_id, start, length))

        return region_id

    def delete(self, region_id):
        """ Call `dmstats` to delete the specified `region_id` from the
            bound device.
        """
        dev_region = (region_id, self.device)
        self.verb = DMS_DELETE
        self.args = DMS_DELETE_ARGS % dev_region

        self._issue_command()

        self._status("Could not delete region_id %d from %s" % dev_region)

        log_verbose("Removed region_id %d from %s" % dev_region)

    def list(self):
        """ Call `dmstats` to list the regions present on the bound device.
        """
        self.verb = DMS_LIST
        self.args = DMS_LIST_ARGS

        out = self._issue_command()

        self._status("Could not retrieve region list for device %s." %
                     self.device)

        return out

    def report(self, fields=_dm_report_fields):
        """ Call `dmstats` to obtain current counter values for the bound
            device.
        """
        self.verb = DMS_REPORT
        self.args = DMS_REPORT_ARGS % (fields, self.device)

        out = self._issue_command()

        self._status("Could not retrieve counter data for regions "
                     "on device %s." % self.device)

        return out

_dmstats_counters = {
    "READS": "read_count",
    "READS_MERGED": "reads_merged_count",
    "READ_SECTORS": "read_sector_count",
    "READ_TIME": "read_time",
    "READ_TICKS": "read_ticks",
    "WRITES": "write_count",
    "WRITES_MERGES": "writes_merged",
    "WRITE_SECTORS": "write_sector_count",
    "WRITE_TIME": "write_time",
    "WRITE_TICKS": "write_ticks",
    "IN_PROGRESS": "in_progress_count",
    "IO_TICKS": "io_ticks",
    "QUEUE_TICKS": "queue_ticks"
}


class DmStatsCounters(object):
    """ Device-mapper statistics counter data.
    """
    counters = None

    def __init__(self, counters):
        """ Initialise a `DmStats` object with the specified list of fields.
        """
        self.counters = set()
        self._parse_counters(counters)

    def __repr__(self):
        """ Return a string represenatation of this `DmStatsCounters` object.
            This representation is identical to that parsed by the class
            constructor and can be used for display, or to initialise a new
            DmStatsCounters object.
        """
        if not self.counters:
            return ""
        return ",".join(list(self.counters))

    def _parse_counters(self, counters):
        for counter in counters.split(","):
            if counter not in _dmstats_counters:
                log_error("Unknown dmstats counter: %s" % counter)
                raise DmstatsException
            if counter in self.counters:
                log_warn("Counter %s specified twice." % counter)
            self.counters.add(counter)

    def fields(self):
        """ Return a comma-separated string list of field names, in the
            format expected by `dmstats report -o`.
        """
        # region_id is always the 0th field
        fields = ["region_id"] + [_dmstats_counters[c] for c in self.counters]
        return ",".join(fields)

    def extract(self, data):
        """ Extract the sum of this `DmStatsCounters` object's counter
            selection from the report row data passed in `data` and return
            the value as an integer. The report must use `--noheadings` with
            the default (':') field delimiter, and with the field
            specification returned by this object's `fields()` method.
        """
        fields = data.split(":")
        region = int(fields[0])
        counters = map(int, fields[1:])
        return (region, sum(counters))


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
    """ Return a list of floating point values that label the value
        (horizontal) axis.
    """
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
    """ Format an integer value to a short floating point string with
        units in binary prefix.
    """
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


class Bin(object):
    """ A histogram bin.

        Stores the current counter value and the starting offset and
        width, and allows splitting the bin into two equal parts (with
        the width and count distributed between the two).
    """
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
        self.count = int(self.count / 2)
        self.width = int(self.width / 2)
        return Bin(self.start + self.width, self.width, count=self.count)


_min_size_fraction = 40

# IOHistogram class constants
RENDER_COUNTS = 1
RENDER_TOTALS = 2


class IOHistogram(object):
    """ An IO distribution histogram.

        Contains current and running total bins for each region in the
        histogram, and provides methods to update, and render the IO
        distributions, and to obtain values derived from them.
    """
    counter = 0  # READS_COUNT
    nr_bins = 0
    dev_size = 0
    device = None
    bounds = None  # current bounds set
    regions = None  # current region_ids
    bins = None
    totals = None
    region_map = dict()

    # Handle to communicate with dmstats
    _dms = None

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

            `nr_regions` is the number of regions (bins) to divide the
            device into.
        """
        if not self.dev_size:
            log_error("Found zero-sized device: %s" % device)
            return None

        step = self.dev_size / nr_regions
        step = int(step)
        bounds = [x for x in range(step, self.dev_size + 1, step)]
        log_verbose(bounds)
        return bounds

    def __init__(self, device, counter, bounds, adapt=True):
        """ Initialise an `IOHistogram` for `device`, using the sector
            boundaries listed in `bounds`, and bind it to the dmstats
            counter specified by `counter`.

            Boundary values are given as an upper bound on the current
            bin, with an implicit lower bound of 0 on the first. The
            final bin should span all remaining space on the device.

            If `adapt` is True, the `IOHistogram` created will attempt
            to adapt bin count and sizes to the observed IO volumes,
            splitting and merging bins as the bin counts reach user
            specified thresholds.

            Bounds offsets are in 512b sectors.

            For e.g.:

            ```
            bounds = [2048, 4096, 6144, 8192]
            ioh = IOHistogram("vg00/lvol0", bounds, READS_COUNT, adapt=False)
            ```

            Will create a non-adaptive histogram with four bins from
            0..1M, 1M..2M,2M..3M, and 3M..4M, tracking the READS_COUNT
            counter.
        """
        self._dms = DmStats(device)
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
        """ Initialise an IOHistogram with `inital_bins` bins, evenly
            spaced across the specified `device`.

            All offsets are in 512b sectors.

            The new IOHistogram is bound to the counter field specified by
            counter.

            If `adapt` is True, the `IOHistogram` created will attempt
            to adapt bin count and sizes to the observed IO volumes,
            splitting and merging bins as the bin counts reach user
            specified thresholds.

            For e.g.:

            ```
            ioh = IOHistogram("vg00/lvol0", READS_COUNT, 8, adapt=True);
            ```

            Will create an adaptive histogram with eight bins and binds
            it to the READS_COUNT counter.
        """
        self._dms = DmStats(device)
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

    def update(self, test=False):
        """ Populate or update the histogram using the string counter
            values in `data`. The current value of the histogram is set
            to the counter values in `data` and this is added to the
            historical totals for the histogram.
        """
        if test:
            return self._test_update()

        data = self._dms.report()

        log_verbose("Updating %d histogram bins" % self.nr_bins)
        # data contains one row per histogram bin
        for line in data.splitlines():
            # python3 requires that both operands to String.split() be either
            # a byte sequence, or a string (but not one or the other).
            line = line.decode('utf8')
            fields = line.split(":")
            counters = fields[1:]
            region = int(fields[0])

            # Skip any non-dmioscope region_id values in the report.
            if region not in self.region_map.keys():
                log_verbose("Skipping foreign region %d" % region)
                continue

            log_verbose("%d: %s" % (region, counters[self.counter]))
            _bin = self.region_map[region]
            value = int(counters[self.counter])
            self.bins[_bin].count = value
            self.totals[_bin].count += value

    def _test_update(self):
        """ Generate simple test data to update an `IOHistogram`.
        """
        data = ""
        for i in range(self.nr_bins):
            data += "%d:%d\n" % (i, i)
        self.update(data)

    def min_width(self):
        """ Return width of the smallest bin in this `IOHistogram`.
        """
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

            If `total` is `True`, return the maximum value from the
            running total histogram instead of the current histogram.
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

            If `total` is `True`, return the maximum value from the
            running total histogram instead of the current histogram.
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

            If `total` is `True`, return the maximum value from the
            running total histogram instead of the current histogram.
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
        """ Print an ASCII representation of an `IOHistogram` and its
            values, suitable for display on a terminal of at least
            'colums' width.

            `render` specifies which values to display. It can include the
            current histogram counts, RENDER_COUNTS, or RENDER_TOTALS.

            FIXME: Plotting both RENDER_COUNTS and RENDER_TOTALS is not
            recommended since values are not normalized: the values of
            the totals distribution swamp the current interval.
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
            for step in range(int(scale)):
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
        """ Update the map of region_id values to `IOHistogram` bins.
        """
        index = range(len(self.regions))
        self.region_map = dict(zip(self.regions, index))

    def update_bin_regions(self, merge=False):
        """ Update `IOHistogram` bin widths and counts by adapting bins to
            current IO levels. Bins with counts falling above the current
            `_threshold` will be split, and their counts and widths
            distributed. If `merge` is `True` then adjacent bins with both
            having a count less than `_merge_threshold` will be merged
            into a single bin.
        """
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

        out = self._dms.list()
        log_verbose(out)

    def _create_bin_region(self, start, length):
        """ Call `dmstats` to create a new region on the bound device,
            beginning at `start`, and `length` sectors in size.

            Returns the region_id of the newly created region on success,
            or raises DmstatsException on error.
        """
        return self._dms.create(start, length)

    def create_bin_regions(self):
        """ Call `dmstats` to create regions on the bound device that
            correspond to the configured bin boundaries for this
            `IOHistogram`.

            Returns nothing on success and raises `DmstatsException` on
            error.
        """
        start = 0
        regions = []
        for bound in self.bounds:
            regions.append(self._create_bin_region(start, bound - start))
            start = bound
        self.regions = regions
        self.update_region_map()

    def _remove_bin_region(self, region_id):
        """ Call `dmstats` to delete the specified `region_id` from the
            bound device.
        """
        return self._dms.delete(region_id)

    def remove_bin_regions(self):
        """ Remove all regions managed by this `IOHistogram` from the
            bound device.
        """
        for region in self.regions:
            self._remove_bin_region(region)
        self.regions = []


_log_commands = True


# threshold at which to split a bin in two.
_threshold = 5000
# threshold at which to merge two bins together.
_merge_threshold = 0


def _parse_options(args):
    """ Parse `dmioscope` command line arguments.
    """
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
_histograms = {}


def _remove_all_regions():
    """ Remove all regions for all `IOHistogram` objects in `_histograms`.
    """
    for dev in _devices:
        _histograms[dev].remove_bin_regions()


def main(argv):
    """ Main `dmioscope` routine.
    """
    global _devices, _histograms, _merge_threshold, _threshold, _verbose

    args = _parse_options(argv)

    interval = args.interval

    if interval < 1:
        log_error("interval must be equal to or greater than 1.")
        return 1

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
        ioh = IOHistogram(dev, READS_COUNT, initial_bins=nr_bins, adapt=adapt)
        ioh.create_bin_regions()
        _histograms[dev] = ioh

    start_time = time.time()

    while count:
        # Sleep at the top of the interval to accumulate some initial
        # data to display. Since there's no direct way to access the
        # timerfd_create() or setitimer() interfaces from python use
        # a simple sleep() but correct the actual duration for each
        # interval based on the accumulated error.
        time.sleep(interval - ((time.time() - start_time) % interval))

        if clear:
            print(CLEAR_SCREEN)

        for dev in _devices:
            ioh = _histograms[dev]
            ioh.update(test=_test)

            if args.summary:
                print("%s: accumulated IO distribution" % dev)
                ioh.print_histogram(columns=width, render=RENDER_TOTALS)
            if args.current:
                print("%s: current IO distribution" % dev)
                ioh.print_histogram(columns=width)

            ioh.update_bin_regions(merge=merge)

            count -= 1

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
