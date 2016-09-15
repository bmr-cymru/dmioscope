# dmioscope - track IO distribution with dmstats

The dmioscope script provides a simple tool to track and visualise
the distribution of IO requests across a device's address space by
maintaining a histogram that counts IO requests by location.

The device can either be split into a fixed number of equally
sized regions (the same set of regions will be reported in each
successive update), or the script can attempt to adapt the number,
and size, of bins to the observed IO volumes: in this case the
set of bins (and their sizes) may change between updates.

The static mode is useful in order to obtain fixed, repeatable
statistics for measurement purposes: adaptive mode is best suited
to interactive use where the user is trying to get a rough idea
of the current IO distribution on the device.

The command accepts count and interval arguments to allow reports
to be continually generated from current data.


# Requirements

To use dmioscope you need a Linux installation with one or more
device-mapper devices to monitor, and an installed copy of the
device-mapper tools and libraries that includes the 'dmstats'
command (at least device-mapper 1.02.104 - included in the lvm2
v2.02.127 tarball and any later distribution packages).

# IOScope histograms

The dmioscope command displays graphical representations of the quantity
of IO reaching different regions of a device-mapper device. This allows
an estimate of the IO distribution to be shown for the device, both for
short and long-term observations.

## Histogram bins

Each bin in an ioscope histogram corresponds to a fixed range of
device sectors, for example, consider an 8MiB device divided into
eight bins:

```
   Bin:  0    1    2    3    4    5    6    7
        +----+----+----+----+----+----+----+----+
 Count: |  0 | 10 |  0 | 50 | 15 | 0  | 90 | 77 |
        +----+----+----+----+----+----+----+----+
```

In this case, each bin is of equal size (1MiB), although this need
not be the case generally.

The count for each bin indicates the number of matching IO events
that were observed in that bin for the current interval.

When plotted graphically, the area of each bin's plotted region is
proportional to the count of IOs observed in that region:

```
        0.0           20.0        40.0        60.0        80.0       100.0
         +-------------------------------------------------------------
0.0B     |
1.1MiB   |########
2.1MiB   |
3.1MiB   |################################
4.1MiB   |##########
5.1MiB   |
6.1MiB   |#######################################################
7.1MiB   |##################################################
```

If all bin sizes are uniform then the height of each bar corresponds to
the number of IOs observed in that region during the last interval.

If bin sizes are non-uniform then the height of the bar represents the
*frequency* of IO observed in that region: the absolute count of IOs is
proportional to the total area of the bar.

For example, a 300MiB region and a 100MiB region observing 300 IOs and
100 IOs respectively will be plotted as:

```
        0.0          100.0       200.0       300.0       400.0       500.0
         +------------------------------------------------------------
0.0B     |##############
         |##############
         |##############
300.1MiB |##############
```

Although the 300MiB region has three times the absolute number of IOs,
it is also three times larger and hence has the same frequency of IOs
per unit disk space:

  100 IOs / 100MiB

# Running ioscope

The command is currently a single python script: there is no need to
install the file, or make it executable, in order to use it (although
adding exec permissions to the file will allow it to be run without
explicitly calling the python interpreter).

The command can be run either as:

```shell
#program python ./dmioscope.py vg00/lvol0
```

Or by first making it executable with `chmod`:

```
# chmod +x dmioscope.py
# ./dmioscope.py vg00/lvol0
```


## Command line options and arguments
=====================================

```
usage: dmioscope.py [-h] [-a] [-b NR_BINS] [--clear] [-c] [-m]
                    [-M MERGE_THRESH] [-n] [-p] [-r ROWS] [-s]
                    [-t THRESH] [-v] [-w WIDTH]
                    [interval] [count] dev [dev ...]

IOScope arguments.

positional arguments:
  interval              Duration of the report interval in seconds.
  count                 Number of intervals to report.
  dev                   Device(s) to monitor.

optional arguments:
  -h, --help            show this help message and exit
  -a, --adaptive        Adapt the number of bins to observed IO volume.
  -b NR_BINS, --bins NR_BINS
                        Divide devices into nr equally sized bins.
  --clear               Clear the screen before each update.
  -c, --current         Show the current interval plot.
  -m, --merge           Allow merging of bins with low IO volume in adaptive
                        mode.
  -M MERGE_THRESH, --merge-threshold MERGE_THRESH
                        Threshold at which to merge adjacent regions with low
                        IO.
  -n, --no-adaptive     Do not adapt the number of bins according to observed
                        IO volume.
  -p, --percent         Show distribution values as percentages.
  -r ROWS, --rows ROWS  Specify the maxumum number of rows to use.
  -s, --summary         Show the accumulated summary plot.
  -t THRESH, --threshold THRESH
                        Threshold at which to split a bin when using
                        --adaptive
  -v, --verbose         Enable verbose debugging messages.
  -w WIDTH, --width WIDTH
                        Specify the maximum terminal width to use.
```

