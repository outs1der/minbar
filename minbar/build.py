"""
Classes to analyse data and add to MINBAR

Prior to running this part of the code, it is assumed that the user has assembled a set of data
in a data directory, organised by source and observation ID. E.g. for JEM-X, the data directory
would be organised by (translated) source name -> obs ID:
IGRJ00291p5934 -> 108300510010
                  108300530010
                        .
                        .
                        .
Within each observation directory, the minimal requirements are for a lightcurve and observation-
averaged spectrum

Below is an example set of commands to run the analysis (for JEM-X)

import minbar.build
jemx = minbar.build.Instrument('JEM-X', '../jemx', 'IJ')
r = minbar.build.get_new(jemx, ignore_unmatched=True)
test = { 'EXO 0748-676': set(list(r['EXO 0748-676'])[:2]) }
# test = {'EXO 0748-676': {'117200690010', '117500700010'}}
minbar.build.analyse(jemx, test)

You can also add any other instrument, provided you've created the directory tree, e.g.

xmm = minbar.build.Instrument('XMM-Newton', '../xmm', 'XN')
"""

import os
import sys
import minbar
import numpy as np

# First group of functions are intended to work on observations and candidate events from
# any instrument; need to pass only enough information e.g. lightcurve path(s), candidate
# times, significance etc. to allow verification

def verify_path(source, path, source_path):
    """
    Generic routine to verify the data path, and try to match up the source names with the
    directories in the tree
    :param path:
    :param source_path:
    :return:
    """

    has_dir = [False] * len(source)
    nonmatched = 0
    with os.scandir(path + '/data/') as entries:
        for entry in entries:
            if os.path.isdir(entry) and not (entry.name in source_path):
                print("** WARNING ** possible non-compliant source path {}".format(entry.name))
                nonmatched += 1
            else:
                # self.has_dir[self.source_path.index(entry.name)] = True
                _i = np.where(source_path == entry.name)[0]
                # print (entry.name,_i)
                assert len(_i) <= 1
                if len(_i) == 1:
                    has_dir[_i[0]] = True

    nmissing = len(np.where(np.logical_not(has_dir))[0])
    if nonmatched > 0:  # or self.nmissing > 0:
        print("""
Possible inconsistencies in the directory tree; 
  {} directories without source matches
  {} sources without data directories
You should check your source to directory name mapping""".format(nonmatched, nmissing))

    return has_dir, nonmatched, nmissing


def get_new(instr, ignore_unmatched=False):
    """
    Method to identify new observations to be analysed. Having defined an instrument, run
    newobs = instrument.get_new()
    Returns a dictionary with source name keys and a set of new obsids for each
    :param instr:
    :return:
    """

    if instr.nonmatched > 0 and ignore_unmatched == False:
        print("** ERROR ** can't continue with unmatched observation directories")
        return

    # Need to fix this hard-wired path prior to deployment, or allow the option to use
    # the text file instead

    obs_db = minbar.Observations('../minbar-obs')

    to_add = {}
    for i, name in enumerate(instr.source_name):
        if instr.has_dir[i]:
            print(name, instr.source_path[i])

            # Filter obs DB on source

            obs_db.clear()
            obs_db.select_like(name)
            # in_db = set([x.strip() for x in obs_db['obsid']])
            in_db = set([x.strip().decode() for x in obs_db['obsid']])

            # Check if the observation is already present

            observations = os.scandir(instr.path + '/data/' + instr.source_path[i])
            # in_dir = set([str.encode(x.name) for x in observations])
            in_dir = set([x.name for x in observations])

            new = in_dir.difference(in_db)
            if len(new) > 0:
                to_add[name] = new
        else:
            print("Skipping {}, no matching source directory".format(name))

    return to_add


def analyse(instr, obs):
    """
    Method to analyse new observations. Example:
    test = { 'EXO 0748-676': set(list(r['EXO 0748-676'])[:2]) }
    test = {'EXO 0748-676': {b'117200690010', b'117500700010'}}
    minbar.build.analyse(jemx, test)
    :param instr:
    :param obs:
    :return:
    """

    # Create blank result tables here and populate them one by one

    pers, bursts = [], []
    for src in obs.keys():
        for obsid in obs[src]:
            print(src, obsid)

            # This routine might be also brought out to be instrument-independent

            pers = instr.analyse_persistent(src, obsid)

            candidates = burst_search(instr, src, obsid)

            confirmed = ()
            if candidates is not None:
                confirmed = burst_verify(candidates)

            for burst in confirmed:
                tmp = instr.analyse_burst(src, obsid, burst)

    return pers, bursts


def burst_verify(candidates):
    """
    Instrument-agnostic function to verify (automatically or by input from the user) burst
    candidates
    :param instr:
    :param candidates:
    :return:
    """

    # Prelimimary behaviour is to accept all candidates

    return candidates


def burst_search(instr, src, obsid):
    """
    Function to search through a lightcurve and identify candidate bursts
    :param instr:
    :param src:
    :param obsid:
    :return:
    """

    # Preliminary behaviour is to fail to find bursts

    return None


def update(instr, obs):
    """
    Method to reanalyse existing observations
    :return:
    """
    pass


# Here's a generic instrument class which can be adapted and/or duplicated for different
# instruments. This class is kept pretty lean to avoid having to replicate lots of code

class Instrument:
    """
    Defines the properties of an instrument with data that we're going to analyse and add to MINBAR
    Example
    import minbar.build
    jemx = minbar.build.Instrument('JEM-X', '../jemx', 'IJ')
    """

    INSTRUMENTS =['XP','SW','IJ']

    def __init__(self, name, path, label,
                 lightcurve=['lc1_3-25keV_1.0s.fits','lc2_3-25keV_1.0s.fits'],
                 spectrum=None):

        self.name = name
        if os.path.isdir(path):
            self.path = path
        else:
            sys.exit(1)
        self.label = label
        if label not in self.INSTRUMENTS:
            print('** WARNING ** new instrument {} not fully implemented'.format(name))

        # Define the paths corresponding to each source here
        # Default is just the source name with spaces removed; this won't work for RXTE

        self.source = minbar.Sources()
        self.source_name = self.source['name']
        self.source_path = self.source_name.replace(" ", "")
        if self.label == 'IJ':
            # for JEM-X, the convention is to also replace '+' with 'p'
            self.source_path = self.source_name.replace("+","p")

        # Define the lightcurve and spectral files

        self.lightcurve = lightcurve
        self.spectrum = spectrum

        # Check that the directories in the data path all correspond to source_paths
        # You don't want to miss observations in some directory because of inconsistencies
        # with the directory names

        self.has_dir, self.nonmatched, self.nmissing = verify_path(self.source, self.path, self.source_path)


    def __str__(self):
        """
        Method to display information about this object
        :return:
        """

        return """
MINBAR instrument definition

Name: {} ({})
Data path: {}
Lightcurve(s): {}
Spectra: {}""".format(self.name, self.label, self.path, self.lightcurve, self.spectrum)


    def analyse_persistent(self, src, obsid):
        """
        Function to analyse the persistent emission (lightcurve and spectrum) for a single
        observation
        :param src:
        :param obsid:
        :return:
        """

    def analyse_burst(self, bursts):
        """
        Function to analyse the persistent emission (lightcurve and spectrum) for a single
        observation
        :param src:
        :param obsid:
        :return:
        """



