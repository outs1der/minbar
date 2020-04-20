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

def get_new(instr, ignore_unmatched=False, sources=None):
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

    # Build the list of source names here

    if sources is None:
        # If you don't pass any sources, we'll do them all
        sources=instr.source_name
    else:
        if isinstance(sources, str):
            sources = [sources]

        # Check the passed list of sources here
        for source in list(sources):
            if source not in instr.source_name:
                print("** ERROR ** {} not present in source name list".format(source))
                return None

    # Need to fix this hard-wired path prior to deployment, or allow the option to use
    # the text file instead

    obs_db = minbar.Observations('../minbar-obs')

    to_add = {}
    for source in sources:
        i = np.where(instr.source_name == source)[0][0]
        if instr.has_dir[i]:
            print(source, instr.source_path[i])

            # Filter obs DB on source

            obs_db.clear()
            obs_db.name_like(source)
            in_db = set([x.strip() for x in obs_db['obsid']])
            # in_db = set([x.strip().decode() for x in obs_db['obsid']])

            # Check if the observation is already present

            observations = os.scandir(instr.path + '/data/' + instr.source_path[i])
            # in_dir = set([str.encode(x.name) for x in observations])
            in_dir = set([x.name for x in observations])

            new = in_dir.difference(in_db)
            if len(new) > 0:
                to_add[source] = new
        else:
            print("Skipping {}, no matching source directory".format(source))

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

    print ('analysing...')
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



