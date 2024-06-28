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
import sqlite3
import pandas as pd
import logging
import requests
import re
import bs4

def connect_db():
    """
    Set up a connection to the (new) database

    :return: database connection
    """

    return sqlite3.connect(minbar.MINBAR_DR2+"/minbar.db")

# Set of functions and tools to operate on the (new) source database

def scrape_jeanslist():
    """
    This routine will scrape Jean's list, and return the table for comparison

    following this advice https://www.storybench.org/how-to-scrape-a-table-from-a-website-with-python/

    :return: pandas table with all data
    """

    url = 'http://www.sron.nl/~jeanz/bursterlist.html'
    req = requests.get(url)

    if req.status_code != 200:
        minbar.logger.error('can''t access {}, check connection'.format(url))
    else:
        # now parse the code
        soup = bs4.BeautifulSoup(req.text, 'html.parser')

        # columns are: Object, Transient?, Orbital period (hr), Accretion rate (% of Edd) and range,
        # Oscillation/pulsation (Hz), Superburst?, Burst only?**, mHz QPOs
        table = soup.find_all('table')[1]
        rows = table.find_all('tr')
        name, name2, host, transient, Porb = [], [], [], [], []
        for row in rows[1:]:
            # format for name includes possible aliases and hosts, for example
            # IGR J18245-2452=AX J18245-2451? (M28)
            cells = row.find_all('td')
            transient.append(cells[1].text.strip())
            Porb.append(cells[2].text.strip())
            _name = cells[0].text.strip()
            x = re.split("[()]", _name)
            if len(x) > 1:
                host.append(x[1].strip())
            else:
                host.append('')
            y = re.split("[=\/]", x[0])
            if len(y) > 1:
                name2.append(y[1].strip())
                name.append(y[0].strip())
            else:
                name2.append('')
                name.append(x[0].strip())
            # print (name[-1],',', name2[-1],',', host[-1])

        minbar.logger.info('got {} burst sources from {}'.format(len(name),url))

        table_data = {'name': name, 'name_alt': name2, 'host': host, 'transient': transient, 'Porb': Porb}

        return pd.DataFrame(table_data)

# now check against my data. First match the names.

def match_name(name, name2, name3=None, comments=None):
    '''
    Function to try to match Jean's names, in my list, or in the alt list, or in
    the comments. Latter two inputs are optional

    :param name: (single) name to match, from Jean's table
    :param name2: MINBAR source names
    :param name3: alternate name listed in MINBAR table (Name_2)
    :param comments: COMMENTS field of MINBAR table, which includes alternate names

    :return: index of match into MINBAR source table
    '''

    assert name != ''
    assert len(name2) > 0

    # some specials here

    name = 'Swift'.join(name.split('SWIFT'))
    if name == '4U 1822-00':
        name += '0'
    elif name == 'AX J18245-2451':
        name = 'AX J1824.5-2451'
    elif name == '4U 2129 + 11':
        name = '4U 2129+11'

    _i = np.where(name.strip() == name2)[0]  # might have problems with trailing blanks
    assert len(_i) <= 1
    if len(_i) > 0:
        return _i[0]

    # try the alternate names
    if name3 is not None:
        _i = np.where(name.strip() == name3)[0]  # might have problems with trailing blanks
        assert len(_i) <= 1
        if len(_i) > 0:
            return _i[0]

    # try the comments
    if comments is not None:
        _i = []
        for i, comment in enumerate(comments):
            # search as lower case to match "Rapid Burster" with "Rapid burster" (groan)
            if comment is not None:
                _pos = comment.lower().find(name.lower())
                if _pos >= 0:
                    _i.append(i)
        # print (name, _i)
        if (len(_i) > 1):
            if _i[0] != 39:
                print('  ---> problem! ', name, _i)

        if len(_i) > 0:
            return _i[0]

    return None

def write_wiki_csv(outfile='minbar_sources_new.csv', con=None):
    """
    This routine converts the source table into a form for use on the wiki, at
    https://burst.sci.monash.edu/wiki/index.php?n=MINBAR.SourceTable
    and saves to a file; defaulting to the DR2 directory

    :param filename: file to save the output as
    :param con: connection to the database, will be created if not supplied
    :return: None
    """

    new_connection = False
    if con is None:
        con = connect_db()
        new_connection = True

    include_columns = ['NAME', 'Type', 'RA_OBJ', 'DEC_OBJ', 'NH', 'Porb', 'nburst', 'exp', 'rate', 'disc_bibcode',
                       'NH_bibcode', 'Porb_bibcode']
    df = pd.read_sql_query("SELECT {} FROM sources ORDER BY RA_OBJ".format(",".join(include_columns)), con)

    # this from source_fits_to_csv.py
    # don't need the str.decode('utf-8') clauses anymore, I think
    # df['NAME'] = df.disc_bibcode.str.decode('utf-8').map('[[http://adsabs.harvard.edu/abs/{}|'.format).str.strip()+df.NAME.str.decode('utf-8').str.strip()+']]'
    df['NAME'] = df.disc_bibcode.map('[[http://adsabs.harvard.edu/abs/{}|'.format) + df.NAME + ']]'

    # nz = df.NH.notnull()
    nz = df.NH > 0.
    df.loc[nz, 'NH'] = df.NH_bibcode[nz].map('[[http://adsabs.harvard.edu/abs/{}|'.format) + df.NH[nz].map(
        '{:.2f}'.format) + ']]'
    df.loc[~nz, 'NH'] = ''

    nz = df.Porb > 0.
    df.loc[nz, 'Porb'] = df.Porb_bibcode[nz].map('[[http://adsabs.harvard.edu/abs/{}|'.format) + df.Porb[nz].map(
        '{:.9g}'.format) + ']]'
    df.loc[~nz, 'Porb'] = ''

    df.exp = (df.exp / 1e3).map('{:.1f}'.format)

    outfile = 'minbar_sources_new.csv'
    print('Writing table to {}...'.format(outfile))

    df[include_columns[:9]].to_csv(minbar.MINBAR_DR2+'/'+outfile, float_format='%.9g', index=False)
    print('...done.')

    if new_connection:
        con.close()

def verify_sources(con):
    """
    Check the source list and see if there are any updates required
    :param con: database connection

    :return: list of sources requiring updates
    """

    new = pd.read_sql_query("""SELECT * FROM sources 
        WHERE ((disc = '' OR pos_method = '') OR (type LIKE '%P%' AND Porb <= 0.0))""", con)
    if len(new) > 0:
        minbar.logger.warning('got {} new sources needing discovery information:\n'.format(len(new)))
        print (new[['NAME', 'Type', 'disc', 'pos_method', 'Porb']])
    else:
        minbar.logger.info('basic information seems to be present for all sources, nothing to do')

    jeans_table = scrape_jeanslist()

    s = pd.read_sql_query("SELECT * from sources", con)

    # Find all the elements in Jean's list that are NOT in mine
    m = []
    for i, _name in enumerate(jeans_table['name']):
        m.append(match_name(_name, s['NAME'].values, s['Name_2'].values, s['COMMENTS'].values))
        # print(_name, m[-1])
    # s.columns

    missing = np.array(m) == None
    if len(np.where(missing)[0]) > 0:
        minbar.logger.info('{} sources missing from my list:'.format(len(np.where(missing)[0])))
        print(np.array(name)[missing])
    else:
        minbar.logger.info("all sources present and accounted for")

    # Should also check the completeness of the transient flags, and the orbital periods

    print('Got {} possible transients, and {} sources with orbital periods'.format(
        len(np.where(np.array(jeans_table['transient']) != '')[0]),
        len(np.where(np.array(jeans_table['Porb']) != '')[0])))

    mismatch = np.where((np.array(jeans_table['transient']) != '') &
                        (~np.array(['T' in x for x in s['Type'][m].values])))[0]
    if len(mismatch) > 0:
        print('\n{:<50} {:<7} {:<7}'.format('Source', 'Jean', 'Type'))
    for i in mismatch:
        print('{:<50} {:<7} {:<7}'.format(s['NAME'][m[i]]+' = '+jeans_table['name'][i], jeans_table['transient'][i], s['Type'][m[i]]))

    print('\n{:<50} {:<7} {:<7}'.format('Source', 'Jean', 'Porb'))

    for i in np.where((np.array(jeans_table['Porb']) != ''))[0]:
        discrepant = True
        try:
            discrepant = np.abs(float(jeans_table['Porb'][i]) - s['Porb'][m[i]]) / s['Porb'][m[i]] > 0.05
        except:
            pass
        if discrepant:
            print('{:<50} {:<7} {:<7}'.format(s['NAME'][m[i]]+' = '+jeans_table['name'][i], jeans_table['Porb'][i], s['Porb'][m[i]]))

    return new

# This group of functions are intended to work on observations and candidate events from
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



