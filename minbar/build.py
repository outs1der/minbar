"""
Functions to analyse data and add to MINBAR

def connect_db():
def scrape_jeanslist():
def match_name(name, name2, name3=None, comments=None):
def write_wiki_csv(outfile='minbar_sources_new.csv', con=None):
def verify_sources(con):
def update_ref(ref, minbar_id, ref_id=None, execute=True, commit=False):
def generate_ref_list(minbar_bursts):
def get_new(instr, ignore_unmatched=False, sources=None):
def analyse(instr, obs):
def burst_verify(candidates):
def burst_search(instr, src, obsid):
def update(instr, obs):
def write_minbar(minbar_obj=None, source='db', savefile=None):

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

from astropy.io import fits
import glob
from astropy.coordinates import SkyCoord
import astropy.units as u

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

    # URL updated 2025 Jun 4
    # url = 'http://www.sron.nl/~jeanz/bursterlist.html'
    url = 'https://www.sronpersonalpages.nl/~jeanz/bursterlist.html'
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

    m = np.array(m)
    missing = m == None
    if len(np.where(missing)[0]) > 0:
        minbar.logger.info('{} source(s) missing from my list:'.format(len(np.where(missing)[0])))
        # print(np.array(name)[missing])
        print(jeans_table['name'][missing])
    else:
        minbar.logger.info("all sources present and accounted for")

    # Check for dupes here

    check = set(m[~missing])
    if len(check) < len(m[~missing]):
        minbar.logger.error("duplicate name matches against Jean's list")

    # Should also check the completeness of the transient flags, and the
    # orbital periods
    # this is only done for the sources in Jean's list that are not
    # missing (and vice versa); careful with the indexing! (both objects
    # are pandas tables)

    minbar.logger.info('Got {} possible transients, and {} sources with orbital periods'.format(
        len(np.where(np.array(jeans_table['transient'][~missing]) != '')[0]),
        len(np.where(np.array(jeans_table['Porb'][~missing]) != '')[0])))

    mismatch = np.where((np.array(jeans_table['transient'][~missing]) != '')&
                        (~np.array(['T' in x for x in s['Type'].iloc[m[~missing]].values])))[0]
    if len(mismatch) > 0:
        print('\n{:<50} {:<7} {:<7}'.format('Source', 'Jean', 'Type'))
    for i in mismatch:
        print('{:<50} {:<7} {:<7}'.format(s['NAME'].iloc[m[~missing][i]]+' = '+jeans_table['name'][~missing].iloc[i], jeans_table['transient'][~missing].iloc[i], s['Type'].iloc[m[~missing][i]]))

    print('\n{:<50} {:<7} {:<7}'.format('Source', 'Jean', 'Porb'))

    for i in np.where((np.array(jeans_table['Porb'][~missing]) != ''))[0]:
        discrepant = True
        try:
            discrepant = np.abs(float(jeans_table['Porb'][~missing].iloc[i]) - s['Porb'].iloc[m[~missing][i]]) / s['Porb'].iloc[m[~missing][i]] > 0.05
        except:
            pass
        if discrepant:
            print('{:<50} {:<7} {:<7}'.format(s['NAME'].iloc[m[~missing][i]]+' = '+jeans_table['name'][~missing].iloc[i], jeans_table['Porb'][~missing].iloc[i], s['Porb'].iloc[m[~missing][i]]))

    return new


def update_ref(ref, minbar_id, ref_id=None, execute=True, commit=False):
    """
    This routine is to update the burst reference table with new entries.
    For example, say Smith et al. (2025) discuss the properties of burst
    2216 (but just that burst, so they don't give it a label or anything),
    in which case you'd want to add this reference as

    build.update_ref('smith25', 2216, '', commit=True)

    :param ref: bibliographic string, e.g. alizai23 (last name & 2-digit year, as listed in all.bib)
    :param minbar_id: list of MINBAR IDs to which the bursts correspond
    :param ref_id: array of IDs as in the reference; could be strings or letters, or nothing
    :param execute: set to True to actually execute the INSERT command; False for testing
    :param commit: set to True to commit

    :return: the connection object IF commit is False, so that you can review and commit later if need be
    """

    if ref_id is None:
        ref_id = np.full(len(minbar_id), '')
    if len(ref_id) != len(minbar_id):
        minbar.logger.error('ref_id needs to be same length as minbar_id (use blanks if missing)')
        return

    assert all(isinstance(x, (int, np.int64)) for x in minbar_id)
    assert (type(ref) == str) & (np.shape(ref) == ())

    # check for existing entries
    con = connect_db()
    existing = pd.read_sql_query("SELECT * from burst_ref WHERE ref_code='{}'".format(ref), con)
    if len(existing) > 0:
        for ex_entry in existing['entry_minbar']:
            if ex_entry in minbar_id:
                minbar.logger.error('ref {} for entry {} already present, #{}'.format(ref, ex_entry,
                        existing[existing['entry_minbar'] == ex_entry]['entry'].values[0]))
                return

    # maybe I don't have to generate the index here
    cur = con.cursor()
    cur.execute('PRAGMA foreign_keys = ON')  # always do this!
    # _command = "INSERT INTO burst_ref (ref_code, entry_minbar, ref_id, entry) VALUES (?, ?, ?, ?)"
    _command = "INSERT INTO burst_ref (ref_code, entry_minbar, ref_id) VALUES (?, ?, ?)"
    if not execute:
        print('Dry run:')
    for i, entry in enumerate(minbar_id):
        if execute:
            cur.execute(_command, (ref, entry, ref_id[i]))#, ind[i]))
        else:
            print (_command, (ref, entry, ref_id[i]))#, ind[i]))

    if commit == False:
        minbar.logger.info('changes not committed!')
        return con
    else:
        con.commit()
        con.close()


def generate_ref_list(minbar_bursts):
    """
    This function goes through the list of bursts, and generates new
    (numbered) reference lists for each burst, along with a running list
    of references

    :param minbar_bursts: MINBAR Bursts object, with appropriate sorting

    :returns: new refs column for all the bursts, and the ordered list of reference biblabels
    """

    if minbar_bursts.order_field != ['ra_obj','time']:
        minbar.logger.error('Bursts object needs to be sorted by ra_obj, time for the ordered output!')
        return None, None

    minbar.logger.info('generating new reference list...')

    con = connect_db()
    # read in the entire burst ref table here

    refs = pd.read_sql_query("SELECT * FROM burst_ref", con)

    ref_col = []
    ord_refs = [] # running array for refs
    for i, entry in enumerate(minbar_bursts['entry']):
        n_refs = [] # array collecting refs for this burst
        match = refs['entry_minbar'] == entry
        if np.any(match):
            for ref in refs[match]['ref_code']:
                if ref.strip() not in ord_refs:
                    ord_refs.append(ref.strip())
                n_refs.append(ord_refs.index(ref.strip())+1)

            # print (entry, b['refs'][i], ''.join(str(n_refs)[1:-1].split(' ')))
            ref_col.append(''.join(str(n_refs)[1:-1].split(' ')))
        else:
            ref_col.append('-')
            if minbar_bursts['refs'][i] != '-':
                # no matches, so ref should be null
                minbar.logger.warning('non-null refs in table but no refs found in db for burst {}'.format(entry))

    minbar.logger.info('... new reference list complete')

    return ref_col, ord_refs


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

    # Load in the existing observations, so we can skip the ones already ingested
    # Need to fix this hard-wired path prior to deployment, or allow the option to use
    # the text file instead

    obs_db = minbar.Observations() # '../minbar-obs')

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

            # observations = os.scandir(instr.path + '/data/' + instr.source_path[i])
            observations = os.scandir('/'.join([minbar.MINBAR_ROOT, instr.path, 'data', instr.source_path[i]]))
            # in_dir = set([str.encode(x.name) for x in observations])

            # need to filter here potentially for matching observation directories, and omit
            # other stuff (e.g. .DS_Store!)

            in_dir = set([x.name for x in observations])

            new = in_dir.difference(in_db)
            if len(new) > 0:
                to_add[source] = new
        else:
            print("Skipping {}, no matching source directory".format(source))

    return to_add


def get_obs_info(test_dir, target_dir, debug=False):
    """Here we want to gather basic information about the observation data contained in the directory of interest:
    the source name(s), pointing coordinates, time range etc.
    minbar-obs columns: name, instr, obsid, tstart, tstop, [ra_obj, dec_obj], angle
    
    :param test_dir: directory name to scan
    :param path: absolute path in which the test_dir resides
    :param debug: set to True to provide debugging info
    
    :returns: (MINBAR) name, instr flag (may be preliminary), obsid, MJD start, stop, pointing RA, Dec (degrees), target separation (arcmin)
    """

    def get_fmi(path):
        """
        Function to return the basic FMI information from an RXTE directory

        :param path: path to a (presumed) RXTE data directory, for example downloaded from the archive

        :returns: tuple with target name, RA & Dec, ObsID, start & stop times (MJD)
        """

        fmi_file = '/'.join((path, 'FMI'))
        if os.path.isfile(fmi_file):
            hdul = fits.open(fmi_file)

            header = hdul[0].header
            # print (header['MJDREFI'], header['MJDREFF'])
            # hdul.info()
            # print (hdul[1].data)#.data[0])
            _name, _ra_obj, _dec_obj = hdul[1].data[0]['Source'], hdul[1].data[0]['RA'], hdul[1].data[0]['Dec']
            _obsid, _metstart, _metstop = hdul[1].data[0]['ObsId'], hdul[1].data[0]['StartMET'], hdul[1].data[0]['StopMET']
            # print(np.array(hdul[1].data[0])[[0,1,2,10,11,12]])#.header
            hdul.close()

            _tstart = minbar.pca_time_to_mjd_tt(_metstart*u.s, header['MJDREFI'], header['MJDREFF'])
            _tstop = minbar.pca_time_to_mjd_tt(_metstop*u.s, header['MJDREFI'], header['MJDREFF'])

            return _name, _ra_obj, _dec_obj, _obsid, _tstart, _tstop

        else:
            return None

    if debug:
        print ('scanning directory {}, subdirectory {}'.format(target_dir, test_dir))
        
    try:
        # test for RXTE INDEX file here
        _name_FMI, _ra_obj, _dec_obj, _obsid, _tstart, _tstop = get_fmi('/'.join((target_dir, test_dir)))
        # print(np.array(hdul[1].data[0])[[0,1,2,10,11,12]])#.header
            
        if debug:
            print ("Got XTE FMI file, obsid {}, target {} (RA={:.5f}, Dec={:.5f}), time range MJD {:.5f}-{:.5f}".format(
                _obsid, _name_FMI, _ra_obj, _dec_obj, _tstart.value, _tstop.value))

        if (len(_obsid) == 15):
            if _obsid[-1] in 'ASZ':
                minbar.logger.warning ('obsid suffix indicates slew/scan data, suggest ignoring')
            elif _obsid[-1] in '0123456789':
                minbar.logger.warning ('part n>1 of a multi-part observation, find the root (presumably {})'.format(
                    _obsid[:-1]))
        else:
            # check for multi-parts
            # print ('/'.join((target_dir, test_dir))+'?')
            additional = glob.glob('/'.join((target_dir, test_dir))+'?')
            if additional != []:
                if debug:
                    print ('Got additional obsID directories: {}'.format(additional))
                _obsid = [_obsid]
                for _dir in additional:
                    _obsid.append(_dir.split('/')[-1:][0])
                    _name2, _ra_obj2, _dec_obj2, _obsid2, _tstart2, _tstop2 = get_fmi(_dir)
                    assert _name2 == _name_FMI
                    _tstart = min([_tstart, _tstart2])
                    _tstop = max([_tstop, _tstop2])
                    # print ("Got XTE FMI file, obsid {}, target {} (RA={}, Dec={}), time range MET {}-{}".format(
                    # _obsid2, _name2, _ra_obj2, _dec_obj2, _tstart2, _tstop2))

            # perform cone search for source in MINBAR.Source, make sure of an unambiguous ID
            s = minbar.Sources()
            c = SkyCoord(ra=s['ra_obj']*u.degree, dec=s['dec_obj']*u.degree, frame='icrs')
            t = SkyCoord(ra=_ra_obj*u.degree, dec=_dec_obj*u.degree, frame='icrs')
            sep = t.separation(c)

            in_fov = sep < 1*u.degree
            if debug:
                print ('measured separation vs. MINBAR sources, got {} within FoV'.format(sum(in_fov)))
            
            # _name_FMI = _name#.copy()
            
            if sum(in_fov) == 0:
                minbar.logger.error ('no MINBAR sources in the FoV, have you got the right observation?')
                return None
            elif sum(in_fov) == 1:
                _name = s['name'][in_fov]
            else:
                _name = s['name'][np.argmin(sep)]
                _sep = min(sep)
                minbar.logger.warning ('more than one source in the FoV, selecting the closest to the aimpoint ({} =? {}, offset by {:.3f})'.format(_name, _name_FMI, min(sep)))

            # set the instrument name, we can't yet know the full label
            _instr = 'PCA'
            # RXTE section complete

        return _name, _instr, _obsid[0], _tstart.value, _tstop.value, _ra_obj, _dec_obj, _sep.to('arcmin').value
    
    except:
        # something went wrong, so presumably this is not an RXTE observation
        minbar.logger.warning ('information gathering for RXTE observation failed, falling back on... nothing for now')
        pass


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


def write_minbar(minbar_obj=None, source='db', savefile=None):
    """
    Method to write a table to a text file, replicating the original
process by which the DR1 tables were created

    :param minbar_obj: MINBAR object to write; if omitted, will write all 3 objects
    :param source: source of data; defaults to 'db'
    :param savefile: name of savefile to write to; if omitted, will write to standard out
    """

    if minbar_obj is None:
        # if you haven't supplied an object, we need to get it from
        # whatever source is specified

        if not (source in ['db']):
            minbar.logger.error('source {} not currently supported'.format(source))
            return None

        if savefile is not None:
            # user needs to supply three savefiles 
            if len(savefile) != 3:
                minbar.logger.error('need to supply a savefile name for each of the three objects: Bursts, Observations and Sources')
                return None
        else:
            savefile=['minbar_dr2.txt', 'minbar-obs_dr2.txt', 'minbar-sources_dr2.txt']

        minbar.logger.info('reading MINBAR data from {}'.format(source))
        o = minbar.Observations(source)
        s = minbar.Sources(source)

        write_minbar(o.bursts, savefile=savefile[0])
        write_minbar(o, savefile=savefile[1])
        write_minbar(s, savefile=savefile[2])

        return

    # now actually do the write part

    if type(minbar_obj) == minbar.Bursts:
        # as of 2025 Nov, this produces a file that (up to the header and
        # whitespace) is identical to the DR1 MRT file. But need more work to
        # future-proof it, including
        # - a way to update the burst references

        if savefile is not None:
            minbar.logger.info('writing object as MINBAR bursts table to file {}'.format(savefile))

        # it is assumed you have already re-generated the burst reference
        # list prior to this step; see build.generate_ref_list
        ascii.write(minbar_obj.records, output=savefile, overwrite=True, 
            # can use this to fill the blank numbers, but it doesn't respect the format
            format='mrt', fill_values=[(ascii.masked, 0.0)], 
            # formats could be part of the table
            # copied from the DR1 MRT file
            formats={'name': '<23', 'instr': '<3', 'obsid': '<15',
                'time': '11.5f', 'entry': '4d', 'entry_obs': '6d', 'bnum': '4d',
                'xref': '4d', 'mult': '1d', 'angle': '7.2f', 'vigcorr': '5.3f',
                'sflag': '<11', 'rexp': '5.2f', 'rise': '5.2f', 'tau': '5.1f',
                'e_tau': '5.1f', 'dur': '6.1f', 'e_dur': '6.1f', 'edt': '6.1f',
                'e_edt': '7.3f', 'tdel': '8.1f', 'trec': '7.1f',
                'perflx': '6.3f', 'e_perflx': '6.3f', 'alpha': '6.1f',
                'e_alpha': '6.1f', 'bc': '5.3f', 'e_bc': '5.3f',
                'gamma': '6.4f', 'sc': '6.3f', 'hc': '6.3f', 's_z': '6.3f',
                'pflux': '6.2f', 'e_pflux': '5.2f', 'fluen': '8.3f',
                'e_fluen': '7.3f', 'bpflux': '6.2f', 'e_bpflux': '5.2f',
                'kT': '4.2f', 'e_kT': '4.2f', 'rad': '6.1f', 'e_rad': '5.1f',
                'bfluen': '6.4f', 'e_bfluen': '6.4f', 'refs': '<20'})

    elif type(minbar_obj) == minbar.Observations:
        minbar.logger.info('writing object as MINBAR observations table')
        minbar.logger.info('(not yet implemented)')
        return

    elif type(minbar_obj) == minbar.Sources:
        minbar.logger.info('writing object as MINBAR sources table')
        minbar.logger.info('(not yet implemented)')
        return

    else:
        minbar.logger.error('unknown type for object {}'.format(minbar))

    return


print('loaded build tools')
