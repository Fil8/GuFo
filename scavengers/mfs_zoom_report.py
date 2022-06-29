# Import modules

import numpy as np
import os
import sys
import json
import argparse
from glob import glob
from matplotlib import pyplot as plt
from datetime import datetime



# Define functions

def create_parser():
    p = argparse.ArgumentParser()
    p.add_argument("-ob", "--obs-dir", type=str, default=os.getcwd(),
                   help="Common parent directory of the single-observation directories"
                        " mfs??/mfs_report_track? . Default is CWD.")
    p.add_argument("-cj", "--create-json", action="store_true",
                   help="Create JSON dictionary with name given by the -jd"
                        " command-line option. Default is false.")
    p.add_argument("-jd", "--json-dict", type=str, default="",
                   help="Full path of JSON dictionary to be created (-cj) or plotted"
                        " (-pj). Default is mfs_zoom_report.json .")
    p.add_argument("-pj", "--plot-json", action="store_true",
                   help="Create plots from JSON dictionary with name given by the -jd"
                        " command-line option. Default is false.")
    p.add_argument("-pd", "--plot-dir", type=str, default=os.getcwd(),
                   help="Directory where to write the JSON dictionary and the report"
                        " plots. Default is CWD.")
    return p


def plot_json_dict(json_dict, plot_dir):
  if plot_dir[-1]=='/':
    plot_dir = plot_dir[:-1]

  if not json_dict:
    print('# ERROR: You forgot to give the full path of the JSON dictionary with the -jd command-line option.')
    sys.exit()
  elif not os.path.exists(json_dict):
    print('# ERROR: Cannot find JSON report {0:s}. Aborting ...'.format(json_dict))
    sys.exit()

  with open(json_dict) as jd:
    zoom_report = json.load(jd)
  print('# Generating plots for {0:d} observations based on JSON report {1:s}'.format(len(zoom_report),json_dict))

  t_start, t_end = [], []
  nr_ant, nr_ant_eff = [], []
  corr_flag_frac, cal_corr_flag_frac = [], []
  line_exp_noise, line_meas_noise = [], []
  scan_date, scan_hour, scan_corr_flag_frac = [], [], []
  cont_noise = []
  for obs in zoom_report:
    obs_rep = zoom_report[obs]
    cal_corr_flag_frac.append(obs_rep['CALIBRATORS_VIS']['CORR_FLAG_FRAC'])
    target_vis = obs_rep['TARGET_VIS']
    t_start.append(datetime.strptime(obs_rep['TIME_START'], '%d-%b-%Y/%H:%M:%S.%f'))
    t_end.append(datetime.strptime(obs_rep['TIME_END'], '%d-%b-%Y/%H:%M:%S.%f'))
    nr_ant.append(obs_rep['NR_ANT'])
    nr_ant_eff.append(obs_rep['NR_ANT'] - target_vis['SPW2']['NR_FLAG_ANT'])
    corr_flag_frac.append(target_vis['SPW2']['CORR_FLAG_FRAC'])
    line_exp_noise.append(obs_rep['TARGET_LINE']['EXP_NOISE_30'])
    line_meas_noise.append(obs_rep['TARGET_LINE']['MEAS_NOISE_30'])
    cont_noise.append(obs_rep['TARGET_CONT']['NOISE'])
    for vv in target_vis:
      if 'SCAN_' in vv:
        t0 = datetime.strptime(target_vis[vv]['SCAN_START'], '%H:%M:%S.%f')
        scan_date.append(t_start[-1])
        scan_hour.append(t0.hour)
        scan_corr_flag_frac.append(target_vis[vv]['CORR_FLAG_FRAC_SPW2'])

  scan_date = np.array(scan_date)
  scan_hour = np.array(scan_hour)
  scan_corr_flag_frac = np.array(scan_corr_flag_frac)
  for hh in range(24):
    select_h = (scan_hour>=hh)*(scan_hour<hh+1)
    if select_h.sum():
      fig = plt.figure(figsize=(8,4))
      plt.plot(scan_date, scan_corr_flag_frac, 'ko', ms=1, zorder=1, alpha=0.5)
      plt.plot(scan_date[select_h], scan_corr_flag_frac[select_h], 'ro', ms=3, zorder=10)
      plt.legend(['all scans', 'scans started between {0:02d}UT and {1:02d}UT'.format(hh,hh+1)])
      plt.axhline(y=np.median(scan_corr_flag_frac), linestyle=':', color='k', linewidth=1)
      plt.ylim(0,1)
      plt.xlabel('Date')
      plt.ylabel('Per-scan target flagged fraction @ 1390 - 1410 MHz')
      plt.savefig('{0:s}/targetflagfrac_{1:02d}h.png'.format(plot_dir, hh), dpi=200)
      plt.close()

  fig = plt.figure(figsize=(8,4))
  plt.plot(t_start, nr_ant, 'ko', mfc='none', ms=3, alpha=0.5)
  plt.plot(t_start, nr_ant_eff, 'ko', ms=2)
  plt.legend(['in MS', '<100% flagged'], loc='lower center')
  plt.xlabel('Date')
  plt.ylabel('Nr antennas')
  plt.savefig('{0:s}/nrant.png'.format(plot_dir), dpi=200)
  plt.close()

  fig = plt.figure(figsize=(8,4))
  plt.plot(t_start, corr_flag_frac, 'ko', ms=2)
  plt.axhline(y=np.median(corr_flag_frac), linestyle=':', color='k', linewidth=1)
  plt.ylim(0,)
  plt.xlabel('Date')
  plt.ylabel('Target flagged fraction @ 1390-1410 MHz')
  plt.savefig('{0:s}/targetflagfrac.png'.format(plot_dir), dpi=200)
  plt.close()

  fig = plt.figure(figsize=(8,4))
  plt.plot(t_start, cal_corr_flag_frac, 'ko', ms=2)
  plt.axhline(y=np.median(cal_corr_flag_frac), linestyle=':', color='k', linewidth=1)
  plt.ylim(0,)
  plt.xlabel('Date')
  plt.ylabel('Calibrators flagged fraction @ 1350-1430 MHz')
  plt.savefig('{0:s}/calflagfrac.png'.format(plot_dir), dpi=200)
  plt.close()

  fig = plt.figure(figsize=(8,4))
  plt.plot(t_start, line_meas_noise, 'ko', ms=2)
  plt.xlabel('Date')
  plt.ylabel('Line noise @ 30", 5 km/s (mJy/beam)')
  plt.savefig('{0:s}/linenoise.png'.format(plot_dir), dpi=200)
  plt.close()

  fig = plt.figure(figsize=(8,4))
  plt.plot(t_start, cont_noise, 'ko', ms=2)
  plt.xlabel('Date')
  plt.ylabel('Cont noise (muJy/beam)')
  plt.savefig('{0:s}/contnoise.png'.format(plot_dir), dpi=200)
  plt.close()

  fig = plt.figure(figsize=(8,4))
  plt.plot(t_start, np.array(line_meas_noise).astype(float)/np.array(line_exp_noise).astype(float), 'ko', ms=2)
  plt.axhline(y=1, linestyle=':', color='k', linewidth=1)
  plt.xlabel('Date')
  plt.ylabel('Line measured/expected noise ratio @ 30", 5 km/s')
  plt.savefig('{0:s}/measexplinenoise.png'.format(plot_dir), dpi=200)
  plt.close()


def create_json_dict(obs_dir, json_dict):
  if obs_dir[-1]=='/':
    obs_dir = obs_dir[:-1]
  observations = glob('{0:s}/mfs??/mfs_report_track?'.format(obs_dir)) # Hardcoded to work for MFS
  observations.sort()
  zoom_report = {}

  print('# Generating JSON report for {0:d} observations found at {1:s}'.format(len(observations), obs_dir))

  for obs in observations:
    ### Extract MFS field ID and track number from directory name
    mfs_id, track = obs.split(obs_dir)[-1].split('/')[1:3]
    mfs_id = mfs_id.replace('mfs','')
    track = track[-1]
    print('# Working on field mfs{0:s}, track {1:s}'.format(mfs_id, track))
  
    ### Extract MS file name and the nr of antennas from Obsinfo report
    nr_ant, ms_file = None, None
    if not os.path.exists('{0:s}/obsinfoReport.txt'.format(obs)):
      print('    ERROR: obsinfoReport.txt not found. Aborting ...')
      sys.exit()

    for ll in open('{0:s}/obsinfoReport.txt'.format(obs)).readlines():
      if 'antennas in MS:' in ll:
        nr_ant = int(ll.split(' ')[0])
      elif '-elevation-tracks.png' in ll:
        ms_file = ll.split('/')[1].split('_')[0]
    if not nr_ant:
      print('    ERROR: Could not find nr_ant in obsinfoReport.txt. Aborting ...')
      sys.exit()
    elif not ms_file:
      print('    ERROR: Could not find ms_file in obsinfoReport.txt. Aborting ...')
      sys.exit()

    ### Extract start time, end time, scan ID, scan start time and scan end time from MS obsinfo file
    ms_name = '{0:s}/mfs{1:s}/msdir/{2:s}_sdp_l0-obsinfo.txt'.format(obs_dir, mfs_id, ms_file)
    scans_info = {}
    if not os.path.exists(ms_name):
      print('    ERROR: {0:s}_sdp_l0-obsinfo.txt not found. Aborting ...'.format(ms_file))
      sys.exit()

    for ll in open(ms_name).readlines():
      if 'Observed from' in ll:
        time_start = ll.split()[2]
        time_end = ll.split()[4]
      if ' mfs{0:s} '.format(mfs_id) in ll:
        scan = ll.split()[:4]
        scans_info[scan[3]] = [scan[0], scan[2], []]
      if ll[:7] == 'Fields:':
        break

    ### Extract calibrators' flag stats from Xcal report
    nr_flag_ant_cal, raw_flag_frac_cal, corr_flag_frac_cal = 0, 0., 0.
    if not os.path.exists('{0:s}/xcalReport.txt'.format(obs)):
      print('    ERROR: xcalReport.txt not found. Aborting ...')
      sys.exit()

    for ll in open('{0:s}/xcalReport.txt'.format(obs)).readlines():
      if ll[:7] == 'antenna' and '100%' in ll:
        nr_flag_ant_cal += 1
      if ll[:13] == 'Total Flagged':
        raw_flag_frac_cal = float(ll.split('(')[-1].split('%')[0])/100
    exp_flag_frac_cal = 1 - (nr_ant - nr_flag_ant_cal) * (nr_ant - nr_flag_ant_cal - 1) / nr_ant / (nr_ant - 1)
    corr_flag_frac_cal = raw_flag_frac_cal - exp_flag_frac_cal
    calibrators_vis = {
                       'NR_FLAG_ANT': nr_flag_ant_cal,
                       'RAW_FLAG_FRAC': raw_flag_frac_cal,
                       'CORR_FLAG_FRAC': corr_flag_frac_cal,
                       }

    ### Extract target's flag stats per SPW and per scan from Calflag report
    nr_spw = 4  # Hardcoded to work for MFS
    nr_flag_ant_trg = np.zeros((nr_spw,))
    raw_flag_frac_trg = np.zeros((nr_spw,))
    corr_flag_frac_trg = np.zeros((nr_spw,))
    spw_id_trg = []
    if not os.path.exists('{0:s}/calflagReport.txt'.format(obs)):
      print('    ERROR: calflagReport.txt not found. Aborting ...')
      sys.exit()

    for ll in open('{0:s}/calflagReport.txt'.format(obs)).readlines():
      if 'output_spw' in ll:
        spw_index = int(ll.split('output_spw')[-1].split('/')[0].replace('_',''))
        spw_id_trg.append('SPW{0:d}'.format(spw_index))
      elif ll[:7] == 'antenna' and '100%' in ll:
        nr_flag_ant_trg[spw_index] += 1
      elif ll[:4] == 'scan':
        scan_id = ll.split()[1]
        scans_info[scan_id][-1].append(float(ll.split('(')[-1].split('%')[0])/100)
      elif ll[:13] == 'Total Flagged':
        raw_flag_frac_trg[spw_index] = float(ll.split('(')[-1].split('%')[0])/100
    exp_flag_frac_trg = 1 - (nr_ant - nr_flag_ant_trg) * (nr_ant - nr_flag_ant_trg - 1) / nr_ant / (nr_ant - 1)
    corr_flag_frac_trg = raw_flag_frac_trg - exp_flag_frac_trg

    target_vis = {}
    for spw in range(len(spw_id_trg)):
      target_vis[spw_id_trg[spw]] = {
                  'NR_FLAG_ANT': int(nr_flag_ant_trg[spw]),
                  'RAW_FLAG_FRAC': raw_flag_frac_trg[spw],
                  'CORR_FLAG_FRAC': corr_flag_frac_trg[spw],
                                       }
    for scan_id in scans_info:
      target_vis['SCAN_{0:s}'.format(scan_id)] = {
                                                  'SCAN_START': scans_info[scan_id][0].split('/')[-1],
                                                  'SCAN_END': scans_info[scan_id][1].split('/')[-1],
                                                  'RAW_FLAG_FRAC': np.array(scans_info[scan_id][2]).mean(),
                                                  'CORR_FLAG_FRAC': (np.array(scans_info[scan_id][2]) - exp_flag_frac_trg).mean(),
                                                  'RAW_FLAG_FRAC_SPW0': np.array(scans_info[scan_id][2])[0],
                                                  'CORR_FLAG_FRAC_SPW0': (np.array(scans_info[scan_id][2]) - exp_flag_frac_trg)[0],
                                                  'RAW_FLAG_FRAC_SPW1': np.array(scans_info[scan_id][2])[1],
                                                  'CORR_FLAG_FRAC_SPW1': (np.array(scans_info[scan_id][2]) - exp_flag_frac_trg)[1],
                                                  'RAW_FLAG_FRAC_SPW2': np.array(scans_info[scan_id][2])[2],
                                                  'CORR_FLAG_FRAC_SPW2': (np.array(scans_info[scan_id][2]) - exp_flag_frac_trg)[2],
                                                  'RAW_FLAG_FRAC_SPW3': np.array(scans_info[scan_id][2])[3],
                                                  'CORR_FLAG_FRAC_SPW3': (np.array(scans_info[scan_id][2]) - exp_flag_frac_trg)[3],
                                                  }

    ### Extract continuum imaging stats from Cont report
    if not os.path.exists('{0:s}/contReport.txt'.format(obs)):
      print('    ERROR: contReport.txt not found. Aborting ...')
      sys.exit()

    for ll in open('{0:s}/contReport.txt'.format(obs)).readlines():
      if 'clean components' in ll:
        nr_cc = int(ll.split()[0])
      elif 'continuum noise' in ll:
        cont_noise = float(ll.split('=')[-1].split()[0])
    target_cont = {
                   'NR_CLEAN_COMP': nr_cc,
                   'NOISE': cont_noise,
                   }

    ### Extract line imaging stars from Line report
    if not os.path.exists('{0:s}/lineReport.txt'.format(obs)):
      print('    ERROR: lineReport.txt not found. Aborting ...')
      sys.exit()

    for ll in open('{0:s}/lineReport.txt'.format(obs)).readlines():
      if 'Expected noise:' in ll:
        exp_noise_30 = float(ll.split('Expected noise:')[-1].split()[0])
      if 'Measured noise:' in ll:
        meas_noise_30 = float(ll.split('Measured noise:')[-1].split()[0])
      if 'low-res cube:' in ll:
        meas_noise_60 = float(ll.split('low-res cube:')[-1].split()[0])
    target_line = {
                   'EXP_NOISE_30': exp_noise_30,
                   'MEAS_NOISE_30': meas_noise_30,
                   'MEAS_NOISE_60': meas_noise_60,
                   }

    # Add current observation to zoom report dictionary
    zoom_report[ms_file] = {
                            'MFS_ID': mfs_id,
                            'TRACK': track,
                            'TIME_START': time_start,
                            'TIME_END': time_end,
                            'NR_ANT': nr_ant,
                            'CALIBRATORS_VIS': calibrators_vis,
                            'TARGET_VIS': target_vis,
                            'TARGET_CONT': target_cont,
                            'TARGET_LINE': target_line,
                            }

  # Save zoom report dictionary to disc
  if len(observations):
    print('# Writing to disc JSON report for {0:d} observations -> {1:s}'.format(len(observations), json_dict))
    with open(json_dict, 'w') as fp:
      json.dump(zoom_report, fp)



# Run code

args = create_parser().parse_args([a for a in sys.argv[1:]])

if args.create_json:
  create_json_dict(args.obs_dir, args.json_dict)
  
if args.plot_json:
  plot_json_dict(args.json_dict, args.plot_dir)