#!/usr/bin/env python3

import argparse
import datetime
import os
import json

parser = argparse.ArgumentParser()
parser.add_argument('--am', required=True, help='the agp file containing the ALLMAPS output')
parser.add_argument('--marker', required=True, help='the bed file containing the RAD markers')
parser.add_argument('--groups', required=True, help='the agp file cotaining scaffold lengths')
parser.add_argument('--gaps', required=True, help='the bed file containing the gaps')
parser.add_argument('--output', required=True, help='name of the agp file to write output to')
#parser.add_argument('--nr', help='print commented list of marker and agp file of the chromosome with the given number to stdout')
parser.add_argument('--cut', action='store_true', help='cut genetic position of marker after three decimals')
parser.add_argument('--log', action='store_true', help='write a tsv file with all interim results')
args = parser.parse_args()

s_link = 100
g_link = 30000
ss_link = 100000
map_base = ""
timer_start = datetime.datetime.now()

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

mapping = {}
stats = {'unused_pre_ss': [], 'unused_post_ss': [], 'unused_pre_s': [], 'unused_post_s': [], 'used_pre_ss': [], 'used_post_ss': [], 'used_pre_s': [], 'used_post_s': [], 'length_unused_pre_ss': 0, 'length_unused_post_ss': 0, 'length_unused_pre_s': 0, 'length_unused_post_s': 0, 'length_used_pre_ss': 0, 'length_used_post_ss': 0, 'length_used_pre_s': 0, 'length_used_post_s': 0, 'n_unused_pre': 0, 'n_unused_post': 0, 'n_used_pre': 0, 'n_used_post': 0, 'np_unused_pre': 0, 'np_unused_post': 0, 'np_used_pre': 0, 'np_used_post': 0, 's_in_gap': 0}
with open(args.am) as m:
	for line in m:
		content = line.strip().split()
		if "chr" in content[0] and content[4] == "W":
			if content[0] in mapping:
				mapping[content[0]].append({'object_beg': int(content[1]), 'object_end': int(content[2]), 'order_in_chr': int(content[3]), 'component_type': content[4], 'component_id': content[5], 'component_beg': int(content[6]), 'component_end': int(content[7]), 'orientation': content[8]})
			else:
				mapping[content[0]] = [{'object_beg': int(content[1]), 'object_end': int(content[2]), 'order_in_chr': int(content[3]), 'component_type': content[4], 'component_id': content[5], 'component_beg': int(content[6]), 'component_end': int(content[7]), 'orientation': content[8]}]
			if "super" in content[5]:
				stats['used_pre_ss'].append(content[5])
				stats['length_used_pre_ss'] += int(content[7])
			else:
				stats['used_pre_s'].append(content[5])
				stats['length_used_pre_s'] += int(content[7])
		elif not "#" in content[0] and content[4] == "W":
			if "super" in content[0]:
				stats['unused_pre_ss'].append(content[5])
				stats['length_unused_pre_ss'] += int(content[7])
			else:
				stats['unused_pre_s'].append(content[5])
				stats['length_unused_pre_s'] += int(content[7])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def sort_marker(m, obj):
	pre_out = sorted(m, key=lambda k: k['genetic_pos'])
	out = [pre_out[0]]
	buf = []
	for i in range(1, len(pre_out)):
		if pre_out[i]['genetic_pos'] == pre_out[i-1]['genetic_pos'] and pre_out[i]['assembly'] == pre_out[i-1]['assembly']:
			if len(buf) == 0:
				buf = [out.pop(-1), pre_out[i]]
			else:
				buf.append(pre_out[i])
		else:
			if len(buf) > 0:
				component_index = next((index for (index, d) in enumerate(mapping[obj]) if d['component_id'] == pre_out[i-1]['assembly']), None)
				if component_index is not None and mapping[obj][component_index]['orientation'] == '-':
					out.extend(sorted(buf, key=lambda k: k['marker_start_pos'], reverse=True))
				else:
					out.extend(sorted(buf, key=lambda k: k['marker_start_pos']))
				out.append(pre_out[i])
				buf = []
			else:
				out.append(pre_out[i])
	return (out)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def sort_buffer(b, obj):
	if len(b) == 1:
		return (b)
	else:
		out = []
		for entry in b:
			entry['i_in_map'] = next((index for (index, d) in enumerate(mapping[obj]) if d['component_id'] == entry['assembly']), 0)
			out.append(entry)
		out.sort(key = lambda entry: (entry['i_in_map'], entry['index']))
		return (out)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def sort_movable_marker(m, obj):
	out = []
	buf = [m[0]]
	for i in range(1, len(m)):
		if m[i]['genetic_pos'] == buf[-1]['genetic_pos']:
			buf.append(m[i])
		else:
			out.extend(sort_buffer(buf, obj))
			buf = [m[i]]
	out.extend(sort_buffer(buf, obj))
	return (out)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def collapse_marker(ma):
	out = []
	indices = {}
	for i in range(0, len(ma)):
		if ma[i]['assembly'] in indices:
			if ma[i]['assembly'] == ma[i-1]['assembly']:
				indices[ma[i]['assembly']][-1]['to'] = i
			else:
				indices[ma[i]['assembly']].append({'from': i, 'to': i})
		else:
			indices[ma[i]['assembly']] = [{'from': i, 'to': i}]
	i = 0
	i2 = 0
	while i < len(ma):
		entry = indices[ma[i]['assembly']].pop(0)
		next_occurences = indices[ma[i]['assembly']]
		out.append({'id': ma[i]['assembly'], 'index_in_collapsed': i2, 'from': entry['from'], 'to': entry['to'], 'next': next_occurences.copy()})
		i = entry['to'] + 1
		i2 += 1
	return (out)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

marker = {}
with open(args.marker) as ma:
	for line in ma:
		content = line.strip().split()
		content[3] = content[3].split(':')
		if not map_base:
			map_base = ''.join(i for i in content[3][0] if not i.isdigit())
		if args.cut:
			genetic_pos = float(int(float(content[3][1])*1000)/1000)
		else:
			genetic_pos = float(content[3][1])
		if content[3][0] in marker:
			marker[content[3][0]]['marker'].append({'genetic_pos': genetic_pos, 'assembly': content[0], 'marker_start_pos': int(content[2])})
		else:
			marker[content[3][0]] = {'marker': [{'genetic_pos': genetic_pos, 'assembly': content[0], 'marker_start_pos': int(content[2])}]}
for x in marker:
	marker[x]['marker'] = sort_marker(marker[x]['marker'], x.replace(map_base, "chr"))
	for i in range(0, len(marker[x]['marker'])):
		marker[x]['marker'][i]['index'] = i
	marker[x]['marker'] = sort_movable_marker(marker[x]['marker'], x.replace(map_base, "chr"))
	for i in range(0, len(marker[x]['marker'])):
		marker[x]['marker'][i]['index'] = i
	marker[x]['collapsed'] = collapse_marker(marker[x]['marker'])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

grs = {}
with open(args.groups) as gr:
	for line in gr:
		content = line.strip().split()
		if args.groups[-4:] == ".agp":
			grs[content[0]] = {'start': int(content[1]), 'end': int(content[2]), 'type': content[4], 'gaps': []}
		else:
			grs[content[0]] = {'start': 1, 'end': int(content[1]), 'type': 'W', 'gaps': []}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

with open(args.gaps) as g:
	for line in g:
		content = line.strip().split()
		grs[content[0]]['gaps'].append({'gap_id': content[3], 'start': int(content[1]), 'end': int(content[2]), 'length': int(content[4])})

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def main(map_identifier, obj, debug_flag, log_file):
	out = []
	used = [[], [], [], [], []]
	into = {}
	delayed = ""
	last_gap = ""
	marks = check_consent(map_identifier, obj, False)
	if debug_flag and args.log:
		log_file.write(obj + "; cycle 0\n")
		write_log(log_file, marks, debug_flag, False, "check_consent")
		log_file.write("\n")
	for pos_marker_collapsed in range(0, len(marks)):
		if debug_flag and args.log:
			log_file.write(obj + "; cycle " + str(pos_marker_collapsed+1) + "\n")
		out.append(marks[pos_marker_collapsed])
#----- check_earlier_later -----#
		r = check_earlier_later(out, used[0], marks, pos_marker_collapsed, obj)
		out = r[0]
		used[0] = r[1]
		write_log(log_file, out, debug_flag, False, "check_earlier_later")
#----- check_gap ---------------#
		r = check_gap(out, into, map_identifier, used[1])
		out = r[0]
		into = r[1]
		used[1] = r[2]
		write_log(log_file, out, debug_flag, False, "check_gap")
#----- check_insertion ---------#
		r = check_insertion(out, map_identifier, pos_marker_collapsed, used[2], used[4], delayed)
		out = r[0]
		used[2] = r[1]
		delayed = r[2]
		write_log(log_file, out, debug_flag, False, "check_insertion")
#----- resolve_unknown ---------#
		r = resolve_unknown(out, map_identifier, pos_marker_collapsed, used[3])
		out = r[0]
		used[3] = r[1]
		if out[-1]['comment'][0] == "consensus":
			used[4].append(out[-1]['id'])
		write_log(log_file, out, debug_flag, False, "resolve_unknown")
#----- recheck_earlier ---------#
		out = recheck_earlier(out)
		write_log(log_file, out, debug_flag, False, "recheck_earlier")
#----- get_gap -----------------#
		r = get_gap(out, map_identifier, last_gap, used[0])
		out = r[0]
		last_gap = r[1]
		used[0] = r[2]
		write_log(log_file, out, debug_flag, True, "get_gap")
	return (out)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def check_consent(map_identifier, obj, recursion):
	out = []
	used = []
	pos_mapping = 0
	flag = False
	for pos_marker_collapsed in range(0, len(marker[map_identifier]['collapsed'])):
		out.append(marker[map_identifier]['collapsed'][pos_marker_collapsed])
		orientation_index = next((index for (index, d) in enumerate(mapping[obj]) if d['component_id'] == out[-1]['id']), None)
		if orientation_index:
			out[-1]['orientation'] = mapping[obj][orientation_index]['orientation']
		else:
			out[-1]['orientation'] = "?"
		if mapping[obj][pos_mapping]['component_id'] == marker[map_identifier]['collapsed'][pos_marker_collapsed]['id']:
			out[-1]['comment'] = ["consensus"]
			used.append(marker[map_identifier]['collapsed'][pos_marker_collapsed]['id'])
			if pos_mapping+1 < len(mapping[obj]):
				pos_mapping += 1
		elif pos_mapping+1 < len(mapping[obj]) and mapping[obj][pos_mapping+1]['component_id'] == marker[map_identifier]['collapsed'][pos_marker_collapsed]['id']:
			oldindices = [index for (index, d) in enumerate(marker[map_identifier]['collapsed']) if d['id'] == mapping[obj][pos_mapping]['component_id']]
			out.pop(-1)
			if oldindices[0] > pos_marker_collapsed:
				marker[map_identifier]['collapsed'].insert(pos_marker_collapsed, marker[map_identifier]['collapsed'].pop(oldindices[0]))
				out.append(marker[map_identifier]['collapsed'][pos_marker_collapsed])
				used.append(marker[map_identifier]['collapsed'][pos_marker_collapsed]['id'])
				if pos_mapping+1 < len(mapping[obj]):
					pos_mapping += 1
			else:
				flag = True
				marker[map_identifier]['collapsed'].insert(pos_marker_collapsed-1, marker[map_identifier]['collapsed'].pop(oldindices[0]))
			out[-1]['comment'] = ["consensus"]
			out[-1]['orientation'] = mapping[obj][pos_mapping]['orientation']
		else:
			out[-1]['comment'] = ["unknown"]
		out[-1]['used'] = marker[map_identifier]['collapsed'][pos_marker_collapsed]['id'] in used
	if flag:
		if recursion:
			print ("\ntried several recursion steps\n")
			return (out)
		else:
			return (check_consent(map_identifier, obj, True))
	else:
		return (out)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def check_earlier_later(out, used, marks, pmc, obj):
	if out[-1]['id'] in used:
		out[-1]['comment'] = ["ignore (earlier)"]
	if out[-1]['comment'][0] == "consensus":
		used.append(out[-1]['id'])
	out[-1]['in_mapping'] = marks[pmc]['id'] in set([x['component_id'] for x in mapping[obj]])
	if out[-1]['comment'][0] == "unknown" and out[-1]['in_mapping'] and not out[-1]['used']:
		out[-1]['comment'] = ["ignore (later)"]
	return ([out, used])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def check_gap(out, into, map_identifier, used):
	if out[-1]['comment'][0] == "unknown" and not out[-1]['in_mapping'] and not out[-1]['used']:
		if  "super" in out[-1]['id']:
			out[-1]['comment'] = ["ignore (new ss)"]
		elif "super" in into['id']:
			if len(out[-1]['next']) > 0:
				used.append(out[-1]['id'])
			out[-1]['used'] = len(out[-1]['next']) > 0 or out[-1]['id'] in used or out[-1]['used']
			out[-1]['comment'] = ["check gaps in ", into['id'], " from ", marker[map_identifier]['marker'][into['to']]['marker_start_pos']]
	else:
		into = out[-1]
	return ([out, into, used])	

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def check_insertion(out, map_identifier, pmc, used, used_consensus, delayed):
	if out[-1]['id'] == delayed:
		out[-1]['comment'] = ["consensus"]
		delayed = ""
	elif "super" in out[-1]['id']:
		if out[-1]['comment'][0] == "ignore (earlier)":
			sub_pos = pmc
			insertions = []
			flag = False
			flag2 = False
			new_ss = []
			while sub_pos >= 0 and not (out[sub_pos]['comment'][0] == "consensus" or "insertion site" in out[sub_pos]['comment'][0]):
				if out[sub_pos]['comment'][0] == "unknown":
					flag = True
				elif "check gaps" in out[sub_pos]['comment'][0]:
					insertions.append(sub_pos)
					flag = True
				elif out[sub_pos]['comment'][0] == "ignore (new ss)":
					new_ss.append(sub_pos)
				elif out[sub_pos]['comment'][0] == "ignore (earlier)" and out[sub_pos]['id'] == out[-1]['id']:
					if flag2:
						out[sub_pos]['comment'] = ["consensus2"]
					flag2 = True
				sub_pos -= 1
			if out[sub_pos]['id'] == out[-1]['id']:
				if flag:
					for x in insertions:
						if out[-1]['id'] not in out[x]['comment'][1]:
							if len(out[x]['next']) > 0:
								used.append(out[x]['id'])
							out[x]['used'] = len(out[x]['next']) > 0 or out[x]['id'] in used or out[x]['used']
							out[x]['comment'] = ["check gaps2 in ", out[-1]['id'], " from ", marker[map_identifier]['marker'][out[sub_pos]['to']]['marker_start_pos']]
							out[x]['orientation'] = "?"
					out[-1]['comment'] = ["insertion site"]
				else:
					back_counter = len(out) - 2
					while not out[back_counter]['id'] == out[-1]['id']:
						if out[back_counter]['comment'][0] == "ignore (new ss)":
							out[back_counter]['comment'] = ["ignore (nested ss)"]
						back_counter -= 1
					out[-1]['comment'] = ["consensus"]
				for ss in new_ss:
					out[ss]['comment'] = ["ignore (nested ss)"]
			elif sub_pos-1 >= 0 and used_consensus.count(out[sub_pos]['id']) == 1 and out[sub_pos-1]['id'] == out[-1]['id'] and out[sub_pos]['to']-out[sub_pos]['from']+1 == 1 and len(out[sub_pos]['next']) > 0 and out[-1]['to']+1 == out[sub_pos]['next'][0]['from']:
				out[sub_pos]['comment'] = ["ignore (delay)"]
				out[-1]['comment'] = ["consensus"]
				used_consensus.pop(-1)
				delayed = out[sub_pos]['id']
			else:
				while sub_pos < pmc:
					if out[sub_pos]['comment'][0] == "consensus2":
						out[sub_pos]['comment'] = ["ignore (earlier)"]
					sub_pos += 1
		elif out[-1]['comment'][0] == "consensus":
			sub_pos = pmc - 1
			while sub_pos >= 0 and "super" not in out[sub_pos]['id']:
				if "check gaps" in out[sub_pos]['comment'][0] and out[-1]['id'] not in out[sub_pos]['comment'][0]:
					out[sub_pos]['comment'] = ["new1"]
				sub_pos -= 1
	return ([out, used, delayed])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def resolve_unknown(out, map_identifier, pmc, used):
	if out[-1]['comment'][0] == "consensus" or "insertion site" in out[-1]['comment'][0]:
		sub_pos = pmc - 1
		unknowns = []
		while sub_pos >= 0 and not (out[sub_pos]['comment'][0] == "consensus" or "insertion site" in out[sub_pos]['comment'][0]):
			if out[sub_pos]['comment'][0] == "unknown":
				unknowns.append(sub_pos)
			sub_pos -= 1
		for i in unknowns:
			message = ""
			if out[sub_pos]['id'] == out[-1]['id']:
				message = ["check gaps3 in ", out[-1]['id'], " from ", marker[map_identifier]['marker'][out[sub_pos]['to']]['marker_start_pos']]
			else:
				message = ["new2"]
				out[i]['used'] = True
			if len(out[i]['next']) > 0:
				used.append(out[i]['id'])
			out[i]['used'] = len(out[i]['next']) > 0 or out[i]['id'] in used or out[i]['used']
			out[i]['comment'] = message
	return ([out, used])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def recheck_earlier(out):
	if out[-1]['comment'][0] == "ignore (earlier)":
		for pos in range(len(out)-2, -1, -1):
			if not "ignore" in out[pos]['comment'][0]:
				if out[pos]['id'] == out[-1]['id']:
					out[-1]['comment'] = ["consensus"]
				break
	return (out)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_gap(out, map_identifier, last_gap, used_new):
	if out[-1]['comment'][0] == "insertion site":
		sub_pos = len(out) - 2
		insertions = []
		insertion_length = 0
		while not out[sub_pos]['id'] == out[-1]['id']:
			if "check gaps" in out[sub_pos]['comment'][0]:
				insertions.append(sub_pos)
				insertion_length += grs[out[sub_pos]['id']]['end']
			sub_pos -= 1
		if len(insertions) >= 1:
			gap = get_fitting_gap(out[-1]['id'], insertion_length, out[insertions[0]]['comment'][-1], marker[map_identifier]['marker'][out[-1]['from']]['marker_start_pos'])
		current_gap = ""
		for x in insertions:
			if out[x]['id'] not in used_new:
				if gap:
					out[x]['comment'] = ["insert in ", out[x]['comment'][1], " in ", gap['gap_id'], ": ", gap['start'], "-", gap['end']]
					stats['s_in_gap'] += 1 #grs[out[x]['id']]['end']
					if out[x]['used'] and not out[x]['in_mapping']:
						used_new.append(out[x]['id'])
					if gap['gap_id'] == last_gap and not current_gap == last_gap:
						out[sub_pos]['comment'] = ["ignore (gap conflict)"]
					current_gap = gap['gap_id']
					last_gap = gap['gap_id']
				else:
					into = out[x]['comment'][1]
					out[x]['comment'] = ["no fitting gaps"]
					index = out[insertions[-1]]['to']
					while can_be_moved(map_identifier, index):
						if not out[-1]['id'] == marker[map_identifier]['marker'][index+1]['assembly']:
							comment = out[next((i for (i, d) in enumerate(out) if d['from'] == index+2), None)]['comment'][0]
							if "consensus" in comment:
								out[x]['comment'] = ["moved out, new"]
								if out[x]['used'] and not out[x]['in_mapping']:
									used_new.append(out[x]['id'])
						else:
							gap = get_fitting_gap(out[-1]['id'], insertion_length, marker[map_identifier]['marker'][index+1]['marker_start_pos'], marker[map_identifier]['marker'][index+2]['marker_start_pos'])
							if gap:
								out[x]['comment'] = ["moved to index ", index+1, ", insert in ", out[-1]['id'], " in ", gap['gap_id'], ": ", gap['start'], "-", gap['end']]
								stats['s_in_gap'] += 1 #grs[out[x]['id']]['end']
								if out[x]['used'] and not out[x]['in_mapping']:
									used_new.append(out[x]['id'])
								if gap['gap_id'] == last_gap and not current_gap == last_gap:
									out[sub_pos]['comment'] = ["ignore (gap conflict)"]
								current_gap = gap['gap_id']
								last_gap = gap['gap_id']
								break
						index += 1
			else:
				out[x]['comment'] = "ignore (earlier)"
	return ([out, last_gap, used_new])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def can_be_moved(map_identifier, index):
	return (index+1 < len(marker[map_identifier]['marker']) and marker[map_identifier]['marker'][index]['genetic_pos'] == marker[map_identifier]['marker'][index+1]['genetic_pos'])

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def get_fitting_gap(gr, l, s, e): # group, length_of_required_gap, start_of_area_in_group_to_be_checked, end_of_area_in_group_to_be_checked
	if s > e:
		b = s
		s = e
		e = b
	gaps = []
	for gap in grs[gr]['gaps']:
		gap['gr'] = gr
		if gap['start'] >= s and gap['length'] >= l:
			gaps.append(gap)
		if gap['start'] > e:
			break
	if gaps:
		return (max(gaps, key=lambda x: x['length']))
	else:
		return (None)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def write_log(log_file, out, flag, f2, step):
	if flag and args.log:
		log_file.write(step + "\nindex_in_collapsed\tid\torientation\tfrom\tto\tin_mapping\tused\tnext\tcomment\n")
		for line in out:
			log_file.write(str(line['index_in_collapsed']) + "\t" + str(line['id']) + "\t" + str(line['orientation']) + "\t" + str(line['from']) + "\t" + str(line['to']) + "\t")
			if 'in_mapping' in line:
				log_file.write(str(line['in_mapping']) + "\t")
			else:
				log_file.write("na\t")
			log_file.write(str(line['used']) + "\t" + str(len(line['next'])) + "\t")
			for x in line['comment']:
				log_file.write(str(x) + " ")
			log_file.write("\n")
		if f2:
			log_file.write("\n")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def create_agp(commented_list, obj):
	collapsed_list = collapse_output(commented_list)
	out = []
	last = ""
	for i in range(0, len(collapsed_list)):
		out.append({'object': obj, 'component_type': grs[collapsed_list[i]['id']]['type'], 'component_id/gap_length': collapsed_list[i]['id'], 'component_beg/gap_type': 1, 'component_end/linkage': grs[collapsed_list[i]['id']]['end'], 'orientation/linkage_evidence': collapsed_list[i]['orientation'], 'indices': collapsed_list[i]['indices']})
		if "super" in collapsed_list[i]['id'] and len(collapsed_list[i]['indices']) > 1:
			if i < collapsed_list[i]['indices'][-1]:
				if collapsed_list[i]['orientation'] == "-":
					out[-1]['component_beg/gap_type'] = collapsed_list[i+1]['comment'][-1] + 1
				else:
					out[-1]['component_end/linkage'] = collapsed_list[i+1]['comment'][-3]
			if i > collapsed_list[i]['indices'][0]:
				if collapsed_list[i]['orientation'] == "-":
					out[-1]['component_end/linkage'] = collapsed_list[i-1]['comment'][-3]
				else:
					out[-1]['component_beg/gap_type'] = collapsed_list[i-1]['comment'][-1] + 1
				sub_pos = len(out) - 2
				insert_length = 0
				while not str(out[sub_pos]['component_id/gap_length']) == str(out[-1]['component_id/gap_length']):
					if "scaffold" in str(out[sub_pos]['component_id/gap_length']):
						insert_length += out[sub_pos]['component_end/linkage']
					else:
						insert_length += out[sub_pos]['component_id/gap_length']
					sub_pos -= 1
				stats['n_used_post'] -= insert_length
				if collapsed_list[i]['orientation'] == "-":
					gap_rest = out[sub_pos]['component_beg/gap_type'] - out[-1]['component_end/linkage'] - 1 - insert_length + 100
				else:
					gap_rest = out[-1]['component_beg/gap_type'] - out[sub_pos]['component_end/linkage'] - 1 - insert_length + 100
				if gap_rest >= g_link:
					out[-2]['component_id/gap_length'] = gap_rest
				else:
					out[-2]['component_id/gap_length'] = g_link
					stats['n_used_post'] += g_link - gap_rest
					stats['length_used_post_ss'] += g_link - gap_rest
		if collapsed_list[i]['comment'][0] == "insert in " or collapsed_list[i]['comment'][0] == "moved to index ":
			out[-2]['component_id/gap_length'] = g_link
		elif "super" in last and "super" in collapsed_list[i]['id']:
			out[-2]['component_id/gap_length'] = ss_link
		out.append({'object': obj, 'component_type': "N", 'component_id/gap_length': s_link, 'component_beg/gap_type': "scaffold", 'component_end/linkage': "yes", 'orientation/linkage_evidence': "map"})
		last = collapsed_list[i]['id']
	out.pop(-1)
	out[0]['object_beg'] = 1
	out[0]['object_end'] = out[0]['component_end/linkage'] - out[0]['component_beg/gap_type'] + 1
	out[0]['order_in_chr'] = 1
	for i in range(2, len(out)+1):
		out[i-1]['order_in_chr'] = i
		out[i-1]['object_beg'] = out[i-2]['object_end'] + 1
		if out[i-1]['component_type'] == "W":
			out[i-1]['object_end'] = out[i-2]['object_end'] + out[i-1]['component_end/linkage'] - out[i-1]['component_beg/gap_type'] + 1
		else:
			out[i-1]['object_end'] = out[i-2]['object_end'] + out[i-1]['component_id/gap_length']
	return (out)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def collapse_output(l):
	out = []
	last = {'id': ""}
	for entry in l:
		if ("consensus" in entry['comment'][0] or "insert" in entry['comment'][0] or entry['comment'][0] == "new1" or entry['comment'][0] == "new2" or entry['comment'][0] == "moved to index ") and not last['id'] == entry['id']:
			out.append(entry)
			last = entry
	for entry in out:
		entry['indices'] = [index for (index, d) in enumerate(out) if d['id'] == entry['id']]
	return (out)

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def write_agp(agps):
	used_groups = []
	with open(args.output, 'w') as f:
		f.write("# COMMENT: Generated by AMrefiner (" + str(timer_start.strftime("%Y-%m-%d_%H-%M")) + ")\n")
		f.write("# COMMAND: parser.py --map " + args.am + " --marker " + args.marker + " --groups " + args.groups + " --gaps " + args.gaps + " --output " + args.output)
		if args.cut:
			f.write(" --cut")
		if args.log:
			f.write(" --log")
		f.write("\n# FIELDS: object, object_beg, object_end, part_number, component_type, component_id/gap_length, component_beg/gap_type, component_end/linkage, orientation/linkage_evidence\n")
		for c in agps:
			for line in c:
					f.write(str(line['object']) + "\t" + str(line['object_beg']) + "\t" + str(line['object_end']) + "\t" + str(line['order_in_chr']) + "\t" + str(line['component_type']) + "\t" + str(line['component_id/gap_length']) + "\t" + str(line['component_beg/gap_type']) + "\t" + str(line['component_end/linkage']) + "\t" + str(line['orientation/linkage_evidence']) + "\n")
					if isinstance(line['component_id/gap_length'], str) and line['component_id/gap_length'] not in used_groups:
						used_groups.append(line['component_id/gap_length'])
		for k in sorted(grs.keys()):
			if k in used_groups:
				if "super" in k:
					stats['used_post_ss'].append(k)
					stats['length_used_post_ss'] += grs[k]['end']
				else:
					stats['used_post_s'].append(k)
					stats['length_used_post_s'] += grs[k]['end']
			else:
				f.write(k + "\t" + str(grs[k]['start']) + "\t" + str(grs[k]['end']) + "\t1\t" + grs[k]['type'] + "\t" + k + "\t" + str(grs[k]['start']) + "\t" + str(grs[k]['end']) + "\t?\n")
				if "super" in k:
					stats['unused_post_ss'].append(k)
					stats['length_unused_post_ss'] += grs[k]['end']
				else:
					stats['unused_post_s'].append(k)
					stats['length_unused_post_s'] += grs[k]['end']

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def print_agp(agps):
	print ("object\tobject_beg\tobject_end\torder_in_chr\tcomponent_type\tcomponent_id/gap_length\tcomponent_beg/gap_type\tcomponent_end/linkage\torientation/linkage_evidence")
	for entry in agps:
		print (str(entry['object']) + "\t" + str(entry['object_beg']), end="")
		if entry['object_beg'] < 10000000:
			print ("\t\t", end ="")
		else:
			print ("\t", end ="")
		print (str(entry['object_end']), end="")
		if entry['object_end'] < 10000000:
			print ("\t\t\t", end ="")
		else:
			print ("\t\t", end ="")
		print (str(entry['order_in_chr']) + "\t\t" + str(entry['component_type']) + "\t" + str(entry['component_id/gap_length']), end="\t\t")
		if "scaffold" not in str(entry['component_id/gap_length']):
			print ("\t", end="")
		if "super" not in str(entry['component_id/gap_length']):
			print ("\t", end="")
		print(str(entry['component_beg/gap_type']), end="\t")
		if not entry['component_beg/gap_type'] == "scaffold" and entry['component_beg/gap_type'] < 10000000:
			print ("\t\t", end="")
		else:
			print ("\t", end="")
		print(str(entry['component_end/linkage']), end="")
		if entry['component_end/linkage'] == "yes" or entry['component_end/linkage'] < 10000000:
			print ("\t\t\t", end="")
		else:
			print ("\t\t", end="")
		print (str(entry['orientation/linkage_evidence']))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def process_stats():
	with open(args.gaps) as g:
		for line in g:
			content = line.strip().split()
			if content[0] in stats['used_pre_ss']:
				stats['n_used_pre'] += int(content[4])
			elif content[0] in stats['used_pre_s']:
				stats['n_used_pre'] += int(content[4])
			elif content[0] in stats['unused_pre_ss']:
				stats['n_unused_pre'] += int(content[4])
			elif content[0] in stats['unused_pre_s']:
				stats['n_unused_pre'] += int(content[4])
			if content[0] in stats['used_post_ss']:
				stats['n_used_post'] += int(content[4])
			elif content[0] in stats['used_post_s']:
				stats['n_used_post'] += int(content[4])
			elif content[0] in stats['unused_post_ss']:
				stats['n_unused_post'] += int(content[4])
			elif content[0] in stats['unused_post_s']:
				stats['n_unused_post'] += int(content[4])
	if (stats['length_used_pre_ss'] + stats['length_used_pre_s']) == 0:
		stats['np_used_pre'] = 0.0000
	else:
		stats['np_used_pre'] = float("{0:.4f}".format((stats['n_used_pre']/(stats['length_used_pre_ss'] + stats['length_used_pre_s']))*100))
	if (stats['length_unused_pre_ss'] + stats['length_unused_pre_s']) == 0:
		stats['np_unused_pre'] = 0.0000
	else:
		stats['np_unused_pre'] = float("{0:.4f}".format((stats['n_unused_pre']/(stats['length_unused_pre_ss'] + stats['length_unused_pre_s']))*100))
	if (stats['length_used_post_ss'] + stats['length_used_post_s']) == 0:
		stats['np_used_post'] = 0.0000
	else:
		stats['np_used_post'] = float("{0:.4f}".format((stats['n_used_post']/(stats['length_used_post_ss'] + stats['length_used_post_s']))*100))
	if (stats['length_unused_post_ss'] + stats['length_unused_post_s']) == 0:
		stats['np_unused_post'] = 0.0000
	else:
		stats['np_unused_post'] = float("{0:.4f}".format((stats['n_unused_post']/(stats['length_unused_post_ss'] + stats['length_unused_post_s']))*100))

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def write_stats():
	with open(os.path.dirname(args.output) + "AMrefiner_" + str(timer_start.strftime("%Y-%m-%d_%H-%M")) + "_stats.txt", 'w') as f:
		f.write("\t\t\t_____________________________________________________________________________\n")
		f.write("\t\t\t| original Allmaps\t\t      | AMrefiner\t\t\t    |\n")
		f.write("\t\t\t| used\t\t   | unused\t      | used\t\t | unused\t    |\n")
		f.write("________________________|__________________|__________________|__________________|__________________|\n")
		f.write("                        |                  |                  |                  |                  |\n")
		f.write("#Superscaffolds\t\t| " + unify_length(len(stats['used_pre_ss'])) + "| " + unify_length(len(stats['unused_pre_ss'])) + "| " + unify_length(len(stats['used_post_ss'])) + "| " + unify_length(len(stats['unused_post_ss'])) + "|\n")
		f.write("#Scaffolds\t\t| " + unify_length(len(stats['used_pre_s'])) + "| " + unify_length(len(stats['unused_pre_s'])) + "| " + unify_length(len(stats['used_post_s'])) + "| " + unify_length(len(stats['unused_post_s'])) + "|\n")
		f.write("length Superscaffolds\t| " + unify_length(stats['length_used_pre_ss']) + "| " + unify_length(stats['length_unused_pre_ss']) + "| " + unify_length(stats['length_used_post_ss']) + "| " + unify_length(stats['length_unused_post_ss']) + "|\n")
		f.write("length Scaffolds\t| " + unify_length(stats['length_used_pre_s']) + "| " + unify_length(stats['length_unused_pre_s']) + "| " + unify_length(stats['length_used_post_s']) + "| " + unify_length(stats['length_unused_post_s']) + "|\n")
		f.write("#N\t\t\t| " + unify_length(stats['n_used_pre']) + "| " + unify_length(stats['n_unused_pre']) + "| " + unify_length(stats['n_used_post']) + "| " + unify_length(stats['n_unused_post']) + "|\n")
		f.write("%N\t\t\t| " + unify_length(stats['np_used_pre']) + "| " + unify_length(stats['np_unused_pre']) + "| " + unify_length(stats['np_used_post']) + "| " + unify_length(stats['np_unused_post']) + "|\n")
		f.write("#Scaffolds in gaps\t| " + unify_length(0) + "| " + unify_length(0) + "| " + unify_length(stats['s_in_gap']) + "| " + unify_length(0) + "|\n")
		f.write("________________________|__________________|__________________|__________________|__________________|\n")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def unify_length(num):
	out = ""
	flag = isinstance(num, int)
	cnt = 0
	for c in str(num)[::-1]:
		out += c
		if flag:
			cnt += 1
			if cnt == 3:
				cnt = 0
				out += ","
		else:
			flag = (c == ".")
	if out[-1] == ",":
		out = out[:-1]
	while len(out) < 16:
		out += " "
	return (out[::-1] + " ")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def print_results(results):
	print ("index\tid\t\t\torientation\tfrom\tto\tin map\tused\tnext\t\tcomment")
	for result in results:
		print (str(result['index_in_collapsed']) + "\t" + result['id'], end="\t")
		if not "super" in result['id']:
			print ("\t", end="")
		print ("\t" + str(result['orientation']) + "\t" + str(result['from']) + "\t" + str(result['to']), end="\t")
		print (str(result['in_mapping']) + "\t" + str(result['used']) + "\t" + str(len(result['next'])), end="\t")
		if len(result['next']) > 0:
			print (str(result['next'][0]['from']), end ="\t")
		else:
			print ("", end="\t")
		for c in result['comment']:
			print (str(c), end="")
		print ("")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

agps = []
log_file = None
#if args.nr:
#	ob = "chr" + str(args.nr)
#	map_identifier = map_base + str(args.nr)
#	r = main(map_identifier, ob, False, log_file)
#	print_results(r)
#	print ("\n")
#	print_agp(create_agp(r, ob))
timer_loaded = datetime.datetime.now()
print ("\nloaded files in:\t" + str((timer_loaded - timer_start).total_seconds()) + "s")
if args.log:
	log_file = open(os.path.dirname(args.output) + "AMrefiner_" + str(timer_start.strftime("%Y-%m-%d_%H-%M")) + ".log", 'w')
for i in range(1, len(mapping)+1):
	map_identifier = map_base + str(i)
	ob = "chr" + str(i)
	agps.append(create_agp(main(map_identifier, ob, True, log_file), ob))
if log_file:
	log_file.close()
timer_processed = datetime.datetime.now()
print ("processed data in:\t" + str((timer_processed - timer_loaded).total_seconds()) + "s")
write_agp(agps)
timer_wrote = datetime.datetime.now()
print ("wrote output file in:\t" + str((timer_wrote - timer_processed).total_seconds()) + "s")
process_stats()
timer_processed_stats = datetime.datetime.now()
print ("processed stats in:\t" + str((timer_processed_stats - timer_wrote).total_seconds()) + "s")
write_stats()
timer_wrote_stats = datetime.datetime.now()
print ("wrote stats file in:\t" + str((timer_wrote_stats - timer_processed_stats).total_seconds()) + "s")
