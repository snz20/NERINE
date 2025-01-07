'''
Author: Sumaiya Nazeen (sumaiya_nazeen@hms.harvard.edu)
This file contains auxiliary methods for NERINE.
'''

import sys
import os
import numpy as np
import pandas as pd
from collections import Counter
import gzip
import re
import pickle
import bz2
import json
from scipy.stats import spearmanr
import sys
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns



#define globals
damaging = ['frameshift_variant','start_lost','stop_gained','stop_lost','splice_donor_variant','splice_acceptor_variant','splice_region_variant','structural_interaction_variant','initiator_codon_variant', 'disruptive_inframe_deletion', 'disruptive_inframe_insertion']
excl = ['3_prime_UTR_variant', '5_prime_UTR_premature_start_codon_gain_variant', '5_prime_UTR_variant', 'bidirectional_gene_fusion', 'downstream_gene_variant', 'gene_fusion',  'intergenic_region', 'intron_variant', 'non_coding_transcript_exon_variant', 'sequence_feature', 'TF_binding_site_variant', 'upstream_gene_variant']
incl = ['missense_variant','synonymous_variant'] + damaging
csq_format = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|SOURCE|SIFT|PolyPhen|LoF|LoF_filter|LoF_flags|LoF_info|CADD_PHRED|CADD_RAW|Ancestral_allele|CADD_phred|GTEx_V8_gene|GTEx_V8_tissue|Interpro_domain|Polyphen2_HVAR_pred|SIFT_pred|aaalt|aapos|aaref|clinvar_hgvs|gnomAD_exomes_AFR_AF|gnomAD_exomes_AMR_AF|gnomAD_exomes_EAS_AF|gnomAD_exomes_FIN_AF|gnomAD_exomes_NFE_AF|gnomAD_exomes_POPMAX_AF|gnomAD_exomes_SAS_AF|gnomAD_genomes_AFR_AF|gnomAD_genomes_AMR_AF|gnomAD_genomes_EAS_AF|gnomAD_genomes_FIN_AF|gnomAD_genomes_NFE_AF|gnomAD_genomes_POPMAX_AF|gnomAD_genomes_SAS_AF|hg19_chr|hg19_pos(1-based)|gnomADgc|gnomADgc_AF_nfe|gnomADgc_AF_fin|gnomADgc_AF_amr|gnomADgc_AF_eas|gnomADgc_AF_sas|gnomADgc_AF_afr|gnomADgc_AF|gnomADec|gnomADec_AF_nfe|gnomADec_AF_fin|gnomADec_AF_amr|gnomADec_AF_eas|gnomADec_AF_sas|gnomADec_AF_afr|gnomADec_AF"
csq_fields = csq_format.split("|")		
dbgnom_g = csq_fields.index("gnomAD_genomes_POPMAX_AF")
dbgnom_e = csq_fields.index("gnomAD_exomes_POPMAX_AF")
dbgnom_ng = csq_fields.index("gnomAD_genomes_NFE_AF")
dbgnom_ne = csq_fields.index("gnomAD_exomes_NFE_AF")
gnom_g = csq_fields.index("gnomADgc_AF")
gnom_ng = csq_fields.index("gnomADgc_AF_nfe")
gnom_e = csq_fields.index("gnomADec_AF")
gnom_ne = csq_fields.index("gnomADec_AF_nfe")
dbpphen_ind = csq_fields.index("Polyphen2_HVAR_pred")
dbsift_ind = csq_fields.index("SIFT_pred")
pphen_ind = csq_fields.index("PolyPhen") 
sift_ind = csq_fields.index("SIFT")
dbaaref = csq_fields.index("aaref")
dbaaalt = csq_fields.index("aaalt")
dbaapos = csq_fields.index("aapos")
ppchg_ind = csq_fields.index("Amino_acids")
pppos_ind = csq_fields.index("Protein_position")
gene_ind = csq_fields.index("SYMBOL")
impact_ind = csq_fields.index("IMPACT")
csq_ind = csq_fields.index("Consequence")
dbcadd_ind = csq_fields.index("CADD_phred")
cadd_ind = csq_fields.index("CADD_PHRED") 

impact_map = {'HIGH':3,'MODERATE':2,'LOW':1,'MODIFIER':0}
mycmap = {4: "#E64B35E5", 3: "#E64B35BF", -4: "#3C5488E5", -3: "#3C5488BF", 2: "#E64B3599", -2: "#3C548899", 1: "#E64B3566", -1: "#3C548866", 0: "#8491B433"}

			
def convert_gzvcf_to_mutations(infile, outprefix):
	'''
	This program prepares the input mutations file from a annotated gzipped vcf file.
	Input:
	- gzipped vcf file: The vcf file must be annotated with VEP and dbNSFP v4.3a database.
	- prefix of the output file
	Output: 
	- outprefix_mutations.tsv
	'''
	of = open(outprefix+"_mutations.tsv", "w")
	genes = []
	mut_rec = []

	of.write("#Gene\tmutid\tchg\tnchg\tatod\tppchg\trsid\ttop_csq\tcsq\timpact\tpolyphen\tsift\tcadd\tcsq2\tgg_af\tgg_nfe_af\tge_af\tge_nfe_af\tac\tan\tPopMAF\tac_case\tan_case\tac_control\tan_control\tmutated_individuals\n")
	with gzip.open(infile,"rt") as fp:
		l = fp.readline().strip()
		counter = 0
		vep_field_names = None
		header = None
		while(l):
			counter += 1
			if l[0:2] == "##":
				if l.find('ID=CSQ') > -1:
					vep_field_names = l.split('Format: ')[-1].strip('">').split('|')
				l = fp.readline().strip()		
				continue
			elif l[0] == '#':
				info_ind = -1
				indiv_start = -1
				header = l.split('\t')
				header = dict(zip(header, range(len(header))))
				indiv_start = header['FORMAT'] + 1
				indivs = list(header.keys())[indiv_start:]
				N = len(indivs)
				l = fp.readline().strip()
				continue
			elif l[0] != '#':
				#print("here")
				fields = l.split('\t')
				#print(x)
				chr = fields[0]
				pos = fields[1]
				if chr+':'+pos not in mut_rec:
					mut_rec.append(chr+':'+pos)
				else:
					#print("here")
					l = fp.readline().strip()
					continue
				rsid = fields[2]
				ref = fields[3]
				alt = fields[4]
				gts = fields[indiv_start:]
				hets = 0
				homs = 0
				missing = 0
				mutated_indiv = []
				ac_case = -1
				ac_control = -1
				m_case = -1
				m_control = -1
				#print(len(gts))
				for i in range(len(gts)):
					#print(gts)
					if gts[i] == '0/1' or gts[i] == '1/0':
						hets += 1
						mutated_indiv.append(indivs[i])
						#if indivs[i] in cases:
						#	ac_case += 1
						#elif indivs[i] in controls:
						#	ac_control += 1
					elif gts[i] == '1/1':
						homs += 1
						mutated_indiv.append(indivs[i])
						mutated_indiv.append(indivs[i])
						#if indivs[i] in cases:
                                                #        ac_case += 2
						#elif indivs[i] in controls:
                                                #        ac_control += 2
					elif gts[i] == './.':
						missing += 1
						#if indivs[i] in cases:
                                                #        m_case += 1
						#elif indivs[i] in controls:
                                                #        m_control += 1
				#print("here")
				af = 0.0
				ac = 2*homs + hets
				an = 2*(N-missing)
				af = ac*1.0/an
				an_case = -1
				an_control = -1
				#an_case = 2*(len(cases)-m_case)
				#an_control = 2*(len(controls)-m_control)
				#if af > 0.5:
				#	af = 1-af
				mutid = chr+':'+pos
				mut = ref + '>' + alt
				if len(ref) == 1 and len(alt) == 1:
					chg = 'SNP'
				else:
					chg = 'INDEL'
				print(chr, pos, ref, alt, af, ac_case, ac_control, len(mutated_indiv))
				info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in 
re.split(r';(?=\w)', fields[header['INFO']])])
				#print(info_field)
				if 'CSQ' not in info_field.keys(): 
					l = fp.readline().strip()
					continue
				annotations = [dict(zip(vep_field_names, x.split('|'))) for x in info_field['CSQ'].split(',') if len(vep_field_names) == len(x.split('|'))]
				#lof_annotations = [x for x in annotations if x['LoF'] == 'HC']
				gnomgAF = 0.0
				gnomgAF_nfe = 0.0
				gnomeAF = 0.0
				gnomeAF_nfe = 0.0
				pphen = "-"
				sift = "-"
				cadd = 0.0
				ann_mane_ind = -1
				ann_imp = -1
				for ai in range(len(annotations)):
					a = annotations[ai]
					#if a['Feature'] == tr_mane:
					#	ann_mane_ind = ai
					#	break
					if a['CANONICAL'] == 'YES' or a['MANE_SELECT'] != '':
						if impact_map[a['IMPACT']] > ann_imp:
							ann_mane_ind = ai
							ann_imp = impact_map[a['IMPACT']]
						elif impact_map[a['IMPACT']] == ann_imp and ann_imp == 0:
							cs = a['Consequence'].split('&')[0]
							if cs in incl:
								ann_mane_ind = ai
				
				if ann_mane_ind == -1:
					print("mane select transcript doesn't contain this variant. Taking top annotation.")
					ann_mane_ind = 0
					
				ann = annotations[ann_mane_ind]
				#print(ann.keys())
				gene = ann['SYMBOL']
				impact = ann['IMPACT']
				csq = ann['Consequence']
				lof_annotations = ann['LoF']
				top_csq = csq.split('&')[0]
				if ann['Polyphen2_HVAR_pred'] != "":
					pphen = ann['Polyphen2_HVAR_pred']
				elif ann['PolyPhen'] != "":
					pphen = ann['PolyPhen']
				if ann['SIFT_pred'] != "":
					sift = ann['SIFT_pred']
				elif ann['SIFT'] != "":						
					sift = ann['SIFT']
				if ann['CADD_phred'] != "":
					cadd = ann['CADD_phred']
				elif ann['CADD_PHRED'] != "":
					cadd = ann['CADD_PHRED']
				print(top_csq)
				#if top_csq not in incl:
				#	l = fp.readline().strip()
				#	continue
				csq2 = 'unknown'
				if top_csq in damaging or 'D' in pphen or 'P' in pphen or 'deleterious' in pphen or 'D' in sift or 'deleterious' in sift:
					csq2 = 'damaging'
				elif ('B' in pphen or 'benign' in pphen) and ('T' in sift or 'tolerated' in sift):
					csq2 = 'neutral'
				elif top_csq == 'synonymous_variant':
					csq2 = 'synonymous'
				anc = ann['Ancestral_allele']
				atod = '-'
				if ref == anc:
					atod = anc + '>' + alt
				elif alt == anc:
					atod = anc + '>' + ref
				ppchg = '-'
				if ann['aaref'] != '':
					ppchg = ann['aaref']+ann['aapos']+ann['aaalt']
				elif ann['Amino_acids'] != '':
					ppchg = ann['Amino_acids'] + "_" + ann["Protein_position"]
				if ann["gnomAD_genomes_POPMAX_AF"] != "":
					dbg = ann["gnomAD_genomes_POPMAX_AF"].split(',')[0]
				else:
					dbg = "-"
				if ann["gnomAD_exomes_POPMAX_AF"] != "":
					dbe = ann["gnomAD_exomes_POPMAX_AF"].split(',')[0]
				else:
					dbe = "-"
				if ann["gnomAD_genomes_NFE_AF"] != "":
					dbng = ann["gnomAD_genomes_NFE_AF"].split(',')[0]
				else:
					dbng = "-"
				if ann["gnomAD_exomes_NFE_AF"] != "":
					dbne = ann["gnomAD_exomes_NFE_AF"].split(',')[0]
				else:
					dbne = "-"
				if ann["gnomADgc_AF"] != "":
					gg = ann["gnomADgc_AF"].split(',')[0].split('&')[0]
				else:
					gg = "-"
				if ann["gnomADec_AF"] != "":
					ge = ann["gnomADec_AF"].split(',')[0].split('&')[0]
				else:
					ge = "-"
				if ann["gnomADgc_AF_nfe"] != "":
					gng = ann["gnomADgc_AF_nfe"].split(',')[0].split('&')[0]
				else:
					gng = "-"
				if ann["gnomADec_AF_nfe"] != "":
					gne = ann["gnomADec_AF_nfe"].split(',')[0].split('&')[0]
				else:
					gne = "-"
				if dbe != '-':
					gnomeAF = float(dbe)
				elif ge != "-" and ge != '.':
					gnomeAF = float(ge.split('&')[0])
				if dbne != '-': 
					gnomeAF_nfe = float(dbne)
				elif gne != '-': 
					gnomeAF_nfe = float(gne.split('&')[0])
				if dbg != '-':
					gnomgAF = float(dbg)
				elif gg != '-':
					gnomgAF = float(gg.split('&')[0])
				if dbng != '-':
					gnomgAF_nfe = float(dbng)
				elif gng != '-':
					gnomgAF_nfe = float(gng)
				#print(mutated_indiv)
				if len(mutated_indiv)>0:
					print(counter, gene, mutid, chg, mut, atod, ppchg, rsid, top_csq, csq, impact, pphen, sift, cadd, csq2, lof_annotations)
					s = gene + '\t' + mutid + '\t' +  chg + '\t' + mut + '\t' + atod + '\t' +  ppchg + '\t' + rsid + '\t' + top_csq + '\t' + csq + '\t' + impact + '\t' + pphen + '\t' + sift + '\t' + str(cadd) + '\t' + csq2 + '\t' + str(gnomgAF) + '\t' + str(gnomgAF_nfe) + '\t' + str(gnomeAF) + '\t' + str(gnomeAF_nfe) + '\t' + str(ac) + '\t' + str(an) + '\t' + str(af) + '\t' + str(ac_case) + '\t' + str(an_case) + '\t' + str(ac_control) + '\t' + str(an_control) +'\t'+';'.join([str(v) for v in mutated_indiv]) + '\n'
					of.write(s)
			l = fp.readline().strip()
	print(counter)	
	of.close()
	return 0

    
def count_table(max_count, Ngene, case_m, control_m):
    max_count = int(max_count)
    case_count = np.zeros((Ngene, max_count))
    con_count = np.zeros((Ngene, max_count))
    for i in range(Ngene):
        z1 = Counter(case_m[i,])
        z2 = Counter(control_m[i,])
        for j in range(max_count):
            case_count[i,j] = z1[j+1]
            con_count[i,j] = z2[j+1]
	
    return(case_count, con_count)
	
	
def write_files(out_prefix, category, case_count, con_count, max_count, genes, cases, controls, max_cutoff):
	z1 = pd.DataFrame(case_count, index=genes, columns=[v for v in range(max_count)])
	z2 = pd.DataFrame(con_count, index=genes, columns=[v for v in range(max_count)])
	z1.to_csv(out_prefix+'_case_freqtable_'+category+'_'+str(max_cutoff)+'.tsv', sep='\t', header=False)
	z2.to_csv(out_prefix+'_control_freqtable_'+category+'_'+str(max_cutoff)+'.tsv', sep='\t', header=False)
	
def write_mut_files(out_prefix, category, case_mut, con_mut, genes, cases, controls, max_cutoff):
	z1 = pd.DataFrame(case_mut, index=genes, columns=cases)
	z2 = pd.DataFrame(con_mut, index=genes, columns=controls)
	z1.to_csv(out_prefix+'_case_mutations_'+category+'_'+str(max_cutoff)+'.tsv', sep='\t')
	z2.to_csv(out_prefix+'_control_mutations_'+category+'_'+str(max_cutoff)+'.tsv', sep='\t')


def write_pickle(out_prefix, category, case_count, con_count, max_count, genes, cases, controls, max_cutoff):
	z1 = pd.DataFrame(case_count, index=genes, columns=[v for v in range(max_count)])
	z2 = pd.DataFrame(con_count, index=genes, columns=[v for v in range(max_count)])
	z1_of = out_prefix+'_case_freqtable_'+category+'_'+str(max_cutoff)+'.pkl'
	z2_of = out_prefix+'_control_freqtable_'+category+'_'+str(max_cutoff)+'.pkl'
	with open(z1_of, 'wb') as f:
		pickle.dump(z1, f)
	with open(z2_of, 'wb') as f:
		pickle.dump(z2, f)


def write_mut_pickle(out_prefix, category, case_mut, con_mut, genes, cases, controls, max_cutoff):
	z1 = pd.DataFrame(case_mut, index=genes, columns=cases)
	z2 = pd.DataFrame(con_mut, index=genes, columns=controls)
	z1_of = out_prefix+'_case_mutations_'+category+'_'+str(max_cutoff)+'.pkl'
	z2_of = out_prefix+'_control_mutations_'+category+'_'+str(max_cutoff)+'.pkl'
	with open(z1_of, 'wb') as f:
		pickle.dump(z1, f)
	with open(z2_of, 'wb') as f:
		pickle.dump(z2, f)
	

def create_freq_table(mutfile, famfile, genefile, min_cutoff, max_cutoff, out_prefix):
	'''
	This program prepares case-control mutation counts tables from mutations file representation of 
	input vcf file.
	Input:
	- tab-separated mutations file
	- prefix of output ftable files
	Output: 
	- for each variant category two files containing case and control mutation counts respectively 
	'''
	
	#lines = [l.strip() for l in open(mutfile)]
	mutdf = pd.read_csv(mutfile,sep='\t',low_memory=False)
	genes = [l.strip().split('\t')[0] for l in open(genefile)]
	Ngene = len(genes)
	seldf = mutdf[mutdf['#Gene'].isin(genes)]

	fl = [l.strip() for l in open(famfile)]
	cases = []
	controls  = []
	for l in fl:
		if "\t" in l:
			x = l.split('\t')
		else:
			x = l.split()
		if x[5] == '1':
			controls.append(x[1])
		elif x[5] == '2':
			cases.append(x[1])

	print(len(cases), len(controls))

	Ncase = len(cases)
	Ncontrol = len(controls)
	case_dam = np.zeros((Ngene, Ncase))
	control_dam = np.zeros((Ngene, Ncontrol))
	case_dmis = np.zeros((Ngene, Ncase))
	control_dmis = np.zeros((Ngene, Ncontrol))
	case_lof = np.zeros((Ngene, Ncase))
	control_lof = np.zeros((Ngene, Ncontrol))
	case_mis = np.zeros((Ngene, Ncase))
	control_mis = np.zeros((Ngene, Ncontrol))
	case_syn = np.zeros((Ngene, Ncase))
	control_syn = np.zeros((Ngene, Ncontrol))
	case_nut = np.zeros((Ngene, Ncase))
	control_nut = np.zeros((Ngene, Ncontrol))
	
	dam_max_count = 0
	dmis_max_count = 0
	lof_max_count = 0
	mis_max_count = 0
	syn_max_count = 0
	nut_max_count = 0
	
	mcut = max_cutoff * 500000
	keeps = cases+controls
	for ind, row in seldf.iterrows():
		gene = row['#Gene']
		if gene not in genes:
			print(gene, "not in genelist")
			continue
		minds = row['mutated_individuals'].split(';')
		if len(minds) >= mcut:
			print(row["mutid"], "too_common")
			continue
		keepminds = [v for v in minds if v in keeps]
		mcase = [v for v in keepminds if v in cases]
		mcon = [v for v in keepminds if v in controls]
		af = len(mcase+mcon)/(2.0*len(cases+controls))
		ggaf = np.max(row['gg_af'])
		geaf = np.max(row['ge_af'])
		if ggaf != "-" and ggaf !="":
			ggaf = float(ggaf)
		else:
			ggaf = -1
		if geaf != "-" and geaf !="":
			geaf = float(geaf)
		else:
			geaf = -1
		print(row["mutid"], len(mcase+mcon), ggaf, geaf)
		if (af > min_cutoff and af <= max_cutoff) or (geaf > min_cutoff and geaf <= max_cutoff) or (ggaf > min_cutoff and ggaf <= max_cutoff):
			print(gene,row['csq2'],len(minds), str(af))
			if af == min_cutoff:
				continue
			r = genes.index(gene)
			if row['csq2'] == 'damaging':
				for v in mcase:
					c = cases.index(v)
					case_dam[r,c] += 1
					if case_dam[r,c] > dam_max_count:
						dam_max_count = case_dam[r,c]
				for v in mcon:
					c = controls.index(v)
					control_dam[r,c] += 1
					if control_dam[r,c] > dam_max_count:
						dam_max_count = control_dam[r,c]
			elif row['csq2'] == 'synonymous':
				for v in mcase:
					c = cases.index(v)
					case_syn[r,c] += 1
					if case_syn[r,c] > syn_max_count:
						syn_max_count = case_syn[r,c]
				for v in mcon:
					c = controls.index(v)
					control_syn[r,c] += 1
					if control_syn[r,c] > syn_max_count:
						syn_max_count = control_syn[r,c]
			elif row['csq2'] == 'neutral':
				for v in mcase:
					c = cases.index(v)
					case_nut[r,c] += 1
					if case_nut[r,c] > nut_max_count:
						nut_max_count = case_nut[r,c]
				for v in mcon:
					c = controls.index(v)
					control_nut[r,c] += 1
					if control_nut[r,c] > nut_max_count:
						nut_max_count = control_nut[r,c]
			if 'missense' in row['top_csq']:
				for v in mcase:
					c = cases.index(v)
					case_mis[r,c] += 1
					if case_mis[r,c] > mis_max_count:
						mis_max_count = case_mis[r,c]
				for v in mcon:
					c = controls.index(v)
					control_mis[r,c] += 1
					if control_mis[r,c] > mis_max_count:
						mis_max_count = control_mis[r,c]
				if row['csq2'] == 'damaging':
					for v in mcase:
						c = cases.index(v)
						case_dmis[r,c] += 1
						if case_dmis[r,c] > dmis_max_count:
							dmis_max_count = case_dmis[r,c]
					for v in mcon:
						c = controls.index(v)
						control_dmis[r,c] += 1
						if control_dmis[r,c] > dmis_max_count:
							dmis_max_count = control_dmis[r,c]
			elif 'synonymous' not in row['top_csq'] and 'missense' not in row['top_csq']:
				if row['csq2'] == 'damaging':
					for v in mcase:
						c = cases.index(v)
						case_lof[r,c] += 1
						if case_lof[r,c] > lof_max_count:
							lof_max_count = case_lof[r,c]
					for v in mcon:
						c = controls.index(v)
						control_lof[r,c] += 1
						if control_lof[r,c] > lof_max_count:
							lof_max_count = control_lof[r,c]


	dam_max_count = int(dam_max_count)
	case_dam_count, con_dam_count = count_table(dam_max_count, Ngene, case_dam, control_dam)
	dmis_max_count = int(dmis_max_count)
	case_dmis_count, con_dmis_count = count_table(dmis_max_count, Ngene, case_dmis, control_dmis)
	lof_max_count = int(lof_max_count)
	case_lof_count, con_lof_count = count_table(lof_max_count, Ngene, case_lof, control_lof)
	mis_max_count = int(mis_max_count)
	case_mis_count, con_mis_count = count_table(mis_max_count, Ngene, case_mis, control_mis)
	syn_max_count = int(syn_max_count)
	case_syn_count, con_syn_count = count_table(syn_max_count, Ngene, case_syn, control_syn)
	nut_max_count = int(nut_max_count)
	case_nut_count, con_nut_count = count_table(nut_max_count, Ngene, case_nut, control_nut)
	
	write_files(out_prefix, 'damaging', case_dam_count, con_dam_count, dam_max_count, genes, cases, controls, max_cutoff)
	#write_mut_files(out_prefix, 'damaging', case_dam, control_dam, genes, cases, controls, max_cutoff)
	write_files(out_prefix, 'damaging_missense', case_dmis_count, con_dmis_count, dmis_max_count, genes, cases, controls, max_cutoff)
	#write_mut_files(out_prefix, 'damaging_missense', case_dmis, control_dmis, genes, cases, controls, max_cutoff)
	write_files(out_prefix, 'LoF', case_lof_count, con_lof_count, lof_max_count, genes, cases, controls, max_cutoff)
	#write_mut_files(out_prefix, 'LoF', case_lof, control_lof, genes, cases, controls, max_cutoff)
	write_files(out_prefix, 'missense', case_mis_count, con_mis_count, mis_max_count, genes, cases, controls, max_cutoff)
	#write_mut_files(out_prefix, 'missense', case_mis, control_mis, genes, cases, controls, max_cutoff)
	write_files(out_prefix, 'synonymous', case_syn_count, con_syn_count, syn_max_count, genes, cases, controls, max_cutoff)
	#write_mut_files(out_prefix, 'synonymous', case_syn, control_syn, genes, cases, controls, max_cutoff)
	write_files(out_prefix, 'neutral', case_nut_count, con_nut_count, nut_max_count, genes, cases, controls, max_cutoff)
	#write_mut_files(out_prefix, 'neutral', case_nut, control_nut, genes, cases, controls, max_cutoff)	
	return 0


def generate_network_phy(infile, resource_dir, resource_prefix, outfile):
	'''
	This program creates the adjacency matrix for given gene list based on physical and genetic interactions from databases.

 	Input:
	- gene_file: list of genes with one gene symbol per line.
	- resource_dir: directory containing the .pkl files of 
    protein interaction databases.
    - resource_prefix (string): Prefix for database files
	- out_file: path to output file
	Output: 
	- a .tsv file containing the network adjacency matrix
	'''
	sel_genes = [l.strip() for l in open(infile)]
	Ngene = len(sel_genes)
	mat = np.zeros((Ngene, Ngene))
	np.fill_diagonal(mat, 2)
	
	#read databases
	res_pre = os.path.join(resource_dir,resource_prefix)
	string_db = pd.read_pickle(res_pre+'_'+"string_0.7.pkl")
	inbio_db = pd.read_pickle(res_pre+'_'+"inbio_sym.pkl")
	huri_db = pd.read_pickle(res_pre+'_'+"huri_sym.pkl")
	sl_db = pd.read_pickle(res_pre+'_'+"synlethal.pkl")
	cell_systems_int_db = pd.read_pickle(res_pre+'_'+"asyn_tnet.pkl")
	ab_int_db = pd.read_pickle(res_pre+'_'+"ab_oe.pkl")
	tdp43_int_db = pd.read_pickle(res_pre+'_'+"tdp_oe.pkl")
	
	#read_mapper_alias
	with open(os.path.join(resource_dir,"gene_aliases.txt")) as fp:
		dat = fp.read()
		
	alias_m = json.loads(dat)
	
	x = []
	for v in sel_genes:
		if v in alias_m.values():
			x.append(v)
		elif v in alias_m.keys():
			x.append(alias_m[v])
		else:
			x.append(v)

	mm = {}
	for idx, v in enumerate(x):
		mm[v] = idx
				
	a = string_db.loc[string_db['protein1'].isin(x) & string_db['protein2'].isin(x)].loc[:,['protein1','protein2']]
	b = huri_db.loc[huri_db['protein1'].isin(x) & huri_db['protein2'].isin(x)].loc[:,['protein1','protein2']]
	c = inbio_db.loc[inbio_db['protein1'].isin(x) & inbio_db['protein2'].isin(x)].loc[:,['protein1','protein2']]
	d = sl_db[sl_db['protein1'].isin(x) & sl_db['protein2'].isin(x)].loc[:,['protein1','protein2']]
	e = cell_systems_int_db[cell_systems_int_db['protein1'].isin(x) & cell_systems_int_db['protein2'].isin(x)].loc[:,['protein1','protein2']]
	f = ab_int_db[ab_int_db['protein1'].isin(x) & ab_int_db['protein2'].isin(x)].loc[:,['protein1','protein2']]
	g = tdp43_int_db[tdp43_int_db['protein1'].isin(x) & tdp43_int_db['protein2'].isin(x)].loc[:,['protein1','protein2']]
	
	cdf = pd.concat([a,b,c,d,e,f,g])
	for v in np.array(cdf):
		id1 = mm[v[0]]
		id2 = mm[v[1]]
		mat[id1,id2] = mat[id2,id1] = 1
		

	mat = pd.DataFrame(mat).astype('int')
	mat.index = sel_genes
	mat.columns = sel_genes
	mat.to_csv(outfile, sep='\t', header=False)
	

def generate_network_coexpression(gctfile, genefile, outfile):
	'''
	This program prepares the adjacency matrix for a given gene list 
    based on co-expression in GTEx tissue

	Input:
	- gctfile: Path to the TPM count matrix in .gct.gz file format from GTEx
	- gene_list_file: list of genes in the network of interest with 
    one gene symbol per line
	- outfile: path to the output file
	Output: 
	- a .tsv file containing the adjacency matrix of the network 
    representing pairwise co-expression of genes
	'''
	df = pd.read_csv(gctfile, sep='\t', skiprows=2, index_col=0)
	genelist = [l.strip() for l in open(genefile)]
	df2 = df[df["Description"].isin(genelist)]
	df3 = df2.groupby('Description').agg('mean',all)
	found = list(df3.index)
	notin = list(set(genelist) - set(found))
	inlist = list(set(genelist)-set(notin))
	corrz = df3.T.corr(method='pearson')
	y = pd.DataFrame(np.eye(len(genelist)))
	y.index = y.columns = genelist
	y.loc[found,found] = corrz
	y.to_csv(outfile, sep='\t', header=False)
    

def generate_network_coessentiality(infile, genefile, outfile):
	'''
	This program prepares the adjacency matrix for a given gene list 
    based on co-essentiality in DepMap cell lines.

	Input:
	- gene_dependency_file: DepMap gene dependency file subsetted to 
    only the cell lines of interest with column 0 being the ModelID 
    and rest of the columns being the genes in DepMap
	- gene_list_file: list of genes in the network of interest 
    with one gene symbol per line
	- outfile: path to the output file
	Output: 
	- a .tsv file containing the adjacency matrix of the network 
    representing pairwise co-essentiality of genes.
	'''

	y = pd.read_csv(infile, sep='\t', low_memory=False, index_col=0)
	colnames = y.columns
	genelist = [l.strip() for l in open(genefile)]
	notin = list(set(genelist) - set(colnames))
	inlist = list(set(genelist)-set(notin))
	z = y.loc[:,inlist]
	dimz = z.shape[1]
	for g in notin:
		z[g] = 0
	corrz = z.corr(method='spearman')
	df = pd.DataFrame(np.eye(len(genelist)))
	df.index = df.columns = genelist
	df.loc[inlist,inlist] = corrz
	print(df)
	df.to_csv(outfile, sep='\t', header=False)


def draw_colored_graph_phy(netfile, alphfile, out_fig, alternative, dim, layout):
	if alternative == "two-sided":
		#cmap = plt.cm.coolwarm
		cmap = sns.color_palette("vlag", as_cmap=True)
	elif alternative == "one-sided":
		cmap = plt.cm.Reds_r
	else:
		cmap = None
	
	df = pd.read_csv(netfile, sep='\t', header=None, index_col=0)
	df.columns = df.index
	print(df)
	x = np.matrix(df)
	for i in range(x.shape[0]):
		x[i,i] = 0
	G = nx.DiGraph(x)
	m = dict(zip(G,df.index))
	G = nx.relabel_nodes(G,m)
    # assumes first line is header
	lines = [l.strip() for l in open(alphfile)]
	m2  = {}
	for l in lines[1:]:
		x= l.split()
		m2[x[0]] = mycmap[int(x[2])]
	
	plt.figure(figsize=(dim,dim))
	if layout == 1:
		pos = nx.circular_layout(G)
	elif layout == 2:
		pos = nx.kamada_kawai_layout(G)
	elif layout == 3:
		pos = nx.random_layout(G)
	elif layout == 4:
		pos = nx.spectral_layout(G)
	elif layout == 5:
		pos = nx.spiral_layout(G)
	pos2 = dict()
	for v in pos.keys():
		pos2[v] = (pos[v][0]+0.02,pos[v][1]+0.05)
	nx.draw(G, pos, with_labels=False, node_size=100, node_color=list(m2.values()), edge_color=(0.8,0.8,0.8), arrows=False)
	nx.draw_networkx_labels(G, pos=pos2)
	plt.savefig(out_fig, dpi=360, format='png', bbox_inches='tight')


def edge_color(v):
	if v <0:
		return '#436EEE'
	elif v>0 and v!=1:
		return '#EE5C42'
	else:
		return "#AAAAAA"

def draw_colored_graph_val(netfile, alphfile, out_fig, alternative, onedgelabel, dim):
	if alternative == "two-sided":
		cmap = plt.cm.bwr
	elif alternative == "one-sided":
		cmap = plt.cm.Reds_r
	else:
		cmap = None
	
	df = pd.read_csv(netfile, sep='\t', header=None, index_col=0)
	df.columns = df.index
	#sdf = df.sort_index()
	x = np.matrix(df)
	for i in range(x.shape[0]):
		x[i,i] = 0
	G = nx.DiGraph(x)
	m = dict(zip(G,df.index))
	G = nx.relabel_nodes(G,m)
	lines = [l.strip() for l in open(alphfile)]
	m2  = {}
    # assumes first line is header
	for l in lines[1:]:
		x= l.split('\t')
		print(x[0],x[2], mycmap[int(x[2])])
		m2[x[0]] = mycmap[int(x[2])]

	print(m2)
	edge_labels = dict([((u,v,), f"{d['weight']:.2f}") for u,v,d in G.edges(data=True)])
	weights = [G[u][v]['weight'] for u,v,d in G.edges(data=True)]
	ecolors = [edge_color(G[u][v]['weight']) for u,v,d in G.edges(data=True)]

	plt.figure(figsize=(dim,dim))
	pos = nx.circular_layout(G)
	pos2 = dict()
	for v in pos.keys():
		pos2[v] = (pos[v][0]+0.05,pos[v][1]+0.05)
	print(G.nodes)
	nc = [m2[n] for n in G.nodes]
	nx.draw(G, pos, with_labels=False, node_size=100, node_color=nc, edge_color=ecolors, 
width=weights,arrows=False)
	nx.draw_networkx_labels(G, pos=pos2)
	if onedgelabel==1:
		nx.draw_networkx_edge_labels(G,pos,edge_labels)
	plt.savefig(out_fig, dpi=360, format='png',bbox_inches='tight')


def subset_coding_mutations(infile, outfile):
	#headline = [v.strip() for v in open(headerfile)]
	incl = ['frameshift_variant','inframe_deletion','inframe_insertion','missense_variant','protein_altering_variant','splice_acceptor_variant','splice_donor_variant','splice_region_variant','start_lost','stop_gained','stop_lost','synonymous_variant']
	x = pd.read_csv(infile, sep='\t', low_memory=False)
	x_coding = x[x['top_csq'].isin(incl)]
	x_coding.to_csv(outfile,sep='\t',index=False,header=False)
