#!/usr/bin/env python2

import time
import optparse
from Bio import SeqIO
from os.path import exists
from subprocess import Popen, PIPE
from collections import defaultdict

################################# Command line options

desc='It takes as input the reference annotation in CAT gff3 format and parse the resulted annotations of the target species in CAT gff3 format to return a single copy orthology table'

parser = optparse.OptionParser(description=desc, version='%prog version 1 - 14-11-2022 - Author: FCicconardi')

parser.add_option('-d', '--Directory', dest='inDir', help='DeltaBLAST hits in outfmt6 of mikado_prepared.fasta. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-g', '--RefGff3', dest='rgff', help='GFF3 file of the reference annotation. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-s', '--Species', dest='sps', help='Comma separated species list to analyse. Mandatory opt.', action='store', metavar='FILE')

(opts, args) = parser.parse_args()

mandatories = ['inDir','rgff']
for m in mandatories:
        if not opts.__dict__[m]:
                print "\nWARNING! One or more options not specified\n"
                parser.print_help()
                exit(-1)

############################## Reading files and parametersfrom sys import argv

file=open(opts.rgff, 'r')
gff=file.readline()

GffOrder=defaultdict(list)
gffDCT=defaultdict(list)
Gorder=0
while gff:
	if gff.startswith('#'): gff=file.readline()
	else:
		el=gff.strip().split('\t')
		feat=el[2]
		if feat == 'gene':
			Gorder+=1
			scf=el[0]
			start=int(el[3])
			end=int(el[4])
			strand=el[6]
			gene=el[8].split(';')[0][3:].replace('gene:','')
			key='%s:%s..%s(%s)' % (scf,start,end,strand)
			GffOrder[Gorder].append(gene)
			gffDCT[gene].append(key)
		gff=file.readline()


OGdict=defaultdict(list)
GeneTransc=defaultdict(list)
PredToRef=defaultdict(list)
GeneToCoords=defaultdict(list)
Paralog=defaultdict(list)
ParalogToCatID=defaultdict(list)
CatIdTORefID=defaultdict(list)
tmapDCT=defaultdict(list)
Alias=defaultdict(list)

ls_output = Popen(['ls', '-1', opts.inDir], stdout=PIPE)
stdout, stderr = ls_output.communicate()
directory = stdout.decode('ascii').splitlines()

AllSpecies=opts.sps.strip().split(',')

#AllSpecies=['Goto','Mcon','Mpol']
#AllSpecies=['Goto','Mcon','Mpol','Herd','Hmel','Eisa','Smor','Diul']
#AllSpecies=['Herd','Hmel','Eisa','Smor','Diul']
#AllSpecies=['Smor','Pdid','Dpha','Diul','Ptel','Djun','Avfl','Avcr','Avpe','Haoe','Herd','Hmel']

NoGoodCC=['x','i','ri','p','u']

for sp in AllSpecies:

	if exists(sp+'.ReftoNovelLoci.tmap.alias'):
		alias=open(sp+'.ReftoNovelLoci.tmap.alias', 'r').readlines()
		for el in alias:
			el=el.strip().split()
			Alias[el[0]].append(el[1])

	gff='%s.gff3' % sp
	species=sp

	gff=open(gff, 'r').readlines()

	for el in gff:
		if el.startswith('#'): continue
		else:
			feat=el.strip().split()[2]
			info=el.strip().split()[8].split(';')
			if feat == 'gene':
				if 'source_gene=None' not in el:
					for i in info:
						if 'source_gene=' in i:
							OG=i.split('=')[1] # The Paralog A
						if 'gene_id=' in i:
							Gid=i.split('=')[1]
						if 'collapsed_gene_names=' in i:
							tmp=i.split('=')[1].split(',') # All other paralogs that map here
							for Par in tmp:
								Paralog[species].append(Par)
					if Gid.replace('_','.') in Alias:
						Gid=Alias[Gid.replace('_','.')][0]
					OGdict[OG].append(Gid)
					#print OG, Gid
				else: Gid='' 
			elif feat == 'transcript':
				if Gid != '':
					for i in info:
						if 'ID=' in i:
							Tid=i.split('=')[1]
							GeneTransc[Gid].append(Tid)
	tmap='Ref%stoCATrun+BlastP.Filtered.tmap' % sp
	refmap=open(tmap, 'r').readlines()

	for el in refmap:
		if el.startswith('ref_id'): continue
		else:
			Loc=el.split()[18]
			if el.split()[2] == 'NA': PredID=species+'_NA'
			else:
				ClassCode=el.split()[2]
				PredID=species+'_'+el.split()[3].split('.')[1]
				if ClassCode.lower() not in NoGoodCC: 
					RefID='.'.join(el.split()[0].replace('transcript:','').split('.')[:-1])
				else: RefID=''
				if RefID == '': RefID=PredID.replace('_','.')
				if PredID.replace('_','.') in Alias:
					PredID=Alias[PredID.replace('_','.')][0]
				if RefID.replace('_','.') in Alias:
					RefID=Alias[RefID.replace('_','.')][0] 
				tmapDCT[PredID].append(el.strip())
				PredToRef[PredID].append(RefID)
				GeneToCoords[PredID].append(Loc)
				#print PredID, RefID, Loc, ClassCode


ToExcude=defaultdict(list)

for sp in AllSpecies:
	for paralog in Paralog[sp]:
		if paralog in OGdict:
			for CATid in OGdict[paralog]:
				if CATid.startswith(sp):
					ToExcude[paralog].append(CATid)


Tsp=GffOrder[Gorder][0].split('.')[0]

print 'OG\t%s\t%s\tOccupancy' % (Tsp,'\t'.join(AllSpecies))
c=0

for pos in GffOrder:
	c+=1
	zeros=5-len(str(c))
	Z=zeros*'0'
	OG='OG_'+Z+str(c)
	og=GffOrder[pos][0]
	GeneCoords=gffDCT[og][0]
	
	if len(OGdict[og]) >= 1:
		SpToGene=defaultdict(list)
		for gene in OGdict[og]:
			sp=gene[:4]
			SpToGene[sp].append(gene)

		for sp in AllSpecies:
			if sp not in SpToGene:
				SpToGene[sp].append('-')

		AllIDs=''
		for sp in AllSpecies:
			PredID=SpToGene[sp][0]
		#	print PredID
			if og in ToExcude:
				if PredID in ToExcude[og]: PredID='-'

			if PredID in GeneToCoords:
				Loc=GeneToCoords[PredID][0]

			if PredID in PredToRef:
				if len(set(PredToRef[PredID])) > 1:
					BestMatch=defaultdict(list)
					for match in tmapDCT[PredID]:
						Bscore=match.split()[23]
						if Bscore == '-': Bscore=0.0
						else: Bscore=float(Bscore)
						BestMatch[Bscore].append(match)
					if max(BestMatch) > 0:
						BEST=BestMatch[max(BestMatch)][0].split()[1].replace('gene:','')
					else:
						BEST='-'
						Loc='-'
						#print sp,og,PredID, set(PredToRef[PredID]), len(set(PredToRef[PredID]))
						#for match in tmapDCT[PredID]:
						#	print match
						#print
					AllIDs+='\t'+BEST+'|'+Loc
				else: AllIDs+='\t'+PredToRef[PredID][0]+'|'+GeneToCoords[PredID][0]
			else: AllIDs+='\t'+'-|-'
		print '%s\t%s|%s\t%s\t%s' % (OG,og,GeneCoords,AllIDs.strip(),len(AllSpecies)+1-AllIDs.strip().count('-|-'))
	else:
		line='-|-\t'*len(AllSpecies)
		print '%s\t%s|%s\t%s\t1' % (OG,og,GeneCoords,line.strip())
		

