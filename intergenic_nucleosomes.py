import re
from Bio import SeqIO
def addword2dict(thedict, key_a, key_b, val):
	if key_a in thedict:
		thedict[key_a].update({key_b: val})
	else:
		thedict.update({key_a:{key_b: val}})
def rev_com(seqa):
	seq_split=list(seqa)
	newseq=list()
	for value in seq_split:
		if value == 'A':
			newseq.append('T')
		elif value == 'T':
			newseq.append('A')
		elif value == 'C':
			newseq.append('G')
		elif value == 'G':
			newseq.append('C')
		else :
			newseq.append('N')
	seqa=''.join(str(newseq.reverse()))
	return(seqa)
def obtain_pos_TATA(this_seq):
	former='A'
	count=0
	cell=[0]*100
	ppp=dict()
	pos=0
	for this in this_seq:
		pos+=1
		if((this == 'A' or this == 'T') and (former == 'A' or former == 'T')):
			cell[count]+=1
			ppp.update({cell[count]:pos})
			former=this
		if((this == 'A' or this == 'T') and (former == 'C' or former == 'G')):
			count+=1
			cell[count]+=1
			former=this
		former=this
	cell.sort(reverse=True)
	if not cell[0] in num_longest:
		num_longest.update({cell[0]:1})
	else:
		num_longest[cell[0]]+=1
	if cell[0]>=10: #This is the limit of length
		o3=ppp[cell[0]]-int(cell[0]/2)
		if not o3 in pos_longest:
			pos_longest.update({o3:1})
		else:
			pos_longest[o3]+=1
		return(o3)
	else:
		return('NA')
xls_inf=open(r"C:\Users\wbzhe\Paramecium2\ptet\pooled\N_dc2_projects_ciliate_Paramecium_tetraurelia_Sample_4V_ptet_mnase.smooth.positions.xls")#nucleosome_position_produced_by_danpos
output1=open(r"C:\Users\wbzhe\Paramecium2\ptet\pooled\intergenic_nucleosomes.txt","w")#table_format_output_file
smt_v_dict=dict()
smt_f_dict=dict()
smt_len_dict=dict()
num_longest=dict()
pos_longest=dict()
for line in xls_inf:
	mobj=re.match(".*sca.*?_(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+\.\d+)\t(\d+\.\d+)",line)
	if mobj:
		chrs=mobj.group(1)
		smt_len=str(int(mobj.group(3))-int(mobj.group(2)))
		smt=mobj.group(4)
		smt_v=mobj.group(5)
		smt_f=mobj.group(6)
		addword2dict(smt_v_dict,chrs,smt,smt_v)
		addword2dict(smt_f_dict,chrs,smt,smt_f)
		addword2dict(smt_len_dict,chrs,smt,smt_len)
print ("tag\ttype\tstart\tend\tlength\tnucleosome_number\tnucleosome_pos\tnucleosome_smt_v\tnucleosome_smt_fuzz\tmotif_L\tL_pos\tmotif_R\tR_pos",file=output1)
for rec in SeqIO.parse(r"C:\Users\wbzhe\Paramecium2\ptet\ptet_intergenic.txt","fasta"):#intergenic_regions
	tag=rec.id
	seq=str(rec.seq)
	length=len(seq)
	mobj=re.match("sca.*?_(\d+)\|.*?\|(..)\|(\d+)\.\.(\d+)",tag)#>scaffold_0001|ig4|-+|10218..10228
	chrs=mobj.group(1)
	ig_type=mobj.group(2)
	start=int(mobj.group(3))
	end=int(mobj.group(4))
	has_motif=0
	has_motif1=0
	has_motif2=0
	nucleosome_number=0
	nucleosome_pos=''
	nucleosome_smt_v=''
	nucleosome_smt_fuzz=''
	if ig_type=='-+':
		pos_TATA_L=obtain_pos_TATA(seq[:100])
		pos_TATA_R=obtain_pos_TATA(rev_com(seq[length-100:]))
		mobj1=re.match("(.*?)(.....AAATCTTT.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=17:
				len_up_seq=len(mobj1.group(1))
				pos_motif1=len_up_seq+10
				motif_type1='AAAGATTT'
				has_motif1=1
		mobj1=re.match("(.*?)(.....AAAGATTT.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=17:
				len_up_seq=len(mobj1.group(1))
				pos_motif1=len_up_seq+10
				motif_type1='AAATCTTT'
				has_motif1=1
		mobj2=re.match(".*(.......AAATCTTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=17:
				len_up_seq=len(mobj2.group(2))
				pos_motif2=len_up_seq+11
				motif_type2='AAATCTTT'
				has_motif2=1
		mobj2=re.match(".*(.......AAAGATTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=17:
				len_up_seq=len(mobj2.group(2))
				pos_motif2=len_up_seq+11
				motif_type2='AAAGATTT'
				has_motif2=1
		mobj1=re.match("(.*?)(.....AAAAACG.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=15:
				len_up_seq=len(mobj1.group(1))
				pos_motif1=len_up_seq+10
				motif_type1='CGTTTTT'
				has_motif1=1
		mobj2=re.match(".*(.......CGTTTTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=15:
				len_up_seq=len(mobj2.group(2))
				pos_motif2=len_up_seq+10
				motif_type2='CGTTTTT'
				has_motif2=1
		if has_motif1==1:
			motif_L=motif_type1
			L_pos=pos_motif1
		else:
			motif_L='NULL'
			L_pos='NA'
		if has_motif2==1:
			motif_R=motif_type2
			R_pos=pos_motif2
		else:
			motif_R='NULL'
			R_pos='NA'
	elif ig_type=='++':
		mobj2=re.match(".*(.......AAATCTTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=17:
				len_up_seq=len(mobj2.group(2))
				pos_motif=len_up_seq+11
				motif_type='AAATCTTT'
				has_motif=1
		mobj2=re.match(".*(.......AAAGATTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=17:
				len_up_seq=len(mobj2.group(2))
				pos_motif=len_up_seq+11
				motif_type='AAAGATTT'
				has_motif=1
		mobj2=re.match(".*(.......CGTTTTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=15:
				len_up_seq=len(mobj2.group(2))
				pos_motif=len_up_seq+10
				motif_type='CGTTTTT'
				has_motif=1
		if has_motif==1:
			motif_R=motif_type
			R_pos=pos_motif
			motif_L='NULL'
			L_pos='NA'
		else:
			motif_R='NULL'
			R_pos='NA'
			motif_L='NULL'
			L_pos='NA'
	elif ig_type=='--':
		mobj1=re.match("(.*?)(.....AAATCTTT.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=17:
				len_up_seq=len(mobj1.group(1))
				pos_motif=len_up_seq+10
				motif_type='AAAGATTT'
				has_motif=1
		mobj1=re.match("(.*?)(.....AAAGATTT.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=17:
				len_up_seq=len(mobj1.group(1))
				pos_motif=len_up_seq+10
				motif_type='AAATCTTT'
				has_motif=1
		mobj1=re.match("(.*?)(.....AAAAACG.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=15:
				len_up_seq=len(mobj1.group(1))
				pos_motif=len_up_seq+10
				motif_type='CGTTTTT'
				has_motif=1
		if has_motif==1:
			motif_L=motif_type
			L_pos=pos_motif
			motif_R='NULL'
			R_pos='NA'
		else:
			motif_R='NULL'
			R_pos='NA'
			motif_L='NULL'
			L_pos='NA'
	#print (chrs)
	if chrs in smt_v_dict:
		for key_x in smt_v_dict[chrs]:
			vk=int(key_x)
			if vk>=start and vk<=end:
				nucleosome_number+=1
				nucleosome_pos+='|'+key_x
				nucleosome_smt_v+='|'+smt_v_dict[chrs][key_x]
				nucleosome_smt_fuzz+='|'+smt_f_dict[chrs][key_x]
			if nucleosome_number==0:
				nucleosome_number=0
				nucleosome_smt_fuzz='NULL'
				nucleosome_pos='NULL'
				nucleosome_smt_v='NULL'
		print ("%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (tag,ig_type,start,end,length,nucleosome_number,nucleosome_pos,nucleosome_smt_v,nucleosome_smt_fuzz,motif_L,L_pos,motif_R,R_pos),file=output1)