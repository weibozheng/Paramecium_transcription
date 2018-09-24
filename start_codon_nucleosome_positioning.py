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
	cell=[0]*20
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
	
	if cell[0]>=10: #This is the limit of length
		o3=ppp[cell[0]]-int(cell[0]/2)
		return(o3,cell[0])
	else:
		return('NA','NA')
def find_cloest_summit(chrs,sc_pos,types,this_gene):
	#print (chr,sc_pos,type)
	np_limit=10
	average_gene_length=1000
	if (chrs in sc_len_dict) and (sc_pos in sc_len_dict[chrs]):
		added_str1=str(sc_len_dict[chrs][sc_pos]) + "\t" + sc_type_dict[chrs][sc_pos]
	else:
		added_str1="NULL\tNULL"
	if (chrs in sc_pos_motif_dict) and (sc_pos in sc_pos_motif_dict[chrs]):
		added_str2=str(sc_pos_motif_dict[chrs][sc_pos]) + "\t" + sc_type_motif_dict[chrs][sc_pos]
	else:
		added_str2="NULL\tNULL"
	if (chrs in sc_TATA_length_dict) and (sc_pos in sc_TATA_length_dict[chrs]):
		added_str3=str(sc_TATA_pos_dict[chrs][sc_pos]) + "\t" + str(sc_TATA_length_dict[chrs][sc_pos])
	else:
		added_str3="NULL\tNULL"
	if types=='+':
		smallest_v=2000000
		for i in range (0,50):
			tss_pos=str(sc_pos-i)
			if chrs in mpile:
				if tss_pos in mpile[chrs]:
					if smallest_v>mpile[chrs][tss_pos]:
						smallest_v=mpile[chrs][tss_pos]
						smallest_pos=i
				else:
					smallest_pos=i
					break
			else:
				smallest_pos=i
				break
	if types=='-':
		smallest_v=2000000
		for i in range (0,50):
			tss_pos=str(sc_pos+i)
			if chrs in mpile:
				if tss_pos in mpile[chrs]:
					if smallest_v>mpile[chrs][tss_pos]:
						smallest_v=mpile[chrs][tss_pos]
						smallest_pos=i
				else:
					smallest_pos=i
					break
			else:
				smallest_pos=i
				break
	if types=='+':
		np=0
		for i in range(1,average_gene_length):
			p_a=str(sc_pos-i)
			if (chrs in smt_v_dict) and (p_a in smt_v_dict[chrs]):
				np+=1
				print ("%s\t%d\t%d\t-%d\t%s\t%s\t-%s\t%s\t%s\t%s\t%s\t%s\t+\t%s\t%s" % (chrs,sc_pos,smallest_pos,np,p_a,smt_len_dict[chrs][p_a],i,smt_v_dict[chrs][p_a],smt_f_dict[chrs][p_a],added_str1,added_str2,exp_dict[this_gene],this_gene,added_str3),file=output1)
				if np>=np_limit:
					break
		np=0
		for i in range(average_gene_length):
			p_a=str(sc_pos+i)
			if (chrs in smt_v_dict) and (p_a in smt_v_dict[chrs]):
				np+=1
				print ("%s\t%d\t%d\t+%d\t%s\t%s\t+%s\t%s\t%s\t%s\t%s\t%s\t+\t%s\t%s" % (chrs,sc_pos,smallest_pos,np,p_a,smt_len_dict[chrs][p_a],i,smt_v_dict[chrs][p_a],smt_f_dict[chrs][p_a],added_str1,added_str2,exp_dict[this_gene],this_gene,added_str3),file=output1)
				if np>=np_limit:
					break
	if types=='-':
		np=0
		for i in range(1,average_gene_length):
			p_a=str(sc_pos+i)
			if (chrs in smt_v_dict) and (p_a in smt_v_dict[chrs]):
				np+=1
				print ("%s\t%d\t%d\t-%d\t%s\t%s\t-%s\t%s\t%s\t%s\t%s\t%s\t-\t%s\t%s" % (chrs,sc_pos,smallest_pos,np,p_a,smt_len_dict[chrs][p_a],i,smt_v_dict[chrs][p_a],smt_f_dict[chrs][p_a],added_str1,added_str2,exp_dict[this_gene],this_gene,added_str3),file=output1)
				if np>=np_limit:
					break
		np=0
		for i in range(average_gene_length):
			p_a=str(sc_pos-i)
			if (chrs in smt_v_dict) and (p_a in smt_v_dict[chrs]):
				np+=1
				print ("%s\t%d\t%d\t+%d\t%s\t%s\t+%s\t%s\t%s\t%s\t%s\t%s\t-\t%s\t%s" % (chrs,sc_pos,smallest_pos,np,p_a,smt_len_dict[chrs][p_a],i,smt_v_dict[chrs][p_a],smt_f_dict[chrs][p_a],added_str1,added_str2,exp_dict[this_gene],this_gene,added_str3),file=output1)
				if np>=np_limit:
					break
xls_inf=open(r"C:\Users\wbzhe\Paramecium2\ptet\pooled\N_dc2_projects_ciliate_Paramecium_tetraurelia_Sample_4V_ptet_mnase.smooth.positions.xls")#nucleosome_position_produced_by_danpos
output1=open(r"C:\Users\wbzhe\Paramecium2\ptet\pooled\start_codon_summit.txt","w")#table_format_output_file
print ("chr\tpos_start_codon\tTSS\tserial_nucleo\tpos_nucleo\tlength_nucleo\trelative_distance\tsmt_value\tsmt_fuzzines\tintergenic_length\tintergenic_type	motif_pos\tmotif_type\texp\tstrand\tgene\tTATA_pos\tTATA_length",file=output1)
pro_inf=open(r"C:\Users\wbzhe\Paramecium2\ptet\trans_veg_bowtie2_ptet.profile") #profile input
mpile=dict()
for line in pro_inf:
	mobj=re.match(">.*scaffold51_(\d+)",line)
	if mobj:
		this_tag=mobj.group(1)
	else:
		mobj_1=re.match("(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)",line)
		
		if mobj_1:
			pos=mobj_1.group(1)
			#print (pos)
			A=int(mobj_1.group(2))
			C=int(mobj_1.group(3))
			G=int(mobj_1.group(4))
			T=int(mobj_1.group(5))
			sums=A+C+G+T
			addword2dict(mpile,this_tag,pos,sums)
pro_inf.close()
smt_v_dict=dict()
smt_f_dict=dict()
smt_len_dict=dict()
for line in xls_inf:
	mobj=re.match(".*sca.*?_(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+\.\d+)\t(\d+\.\d+)",line)
	if mobj:
		chr=mobj.group(1)
		smt_len=str(int(mobj.group(3))-int(mobj.group(2)))
		smt=mobj.group(4)
		smt_v=mobj.group(5)
		smt_f=mobj.group(6)
		addword2dict(smt_v_dict,chr,smt,smt_v)
		addword2dict(smt_f_dict,chr,smt,smt_f)
		addword2dict(smt_len_dict,chr,smt,smt_len)
gff_inf=open(r"C:\Users\wbzhe\Paramecium\annotations\ptetraurelia.gff3");
sc_len_dict=dict()
sc_type_dict=dict()
sc_pos_motif_dict=dict()
sc_type_motif_dict=dict()
sc_TATA_pos_dict=dict()
sc_TATA_length_dict=dict()
for rec in SeqIO.parse(r"C:\Users\wbzhe\Paramecium2\ptet\ptet_intergenic.txt","fasta"):#intergenic_regions
	tag=rec.id
	seq=str(rec.seq)
	length=len(seq)
	mobj=re.match("sca.*?_(\d+)\|.*?\|(..)\|(\d+)\.\.(\d+)",tag)#>scaffold_0001|ig4|-+|10218..10228
	chr=mobj.group(1)
	ig_type=mobj.group(2)
	start=int(mobj.group(3))
	end=int(mobj.group(4))
	if ig_type=='-+':
		sc1=start-1
		sc2=end+1
		(pos_TATA_L,length_TATA_L)=obtain_pos_TATA(seq[:50])
		(pos_TATA_R,length_TATA_R)=obtain_pos_TATA(rev_com(seq[length-50:]))
		addword2dict(sc_len_dict,chr,sc1,length)
		addword2dict(sc_type_dict,chr,sc1,'2')
		addword2dict(sc_len_dict,chr,sc2,length)
		addword2dict(sc_type_dict,chr,sc2,'2')
		addword2dict(sc_TATA_pos_dict,chr,sc1,pos_TATA_L)
		addword2dict(sc_TATA_length_dict,chr,sc1,length_TATA_L)
		addword2dict(sc_TATA_pos_dict,chr,sc2,pos_TATA_R)
		addword2dict(sc_TATA_length_dict,chr,sc2,length_TATA_R)
		mobj1=re.match("(.*?)(.....AAATCTTT.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=17:
				len_up_seq=len(mobj1.group(1))
				pos_motif=len_up_seq+10
				addword2dict(sc_pos_motif_dict,chr,sc1,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc1,'AAAGATTT')
		mobj1=re.match("(.*?)(.....AAAGATTT.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=17:
				len_up_seq=len(mobj1.group(1))
				pos_motif=len_up_seq+10
				addword2dict(sc_pos_motif_dict,chr,sc1,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc1,'AAATCTTT')
		mobj2=re.match(".*(.......AAATCTTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=17:
				len_up_seq=len(mobj2.group(2))
				pos_motif=len_up_seq+11
				addword2dict(sc_pos_motif_dict,chr,sc2,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc2,'AAATCTTT')
		mobj2=re.match(".*(.......AAAGATTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=17:
				len_up_seq=len(mobj2.group(2))
				pos_motif=len_up_seq+11
				addword2dict(sc_pos_motif_dict,chr,sc2,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc2,'AAAGATTT')
		mobj1=re.match("(.*?)(.....AAAAACG.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=15:
				len_up_seq=len(mobj1.group(1))
				pos_motif=len_up_seq+10
				addword2dict(sc_pos_motif_dict,chr,sc1,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc1,'CGTTTTT')
		mobj2=re.match(".*(.......CGTTTTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=15:
				len_up_seq=len(mobj2.group(2))
				pos_motif=len_up_seq+10
				addword2dict(sc_pos_motif_dict,chr,sc2,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc2,'CGTTTTT')
	if ig_type=='++':
		sc=end+1
		(pos_TATA,length_TATA)=obtain_pos_TATA(rev_com(seq[length-50:]))
		addword2dict(sc_TATA_pos_dict,chr,sc,pos_TATA)
		addword2dict(sc_TATA_length_dict,chr,sc,length_TATA)
		addword2dict(sc_len_dict,chr,sc,length)
		addword2dict(sc_type_dict,chr,sc,'1')
		mobj2=re.match(".*(.......AAATCTTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=17:
				len_up_seq=len(mobj2.group(2))
				pos_motif=len_up_seq+11
				addword2dict(sc_pos_motif_dict,chr,sc,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc,'AAATCTTT')
		mobj2=re.match(".*(.......AAAGATTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=17:
				len_up_seq=len(mobj2.group(2))
				pos_motif=len_up_seq+11
				addword2dict(sc_pos_motif_dict,chr,sc,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc,'AAAGATTT')
		mobj2=re.match(".*(.......CGTTTTT.....)(.*)",seq[length-50:])
		if mobj2:
			AT_num=len(re.findall("A|T",mobj2.group(1)))
			if AT_num>=15:
				len_up_seq=len(mobj2.group(2))
				pos_motif=len_up_seq+10
				addword2dict(sc_pos_motif_dict,chr,sc,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc,'CGTTTTT')
	if ig_type=='--':
		sc=start-1
		(pos_TATA,length_TATA)=obtain_pos_TATA(seq[:50])
		addword2dict(sc_TATA_pos_dict,chr,sc,pos_TATA)
		addword2dict(sc_TATA_length_dict,chr,sc,length_TATA)
		addword2dict(sc_len_dict,chr,sc,length)
		addword2dict(sc_type_dict,chr,sc,'1')
		mobj1=re.match("(.*?)(.....AAATCTTT.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=17:
				len_up_seq=len(mobj1.group(1))
				pos_motif=len_up_seq+10
				addword2dict(sc_pos_motif_dict,chr,sc,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc,'AAAGATTT')
		mobj1=re.match("(.*?)(.....AAAGATTT.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=17:
				len_up_seq=len(mobj1.group(1))
				pos_motif=len_up_seq+10
				addword2dict(sc_pos_motif_dict,chr,sc,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc,'AAATCTTT')
		mobj1=re.match("(.*?)(.....AAAAACG.......)",seq[:50])
		if mobj1:
			AT_num=len(re.findall("A|T",mobj1.group(2)))
			if AT_num>=15:
				len_up_seq=len(mobj1.group(1))
				pos_motif=len_up_seq+10
				addword2dict(sc_pos_motif_dict,chr,sc,pos_motif)
				addword2dict(sc_type_motif_dict,chr,sc,'CGTTTTT')
exp_inf=open(r"C:\Users\wbzhe\Paramecium2\ptet\pooled\CDS_rsem_out.genes.results")#expression_data
exp_dict=dict()
for line in exp_inf:
	mobj=re.match("(.*?)\t.*?\t.*?\t.*?\t.*?\t(.*?)\t.*",line)
	tag=mobj.group(1)
	TPM=mobj.group(2)
	exp_dict.update({tag:TPM})
name_former='NO'
for line in gff_inf:
	mobj=re.match("(.*?)\t.*?\t(.*?)\t(.*?)\t(.*?)\t.*?\t(.*?)\t.*?\t(.*?)\n",line)
	name_this=mobj.group(1)
	tag_this=mobj.group(2)
	start_this=mobj.group(3)
	end_this=mobj.group(4)
	frame_this=mobj.group(5)
	att=mobj.group(6)
	mobj_x=re.search("Parent=(.*)",att)
	if mobj_x:
		gene=mobj_x.group(1)
	if (name_this == name_former):
		if (tag_this == 'CDS'):
			if (tag_former == 'CDS'):
				intergene_start=end_this
				intergene_type_f=frame_this
				gene_f=gene
				tag_former=tag_this;
			else:
				mobj2=re.match(".*?_(\d+)",name_this)
				chr=mobj2.group(1)
				intergene_end=start_this
				intergene_type_l=frame_this
				gene_l=gene
				if intergene_type_f == '-':
					start_codon_pos=int(intergene_start)
					find_cloest_summit(chr,start_codon_pos,intergene_type_f,gene_f)
				if intergene_type_l == '+':
					start_codon_pos=int(intergene_end)
					find_cloest_summit(chr,start_codon_pos,intergene_type_l,gene_l)
				intergene_start=end_this
				intergene_type_f=frame_this
				gene_f=gene
				tag_former=tag_this
		else:
			tag_former=tag_this
		name_former=name_this
	else:
		if (tag_this == 'CDS'):
			name_former=name_this
			intergene_start=end_this
			intergene_type_f=frame_this
			gene_f=gene
			tag_former=tag_this
gff_inf.close()

