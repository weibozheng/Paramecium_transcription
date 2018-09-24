from Bio import SeqIO
import re
output=open(r"C:\Users\wbzhe\Paramecium2\ptet\intergenic_TATA_sliding_PP.v2.matrix","w")
store_dict=dict()
for rec in SeqIO.parse(r"C:\Users\wbzhe\Paramecium2\ptet\ptetraurelia_mac_51_renamed.fa","fasta"): #intergenic region input
	id=rec.id
	mobj=re.match("ptetraurelia_mac_51_scaffold51_(\d+)",id)
	chrs=mobj.group(1)
	store_dict.update({chrs:str(rec.seq)})
	

seqs_dict=SeqIO.to_dict(SeqIO.parse(r"C:\Users\wbzhe\Paramecium\annotations\++\ptetraurelia_++.txt","fasta")) #intergenic region input
seqs_dict=sorted(seqs_dict.items(),key=lambda x:(len(str(x[1].seq))))
for j in range(0,550):
	print ("\t%d" % (j),file=output,end="")
print ("",file=output)
for key in seqs_dict:
	this_tag=key[1].id
	mobj=re.match("scaffold51_(\d+)\|(.*?)\|(\d+)\.\.(\d+)",this_tag)
	chr=mobj.group(1)
	start=mobj.group(3)
	end=mobj.group(4)
	length_this=len(str(key[1].seq))
	tag=chr + "_" + mobj.group(2) + "_" + start + "_" + end + "_" + str(length_this);
	seq_this=store_dict[chr][int(end)-500:int(end)+50] #++
	#seq_this=store_dict[chr][int(start)-50:int(start)+500] #--,-+,+-
	print ("%s\t" % (tag),file=output,end="")
	for i in range(len(seq_this)):
		A=0
		T=0
		C=0
		G=0
		seq_cal=seq_this[i:i+5]
		A=seq_cal.count('A')
		T=seq_cal.count('T')
		C=seq_cal.count('C')
		G=seq_cal.count('G')
		AT=(A+T)/5
		if AT<=0.4:
			AT=0.4
		print ("%s\t" % (AT),file=output,end="")
	print ("",file=output)