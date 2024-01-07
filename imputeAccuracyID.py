#coding:utf-8
import sys
from math import sqrt
import gzip

#python2 ../imputeAccuracyID.py allbam.chr10.phased.vcf.gz allbam.chr10.pass.vcf.gz ../allbam.chr10.pass.vcf.gz ../id.txt bcf+bg4.txt

phased=gzip.open(sys.argv[1]) #vcf file with after imputation: SNP2
passsnp=gzip.open(sys.argv[2]) #vcf file with after pass:  SNP1
highdepth=gzip.open(sys.argv[3]) #vcf file with high depth:  SNP3
selectID=open(sys.argv[4])    #ID file used for validated imputation accuracy

q=open(sys.argv[5],'w')
q.write("ID"+'\t'+'Nsnps'+'\t'+'Concordance'+'\t'+'Correlation'+'\n')

def multiply(a,b):
    sum_ab=0.0
    for i in range(len(a)):
        temp=a[i]*b[i]
        sum_ab+=temp
    return sum_ab

def calc_corr(x, y):
    n=len(x)
    sum_x=sum(x)
    sum_y=sum(y)
    sum_xy=multiply(x,y)
    sum_x2 = sum([pow(i,2) for i in x])
    sum_y2 = sum([pow(j,2) for j in y])
    molecular=sum_xy-(float(sum_x)*float(sum_y)/n)
    denominator=sqrt((sum_x2-float(sum_x**2)/n)*(sum_y2-float(sum_y**2)/n))
    return molecular/denominator

def concordance(a,b):
	p=0.0
	for i in range(0,len(a)):
		if a[i]==b[i]:
			p+=1
	return p/len(a)
	
def common(a,b):
	commonEle=[val for val in a if val in b]
	return commonEle

def recode1(list):
	if list=='0|1' or list=='1|0':
		list='1'
	elif list=='0|0':
		list='0'	
	elif list=='1|1':
		list='2'	
	return list

def recode2(list):
	if list=='0/1':
		list='1'
	elif list=='0/0':
		list='0'	
	elif list=='1/1':
		list='2'
	elif list=='./.':
        	list='NaN'
	return list
	
#--提取SNP和id
line=passsnp.readlines()
ID1=[]
SNP1={}
tt=[]
for i in line[:]:
	if i.startswith(b'#CHROM'):
		f=i.split()
		for e in range(9,len(f)):		
			ID1.append(f[e])
	elif i.startswith(b'##'):
		pass
	else:
		f=i.split()
		for e in range(9,len(f)):	
			SNP1[f[1]]=tt
			tt.append(f[e].split(':')[0])
		tt=[]		
passsnp.close()

#--提取SNP和id
line=highdepth.readlines()
ID3=[]
SNP3={}
tt=[]
for i in line[:]:
	if i.startswith(b'#CHROM'):
		f=i.split()
		for e in range(9,len(f)):		
			ID3.append(f[e])
	elif i.startswith(b'##'):
		pass
	else:
		f=i.split()
		for e in range(9,len(f)):	
			SNP3[f[1]]=tt
			tt.append(f[e].split(':')[0])
		tt=[]

dick3={}
for i,j in enumerate(ID3):
	dick3[j]=i
# print(ID3)
		
highdepth.close()
# print(SNP3)

#--提取SNP和id
ID2=[]
ID={}
tt=[]
SNP2={}
for i in phased.readlines()[:]:
	if i.startswith(b'#CHROM'):
		f=i.split()
		for e in range(9,len(f)):		
			ID[f[e]]=[]
			ID2.append(f[e])
	elif i.startswith(b'##'):
		pass
	else:
		f=i.split()
		for e in range(0,len(ID2)):	
			if f[1] in SNP1:
				if f[e+9].split(':')[0]!='./.' and SNP1[f[1]][e]=='./.':
					ID[ID2[e]].append(f[1])
		for e in range(9,len(f)):	
			SNP2[f[1]]=tt
			tt.append(f[e].split(':')[0])
		tt=[]		
		
		
phased.close()

sid={}
for i in selectID:
	f=i.split()
	sid[f[0]]=f[0]
selectID.close()
	
geno1=[]
geno2=[]
for i,j in enumerate(ID2):
	if j in sid:
		# print(j,dick3[j])
		for e in ID[j]:
			if e in SNP3:
				if SNP3[e][dick3[j]]!='./.':
					if SNP2[e][i][1]=='|':  #for beagle
						geno1.append(float(recode1(SNP2[e][i])))
						geno2.append(float(recode2(SNP3[e][dick3[j]])))
					elif SNP2[e][i][1]=='/': #for stitch
						geno1.append(float(recode2(SNP2[e][i])))
						geno2.append(float(recode2(SNP3[e][dick3[j]])))					
		print(j)
		# print(geno1)
		# print(geno2)
		cor=calc_corr(geno1,geno2)
		acc=concordance(geno1,geno2)
		print (cor)
		print (acc)	
		q.write(j+'\t'+str(len(geno1))+'\t'+str(acc)+'\t'+str(cor)+'\n')
		geno1=[]
		geno2=[]
q.close()
