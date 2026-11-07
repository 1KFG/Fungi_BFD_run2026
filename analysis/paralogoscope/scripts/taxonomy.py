#!/usr/bin/python3.12
"""Taxonomy table for the genome keys of the paralogoscope run: samples.csv columns joined by makeSampleTag(SPECIES, STRAIN) and the representative flag.
Usage: taxonomy.py samples_all.csv repr_assignments.tsv out.tsv"""
import sys, re, csv
import pandas as pd
def clean(st):
    st=(st or "").strip(); st=re.sub(r"['\"]","",st).split(";")[0].strip().replace(":"," ")
    st=re.sub(r"^\s*\*+","",st); st=re.sub(r"\*+\s*$","",st); st=re.sub(r"\s*\*+\s*","-",st); st=re.sub(r"[&`$|<>()]+"," ",st).strip(); return st
def tag(sp,st):
    sp=re.sub(r"['\"]","",(sp or "").strip()); j="_".join(x for x in [sp,clean(st)] if x); return re.sub(r"[\s/#\[\]?{}]+","_",j)
s=pd.read_csv(sys.argv[1],dtype=str,keep_default_na=False)
s["genome"]=[tag(a,b) for a,b in zip(s.SPECIES,s.STRAIN)]
r=pd.read_csv(sys.argv[2],sep="\t"); rep=set(r[r.is_representative].out)
s["is_representative"]=s.genome.isin(rep).astype(int)
cols=["genome","SPECIES","STRAIN","PHYLUM","SUBPHYLUM","CLASS","SUBCLASS","ORDER","FAMILY","GENUS","TRANSL_TABLE","is_representative","NCBI_TAXONID","ASMID"]
t=s[cols].copy()
for c in ["PHYLUM","SUBPHYLUM","CLASS","SUBCLASS","ORDER","FAMILY","GENUS"]: t[c]=t[c].replace("", pd.NA)
t.to_csv(sys.argv[3],sep="\t",index=False)
print(len(t),"genomes;",{c:int(t[c].notna().sum()) for c in ["PHYLUM","CLASS","ORDER","FAMILY","GENUS"]},"| representatives",int(t.is_representative.sum()))
