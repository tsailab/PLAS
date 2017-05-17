#!/usr/bin/python
# -*- coding: utf-8 -*-

from sys import argv
import sys
import re



class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


def load_annotation():
    hash_desc = AutoVivification()
    hash_evalue = AutoVivification()
    #get_gene= re.compile("(.+)\.\d+$")
    file_list = ["ath","rice","ptri"]
    for file_in in file_list :
        print("Read "+file_in)
        file_read = "Peuph_"+file_in+"_desc_out.tab"   #Shanmu_ptri_desc_out.tab
        with open(file_read,'r') as FIN:
            for line in FIN:
                record = line.rstrip().split("\t")
                #ggg = get_gene.search(record[0])
                #gene=ggg.group(1)
                gene = record[0]
                hit = record[5]
                evalue = float(record[10])
                if file_in == 'ptri' :
                    desc = ''
                else :
                    desc = record[12]
                gene_out=gene+file_in
                if gene_out in hash_evalue.keys():
                    if evalue < hash_evalue[gene_out]:
                        hash_evalue[gene_out] = evalue
                        hash_desc[gene_out] = hit+"\t"+desc
                else :
                    hash_evalue[gene_out] = evalue
                    hash_desc[gene_out] = hit + "\t" + desc
    return hash_desc


def read_gene(file_in, file_out, hash_in):
    with open(file_in, 'r') as FIN:
        file_list = ["ath", "rice",  "ptri"]
        FOUT = open(file_out,'w')
        FOUT.write("gene"+"\t"+"\t".join(["ath","ath_desc","rice","rice_desc","ptri","ptri_desc"])+"\n")
        tracking = 0
        gene2trans = AutoVivification()
        for line in FIN:
            tracking+=1
            if int(tracking /1000)==tracking /1000:
                print("Read "+str(tracking)+" lines")
            record = line.rstrip().split("\t")
            gene = record[0]
            if len(record) >1 :
                trans = record[1]
                gene2trans[gene][trans] = 1
        # end of for

        for gene in gene2trans.keys():
            out_line = [gene]
            for file_in in file_list:
                #print ("Check "+file_in)
                find_str = 0
                for trans in gene2trans[gene].keys():
                    gene_out = trans + file_in
                    if gene_out in hash_in.keys():
                        if find_str == 0:
                            out_line.append(hash_in[gene_out])
                            find_str = 1
                if  find_str == 0:
                    out_line.append('NA' + "\t" + 'NA')
            FOUT.write("\t".join(out_line) + "\n")
        FOUT.close()


 # main

if __name__ == "__main__":
    script, file_out = argv
    anno_dict = load_annotation()
    print ("Load Annotation done")
    read_gene("GCF_000495115.1_gene_list.tab", file_out , anno_dict)

## Author : lxue@uga.edu