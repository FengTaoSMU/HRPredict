#!/usr/bin/python3
#-----------------------------------------------
# Software: HRPredict
# Author: FengTaoSMU
# Function: Reference distance matrix generation
#-----------------------------------------------

# 'prokka /data3/Group7/fengtao/2.MOBFinder/1.plasmid/01.complete_fasta/CP100462.1.fasta --outdir ./test --prefix NZ_CP050157.1 --kingdom Bacteria'

from Bio import SeqIO
import argparse
import re, os, copy
import numpy as np
import pandas as pd

def parameter_passing():
    usage = 'python3 %(prog)s [-h] [-i input_fasta] [-m model_dir] [-r reference_dir] [-o output_dir]'
    main = 'Generating reference distance matrix for input plasmid genomes'
    parser = argparse.ArgumentParser(usage = usage, description = main)
    parser.add_argument('-i', type = str, metavar = '', default = '', help = "input plasmid fasta file")
    parser.add_argument('-m', type = str, metavar = '', default = '', help = "model directory")
    parser.add_argument('-r', type = str, metavar = '', default = '', help = "reference directory")
    parser.add_argument('-o', type = str, metavar = '', default = '', help = "output directory")
    parser.add_argument('-v', action="version", version = "%(prog)s version 1.0", help = "show script's version and exit")
    parameter = parser.parse_args()
    return parameter.i, parameter.m, parameter.r, parameter.o

def dir_check(input_dir):
    if not os.path.exists(input_dir):
        os.mkdir(input_dir)
    if not input_dir.endswith('/'):
        input_dir = input_dir + '/'
    return input_dir

def accession_id_get(fasta, fa_dir):
    num, plist = 0, []
    with open(fasta) as fa:
        for line in fa:
            if '>' in line:
                if num != 0:
                    output.write(">%s\n%s\n" % (pid, seq))
                    output.close()
                pid = line.strip().split()[0][1:]
                plist.append(pid)
                output = open("{}{}.fasta".format(fa_dir, pid), 'w')
                seq = ''
                num += 1
            else:
                seq += line.strip()
        output.write(">%s\n%s\n" % (pid, seq))
        output.close()
    return plist

def word_vector_get(word_vector):
    word_dic = {}
    with open(word_vector) as word_vec:
        for line in word_vec:
            word = line.strip().split()
            word_dic[word[0]] = [float(num) for num in word[1:]]
    return word_dic

def fasta_seq_get(fasta):
    seq = ''
    with open(fasta) as fa:
        for line in fa:
            if '>' not in line:
                seq += line.strip()
    return seq

def chain_word_embedding(word_dic, chain, word_len):
    kmer_num = 0
    kmer_sum = np.array([0]*100)
    for i in range(0, len(chain)- word_len + 1):
        kmer = chain[i: i + word_len]
        if kmer in word_dic:
            kmer_arr = np.array(word_dic[kmer])
            kmer_sum = kmer_sum + kmer_arr
            kmer_num += 1
        else:
            continue
    if kmer_num != 0:
        kmer_sum = kmer_sum / kmer_num
    return kmer_sum

def ref_prot_dic_get(prot, prot_vec):
    dic, num = {}, 0
    with open(prot) as ref_prot:
        for line in ref_prot:
            if '>' in line:
                num += 1
                tag = 'prot' + str(num)
            else:
                seq = line.strip()
                vec = chain_word_embedding(prot_vec, seq, 3)
                dic[tag] = vec
    return dic

def marker_word_embedding(seq, dna_word, kmer_len):
    seq_com = dna_complement_chain(seq)
    seq_word = chain_word_embedding(dna_word, seq, kmer_len)
    seq_word_com = chain_word_embedding(dna_word, seq_com, kmer_len)
    seq_vec = (seq_word + seq_word_com) / 2
    return seq_vec

def genebank_protein_dic_get(genebank, prot_vector):
    dic = {}
    num = 0
    for record in SeqIO.parse(genebank, 'gb'):
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                num += 1
                cds_seq = feature.qualifiers.get('translation', [''])[0]
                cds_vec = chain_word_embedding(prot_vector, cds_seq, 3)
                cds_tag = 'prot_' + str(num)
                dic[cds_tag] = cds_vec
    if num == 0:
        dic['prot_1'] = np.array([float(0)]*100)
    return dic

def prot_dist_calculate(plasmid_id, fasta, ref_prot_dic, prot_vec, tmp_dir):
    prokka_tmp = "{}{}".format(tmp_dir, plasmid_id)
    prokka_run = "prokka {} --outdir {} --prefix {} --kingdom Bacteria".format(fasta, prokka_tmp, plasmid_id)
    prokka_gb = "{}/{}.gbk".format(prokka_tmp, plasmid_id)
    os.system(prokka_run)
    if not os.path.getsize(prokka_gb):
        dic = {}
        dic['prot_1'] = np.array([float(0)]*100)
    else:
        dic = genebank_protein_dic_get(prokka_gb, prot_vec)
    result = np.array([float(0)]*len(ref_prot_dic))
    num = 0
    for prot in ref_prot_dic:
        dist_list = [np.linalg.norm(vec - ref_prot_dic[prot]) for vec in dic.values()]
        result[num] = min(dist_list)
        num += 1
    os.system("rm -rf {}".format(prokka_tmp))
    return result

def plasmid_embedding_generate(acc_list, ref_dir, prot_vector, fa_dir, output_dir):
    output = open("{}protein_feature_matrix.tsv".format(output_dir), 'w')
    tmp_dir = "{}tmp".format(output_dir)
    tmp_dir = dir_check(tmp_dir)
    ref_prot = "{}reference_protein.fasta".format(ref_dir)
    ref_prot_dic = ref_prot_dic_get(ref_prot, prot_vector)
    for acc in acc_list:
        fasta = "{}{}.fasta".format(fa_dir, acc)
        dist_total = prot_dist_calculate(acc, fasta, ref_prot_dic, prot_vector, tmp_dir)
        output.write("{}".format(acc))
        for edist in dist_total:
            output.write("\t{}".format(edist))
        output.write("\n")
        os.system("rm -rf {}".format(fasta))
    output.close()

def lineage_get(input_lineage):
    lineage = []
    with open(input_lineage) as lineage_input:
        for line in lineage_input:
            lineage.append(line.strip().split(','))
    return lineage

def host_range_output(result_dic, output, prob_F, prob_G, prob_S, lineageF, lineageG, lineageS, result_level):
    FP_df = pd.read_csv(prob_F, delimiter=',')
    GP_df = pd.read_csv(prob_G, delimiter=',')
    SP_df = pd.read_csv(prob_S, delimiter=',')

    g_dic = {}
    for h_index in range(len(lineageG)):
        g_dic[lineageG[h_index][1]] = lineageG[h_index][0]

    s_dic = {}
    for h_index in range(len(lineageS)):
        s_dic[lineageS[h_index][2]] = [lineageS[h_index][0], lineageS[h_index][1]]
        g_dic[lineageS[h_index][1]] = lineageS[h_index][0]

    if result_level == 'F':
        for pid in result_dic:
            output.write("{}".format(pid))
            if result_dic[pid]:
                hr_list = sorted(result_dic[pid])
                for host in hr_list:
                    output.write("\tf__{}({})".format(host, round(FP_df.loc[pid, host], 2)))
            else:
                output.write("\tNo family was predicted")
            output.write("\n")

    elif result_level == 'G':
        for pid in result_dic:
            output.write("{}".format(pid))
            if result_dic[pid]:
                hr_list = sorted(result_dic[pid])
                for host in hr_list:
                    output.write("\tf__{}({})".format(g_dic[host], round(FP_df.loc[pid, g_dic[host]], 2)))
                    output.write(";g__{}({})".format(host, round(GP_df.loc[pid, host], 2)))
            else:
                output.write("\tNo genus was predicted")
            output.write("\n")

    else:
        for pid in result_dic:
            output.write("{}".format(pid))
            if result_dic[pid]:
                hr_list = sorted(result_dic[pid])
                for host in hr_list:
                    output.write("\tf__{}({})".format(s_dic[host][0], round(FP_df.loc[pid, s_dic[host][0]], 2)))
                    output.write(";g__{}({})".format(s_dic[host][1], round(GP_df.loc[pid, s_dic[host][1]], 2)))
                    output.write(";s__{}({})".format(host, round(SP_df.loc[pid, host], 2)))
            else:
                output.write("\tNo species was predicted")
            output.write("\n")

def host_range_predict(model_dir, ref_dir, tmp_dir, output_dir):
    input_d = "{}{}".format(tmp_dir, "protein_feature_matrix.tsv")
    rscirpt = "{}{}".format(model_dir, "HRPredict.R")

    cutoff_F = "{}{}".format(model_dir, "cutoff_Family.tsv")
    cutoff_G = "{}{}".format(model_dir, "cutoff_Genus.tsv")
    cutoff_S = "{}{}".format(model_dir, "cutoff_Species.tsv")

    os.system("Rscript {} {} {} {} {} {} {}".format(rscirpt, model_dir, cutoff_F, cutoff_G, cutoff_S, input_d, tmp_dir))

    lineageF = lineage_get("{}{}".format(model_dir, "lineage_Family.tsv"))
    lineageG = lineage_get("{}{}".format(model_dir, "lineage_Genus.tsv"))
    lineageS = lineage_get("{}{}".format(model_dir, "lineage_Species.tsv"))

    prob_F = "{}{}".format(tmp_dir, "prob_result_family.tsv")
    prob_G = "{}{}".format(tmp_dir, "prob_result_genus.tsv")
    prob_S = "{}{}".format(tmp_dir, "prob_result_species.tsv")

    result_F = "{}{}".format(tmp_dir, "result_family.tsv")
    result_G = "{}{}".format(tmp_dir, "result_genus.tsv")
    result_S = "{}{}".format(tmp_dir, "result_species.tsv")

    result_F_dic = {}
    result_G_dic = {}
    result_S_dic = {}

    with open(result_F) as Family_result:
        num, header = 0, []
        for line in Family_result:
            line_id = ''
            line = line.strip().split(',')
            if num == 0:
                header = line
                num += 1
            else:
                for item in range(len(line)):
                    if item == 0:
                        line_id = line[item]
                        result_F_dic[line_id] = []
                    else:
                        if line[item] == "Yes":
                            result_F_dic[line_id].append(header[item])

    with open(result_G) as Genus_result:
        num, header = 0, []
        for line in Genus_result:
            line_id = ''
            line = line.strip().split(',')
            if num == 0:
                header = line
                num += 1
            else:
                for item in range(len(line)):
                    if item == 0:
                        line_id = line[item]
                        result_G_dic[line_id] = []
                    else:
                        if line[item] == "Yes":
                            genus_list_1, genus_list_2 = [], []
                            genus_list_1 = [glist for glist in lineageG if header[item] in glist]
                            genus_list_2 = [glist for glist in lineageS if header[item] in glist]
                            if genus_list_1:
                                if genus_list_1[0][0] in result_F_dic[line_id]:
                                    result_G_dic[line_id].append(header[item])
                            if genus_list_2:
                                if genus_list_2[0][0] in result_F_dic[line_id]:
                                    result_G_dic[line_id].append(header[item])

    with open(result_S) as Species_result:
        num, header = 0, []
        for line in Species_result:
            line_id = ''
            line = line.strip().split(',')
            if num == 0:
                header = line
                num += 1
            else:
                for item in range(len(line)):
                    if item == 0:
                        line_id = line[item]
                        result_S_dic[line_id] = []
                    else:
                        if line[item] == "Yes":
                            species_list = []
                            species_list = [slist for slist in lineageS if header[item] in slist]
                            if species_list:
                                if species_list[0][1] in result_G_dic[line_id]:
                                    result_S_dic[line_id].append(header[item])

    output_F = open("{}{}".format(output_dir, "HostRange_Family.tsv"), 'w')
    output_G = open("{}{}".format(output_dir, "HostRange_Genus.tsv"), 'w')
    output_S = open("{}{}".format(output_dir, "HostRange_Species.tsv"), 'w')

    host_range_output(result_F_dic, output_F, prob_F, prob_G, prob_S, lineageF, lineageG, lineageS, 'F')
    host_range_output(result_G_dic, output_G, prob_F, prob_G, prob_S, lineageF, lineageG, lineageS, 'G')
    host_range_output(result_S_dic, output_S, prob_F, prob_G, prob_S, lineageF, lineageG, lineageS, 'S')

    output_F.close()
    output_G.close()
    output_S.close()

def plasmid_host_range_predict(fasta, model_dir, ref_dir, output_dir):
    tmp_dir = "{}{}/".format(output_dir, "tmp")
    tmp_dir = dir_check(tmp_dir)
    prot_vector = "{}{}".format(ref_dir, "protvec.txt")
    prot_vec = word_vector_get(prot_vector)
    fa_dir = "{}fasta".format(tmp_dir)
    fa_dir = dir_check(fa_dir)
    acc_list = accession_id_get(fasta, fa_dir)
    plasmid_embedding_generate(acc_list, ref_dir, prot_vec, fa_dir, tmp_dir)
    host_range_predict(model_dir, ref_dir, tmp_dir, output_dir)
    os.system("rm -rf {}".format(tmp_dir))

def main():
    input_fasta, model_dir, ref_dir, output_dir = parameter_passing()
    model_dir = dir_check(model_dir)
    ref_dir = dir_check(ref_dir)
    output_dir = dir_check(output_dir)
    plasmid_host_range_predict(input_fasta, model_dir, ref_dir, output_dir)

if __name__ == '__main__':
    main()

