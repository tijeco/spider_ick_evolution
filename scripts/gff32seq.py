from Bio.Seq import Seq
import sys

def getIntron(geneList,seq,strand):

    for exon_num in range(len(geneList)):

        exon = geneList[exon_num]
        # print(1,exon,exon_num)

        if (len(geneList) - exon_num ) > 1:
            # print(exon)
            intron = seq[geneList[exon_num +0][1]+1+1:geneList[exon_num +1][0]]
            if strand == "-":
                intron = Seq(intron).reverse_complement()

                # print(intron[:5],intron[-5:])
            else:

                print(intron[:5],intron[-5:])


gene_dict = {}
seq_dict = {}

# with open("data/species/Stegodyphus_mimosarum/KK114158.fa") as f:
with open(sys.argv[1]) as f:
    header = ""
    for line in f:
        if line[0] == ">":
            header = line.strip().strip(">")
            seq_dict[header] = ""
        else:
            seq_dict[header] += line.strip()
with open(sys.argv[2]) as f:
# with open("data/species/Stegodyphus_mimosarum/KK114158.gff3") as f:
    contig = ""
    gene_id = ""
    for line in f:
        row = line .strip().split()
        if "#" in line:
            continue
        # if "##sequence-region" in line:
            # contig = row[1]
            # gene_dict[contig] = {}
            # gene_dict[contig]["length"] = int(row[3])

            # gene_dict[row[1]]["exons"] = []
            # gene_dict[row[1]]["introns"] = []
            # continue
        else:
            contig = row[0]
            if contig not in gene_dict:
                gene_dict[contig] = {}
        # if row[0] not in gene_dict:
            # contig = row[0]
        if row[2] == "mRNA":
            # contig = row[1]
            gene_id = row[8].split(";")[0].split(":")[1]
            gene_dict[contig][gene_id] = {"exons":[],"introns":[],"strand":row[6]}

        if row[2] == "CDS":
            # rf = int(row[7])
            # rf_adjust_left = int(-1.5*(rf**2) + 3.5*rf)
            # rf_adjust_right = int(-1.5*(rf**2) + 2.5*rf +1)
            # print(rf,rf_adjust_left)
            gene_dict[contig][gene_id]["exons"].append((int(row[3])-1,int(row[4])-1))
# print(gene_dict)

# print(gene_dict["KK114158"]["KFM61793"]["exons"])
# gene_dict["KK114158"].pop("length")
for contig in gene_dict:
    for gene in gene_dict[contig]:
        # print(gene)
        gene_posList = gene_dict[contig][gene]["exons"]
        # print(gene_posList,seq_dict[contig][0],gene_dict[contig][gene]["strand"])
        getIntron(gene_posList,seq_dict[contig],gene_dict[contig][gene]["strand"])


# KFM61792 = gene_dict["KK114158"]["KFM61792"]["exons"]
# KFM61793 = gene_dict["KK114158"]["KFM61793"]["exons"]
# print(KFM61793)
# KFM61794 = gene_dict["KK114158"]["KFM61794"]["exons"]
# KFM61795 = gene_dict["KK114158"]["KFM61795"]["exons"]

# print(gene_dict["KK114158"]["KFM61793"]["exons"])
# print(seq_dict["KK114158"][KFM61793[0][0]:KFM61793[-1][1]])
# gene = seq_dict["KK114158"][KFM61793[0][0]:KFM61793[-1][1]]

# print(KFM61793)
# ATG CAT TAA
# XTG CAT TAA x
# XXG CAT TAA xx
########## correct #########
# ATGAAGCTGCTGATTTTTGCTGTAGTTCTTTTCGTAGTTGCAACTTCTATCGCTTCTCCTGCGAAA
# GATGAATCCCCGATGGAGCCTCAAGAA
#    AT GAA TCCCCGATGGAGCCTCAAGAAA
#    -1,-1
# # AGAGGCTCATTTTGTACAGATTATAAGGAAAGCTGCATTTCTAAAACGGATTGCTGCAGGGGTCAGTGTCTGTGTGACTTTGCTGGTAATAACTGCGTCTGCGCC
#   XXA GGC TCATTTTGTACAGATTATAAGGAAAGCTGCATTTCTAAAACGGATTGCTGCAGGGGTCAGTGTCTGTGTGACTTTGCTGGTAATAACTGCGTCTGCGCX
#  -2,+1
# MKLLIFAVVLFVVATSIASPAK
# DESPMEPQE
# RGSFCTDYKESCISKTDCCRGQCLCDFAGNNCVCA


# MKLLIFAVVLFVVATSIASPAK
# ESPMEPQE
# GSFCTDYKESCISKTDCCRGQCLCDFAGNNCVC
###########################

# exon_1 = seq_dict["KK114158"][KFM61793[0][0]:KFM61793[0][1]+1] # exon_1
# print(exon_1)
# intron_1 = seq_dict["KK114158"][KFM61793[0][1]+1+1:KFM61793[1][0]]
# print(intron_1[:5],intron_1[-5:])
# exon_2 = seq_dict["KK114158"][KFM61793[1][0]:KFM61793[1][1]+1] # exon_2
# print(exon_2)
# intron_2 = seq_dict["KK114158"][KFM61793[1][1]+1+1:KFM61793[2][0]]
# print(intron_2[:5],intron_2[-5:])
# exon_3 = seq_dict["KK114158"][KFM61793[2][0]:KFM61793[2][1]+1] # exon_3
# print(exon_3)
# print()
# print(exon_1+exon_2+exon_3)
# print("\nKFM61794")
# exon_1 = seq_dict["KK114158"][KFM61794[0][0]:KFM61794[0][1]+1] # exon_1
# print(exon_1)
# intron_1 = seq_dict["KK114158"][KFM61794[0][1]+1+1:KFM61794[1][0]]
# print(intron_1[:5],intron_1[-5:])
# exon_2 = seq_dict["KK114158"][KFM61794[1][0]:KFM61794[1][1]+1] # exon_2
# print(exon_2)
# intron_2 = seq_dict["KK114158"][KFM61794[1][1]+1+1:KFM61794[2][0]]
# print(intron_2[:5],intron_2[-5:])
# exon_3 = seq_dict["KK114158"][KFM61794[2][0]:KFM61794[2][1]+1] # exon_3
# print(exon_3)
# print()
# print(exon_1+exon_2+exon_3)
#
# print("KFM61795")
# print(KFM61795)
# exon_1 = Seq(seq_dict["KK114158"][KFM61795[0][0]:KFM61795[0][1]+1]).reverse_complement()
# print(exon_1)
# intron_1 = Seq(seq_dict["KK114158"][KFM61795[0][1]+1+1:KFM61795[1][0]])
# print(intron_1[:5],intron_1[-5:])
# exon_2 = Seq(seq_dict["KK114158"][KFM61795[1][0]:KFM61795[1][1]+1]).reverse_complement()
# print(exon_2)
# intron_2 = Seq(seq_dict["KK114158"][KFM61795[1][1]+1+1:KFM61795[2][0]])
# print(intron_2[:5],intron_2[-5:])
# exon_3 = Seq(seq_dict["KK114158"][KFM61795[2][0]:KFM61795[2][1]+1]).reverse_complement()
# print(exon_3)
# print()
# print(exon_3+exon_2+exon_1)
#
# print("KFM61792")
# print(KFM61792)
# exon_1 = Seq(seq_dict["KK114158"][KFM61792[0][0]:KFM61792[0][1]+1]).reverse_complement()
# print(exon_1)
# intron_1 = Seq(seq_dict["KK114158"][KFM61792[0][1]+1+1:KFM61792[1][0]])
# print(1,intron_1[:5],intron_1[-5:])
# exon_2 = Seq(seq_dict["KK114158"][KFM61792[1][0]:KFM61792[1][1]+1]).reverse_complement()


# print(exon_2)
# intron_2 = Seq(seq_dict["KK114158"][KFM61792[1][1]+1+1:KFM61792[2][0]])
# print(intron_2[:5],intron_2[-5:])
# exon_3 = Seq(seq_dict["KK114158"][KFM61792[2][0]:KFM61792[2][1]+1]).reverse_complement()
# print(exon_3)

# print()
# print(exon_2+exon_1)

# print(KFM61794)
# print(seq_dict["KK114158"][KFM61794[0][0]:KFM61794[0][1]]) # exon_1
# print(seq_dict["KK114158"][KFM61794[0][1]:KFM61794[0][1]+5])
# print(seq_dict["KK114158"][KFM61794[1][0]-6:KFM61794[1][0]-1])
# print(seq_dict["KK114158"][KFM61794[1][0]-1:KFM61794[1][1]-1]) # exon_2
# print(seq_dict["KK114158"][KFM61794[1][1]-1:KFM61794[1][1]+6])
# print(seq_dict["KK114158"][KFM61794[2][0]-7:KFM61794[2][0]-2])
# print(seq_dict["KK114158"][KFM61794[2][0]-2:KFM61794[2][1]+1]) # exon_3

# print(KFM61794)
# print(seq_dict["KK114158"][KFM61794[0][0]:KFM61794[0][1]])
# print(seq_dict["KK114158"][KFM61794[1][0]-1:KFM61794[1][1]-1])
# print(seq_dict["KK114158"][KFM61794[2][0]-2:KFM61794[2][1]+1])
# MKFLSPLVFLALIVFVISNVMASPYS
# EESLKGNNAAE
# RALKGDCIRKEKDCLPDPVDDKCCEGLQCKCTYFGGYCVCMEKTIKTR
#
# MKFLSPLVFLALIVFVISNVMASPYSEESLKGNNAAERALKGDCIRKEKDCLPDPVDDKCCEGLQCKCTYFGGYCVCMEKTIKTR
# print(seq_dict["KK114158"][47246-1:47312-1])
# print(seq_dict["KK114158"][47312-1:47312+5])
# exon_1 = Seq(seq_dict["KK114158"][KFM61793[0][0]:KFM61793[0][1]])
# exon_2 = Seq(seq_dict["KK114158"][KFM61793[1][0]:KFM61793[1][1]])
# exon_3 = Seq(seq_dict["KK114158"][KFM61793[2][0]:KFM61793[2][1]])
#
# intron_1 = Seq(seq_dict["KK114158"][KFM61793[0][1]:KFM61793[1][0]])
# intron_2 = Seq(seq_dict["KK114158"][KFM61793[1][1]:KFM61793[2][0]])
# print(len(gene))
# print(len(exon_1)+len(exon_2)+len(exon_3)+len(intron_1)+len(intron_2))
# print(exon_1,len(exon_1)%3,"\n",exon_2,len(exon_2)%3,"\n",exon_3,len(exon_3)%3,"\n")

# print(exon_1,len(exon_1),len(exon_1)%3, exon_1.translate())
# print(exon_2,len(exon_2),len(exon_2)%3, exon_2.translate()) #GAATCCCCGATGGAGCCTCAAGAAAG
# print(exon_3,len(exon_3),len(exon_3)%3, exon_3.translate())
# print(exon_1.translate()+exon_2.translate()+exon_3.translate())
# print(intron_1[:5]+"...N"+str(len(intron_1)-10)+"..."+intron_1[-5:])
# print(intron_2[:5]+"...N"+str(len(intron_2)-10)+"..."+intron_2[-5:])
# print(exon_1,intron_1[:5]+"...N"+str(len(intron_1)-10)+"..."+intron_1[-5:], exon_2,intron_2[:5]+"...N"+str(len(intron_2)-10)+"..."+intron_2[-5:],exon_3)



# print(gene)
# print(seq_dict["KK114158"][56645:56672])
# print(seq_dict["KK114158"][56644:56672])

# print(seq_dict["KK114158"][KFM61793[2][0]:KFM61793[2][1]])


# KFM61792 = gene_dict["KK114158"]["KFM61792"]["exons"]
# exon_1 = Seq(seq_dict["KK114158"][KFM61792[0][0]:KFM61792[0][1]]).reverse_complement()
# exon_2 = Seq(seq_dict["KK114158"][KFM61792[1][0]:KFM61792[1][1]]).reverse_complement()
# # exon_3 = Seq(seq_dict["KK114158"][KFM61792[2][0]:KFM61792[2][1]]).reverse_complement()
#
# intron_1 = Seq(seq_dict["KK114158"][KFM61792[0][1]:KFM61792[1][0]]).reverse_complement()
# # intron_2 = Seq(seq_dict["KK114158"][KFM61792[1][1]:KFM61792[2][0]]).reverse_complement()
# print(len(gene))
# print(len(exon_1)+len(exon_2)+len(exon_3)+len(intron_1)+len(intron_2))
# print(exon_1, exon_1.translate())
# print(exon_2, exon_2.translate()) #GAATCCCCGATGGAGCCTCAAGAAAG
# print(exon_3, exon_3.translate())
# print(exon_2.translate()+exon_1.translate())#+exon_3.translate())
# print(intron_1[:5]+"...N"+str(len(intron_1)-10)+"..."+intron_1[-5:])
# # print(intron_2[:5]+"...N"+str(len(intron_2)-10)+"..."+intron_2[-5:])
# print(exon_2,intron_1[:5]+"...N"+str(len(intron_1)-10)+"..."+intron_1[-5:], exon_1)#,intron_2[:5]+"...N"+str(len(intron_2)-10)+"..."+intron_2[-5:],exon_3)

# print(seq_dict["KK114158"][KFM61792[0][0]:KFM61792[-1][1]])
# print(Seq(seq_dict["KK114158"][KFM61792[0][0]:KFM61792[0][1]]).reverse_complement())
# print(Seq(seq_dict["KK114158"][KFM61792[1][0]:KFM61792[1][1]]).reverse_complement())
# print(gene_dict["KK114158"]["KFM61792"]["exons"])
# print(gene_dict)
# ATGAAGCTGCTGATTTTTGCTGTAGTTCTTTTCGTAGTTGCAACTTCTATCGCTTCTCCTGCGAAAG
# print(KFM61793[0])
# MKLLIFAVVLFVVATSIASPAK
# ESPMEPQE
# GSFCTDYKESCISKTDCCRGQCLCDFAGNNCVCA
# MKLLIFAVVLFVVATSIASPAKESPMEPQEGSFCTDYKESCISKTDCCRGQCLCDFAGNNCVCA
#
# MKLLIFAVVLFVVATSIASPAK
# DESPMEPQE
# RGSFCTDYKESCISKTDCCRGQCLCDFAGNNCVCA
#
# MKLLIFAVVLFVVATSIASPAKDESPMEPQERGSFCTDYKESCISKTDCCRGQCLCDFAGNNCVCA
