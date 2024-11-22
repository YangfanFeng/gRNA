import re
import sys
from PyQt5 import QtWidgets, uic, QtCore
from fuzzysearch import find_near_matches

qtCreatorFile = "gRNA.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
#_translate=QtCore.QCoreApplication.translate
class MainUi(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.sel_goalseq_button.clicked.connect(self.SelGoal)
        self.sel_simseq_button.clicked.connect(self.SelSimseq)
        self.sel_gfffile_button.clicked.connect(self.Selanno)
        self.run_button.clicked.connect(self.Run)
        self.fileDialog = QtWidgets.QFileDialog(self)

    def SelGoal(self):
        self.goalfile_input,_filter= self.fileDialog.getOpenFileName()
        self.goalseq_input.setText("".join(self.goalfile_input) + "\n")
        print("input"+"".join(self.goalfile_input)+"\n")

    def SelSimseq(self):
        self.fasta_input, _filter = self.fileDialog.getOpenFileName()
        self.simseq_input.setText("".join(self.fasta_input) + "\n")
        print("input"+self.fasta_input+"\n")

    def Selanno(self):
        self.gfffile_input, _filter= self.fileDialog.getOpenFileName()
        self.gff_input.setText("".join(self.gfffile_input) + "\n")
        print("input"+self.gfffile_input+"\n")

    def Run(self):
        self.file_output,_filter=self.fileDialog.getSaveFileName()
        out = "Seq Searching Job Start"+"\n"
        self.output.setText(out)
        self.output.repaint()
        AR=AnnoResult
        gRNA_list=AR.read_file(self.goalfile_input)
        seq_gene=AR.readFasta(self.fasta_input)
        seq_list=[]
        for i in seq_gene:
            seq_list.append(i)
        out="Seq Scoring Job Start"+"\n"
        self.output.setText(out)
        self.output.repaint()
        match_list=AR.matchesRNA(seq_list,gRNA_list)
        site_list=AR.seq2mismatch(match_list)
        site_score=AR.sites2score(site_list)
        score_output=AR.end_score(site_score,gRNA_list)
        out += "Annotate Job Start"+"\n"
        self.output.setText(out)
        self.output.repaint()
        dict=AR.readgff(self.gfffile_input)
        anno_output=AR.exon_pos(dict,gRNA_list)
        AR.write_file(self.file_output,score_output,anno_output)
        out += "All Job Done""\n"
        self.output.setText(out)
        self.output.repaint()

class AnnoResult:
    def read_file(file_in):
        with open(file_in, "r") as f:
            li = f.read()
            line=li.strip().split()
            seq=line[1]
            out_put = []
            gg_site_list = []
            cc_site_list=[]
            gg_starts = []
            cc_starts = []
            output_lines = []
            for g in re.compile("(?=GG)").finditer(seq):
                gg_start = g.start()
                gg_starts.append(gg_start)  # 将gg开始位点构成列表
            for gg_start in gg_starts:
                if gg_start < 21 :
                    pass
                else:
                    output_line = seq[(gg_start - 21):(gg_start + 2)]  # 切片
                    if "N" in output_line:
                        pass
                    else:
                        gg_site_list.append([gg_start - 21,output_line])
            for c in re.compile("(?=CC)").finditer(seq):
                cc_start = c.start()
                cc_starts.append(cc_start)
            for cc_start in cc_starts:
                if len(seq) - cc_start < 21:
                    pass
                else:
                    output_line = seq[(cc_start):(cc_start + 23)]
                    if "N"in output_line:
                        pass
                    else:
                        cc_site_list.append([cc_start,output_line])

        return (gg_site_list,cc_site_list)
    def readFasta(f):
        with open(f, 'r') as FA:
            seqName, seq = '', ''
            while 1:
                line = FA.readline()
                line = line.strip('\n')
                if (line.startswith('>') or not line) and seqName:
                    yield ((seqName, seq))
                if line.startswith('>'):
                    seqName = line[1:]
                    seq = ''
                else:
                    seq += line
                if not line:
                    break

    def matchesRNA(sq,list):
        def complement(letters1):
            letters2=[]
            basecomplemt = {"A": "T", "T": "A", "G": "C", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g"}
            for base in letters1:
                if base in basecomplemt:
                    base = basecomplemt[base]
                else:
                    base=base
                letters2.append(base)
            return "".join(letters2)

        def revcomp(seq):
            return complement(seq)[::-1]
        gg_fuzzysearch=[]
        cc_fuzzysearch=[]
        gg_list=list[0]
        cc_list=list[1]
        for seq in gg_list:
            patten_seq=seq[1]
            find_string=[]
            ant_patten_seq=revcomp(patten_seq)
            for item in sq:
                dna_seq=item[1]
                matches = find_near_matches(patten_seq,dna_seq, max_substitutions=4, max_insertions=0, max_deletions=0)
                for c in matches:
                    find_string.append(dna_seq[c.start:c.end])
                ant_matches=find_near_matches(ant_patten_seq,dna_seq,max_substitutions=4,max_insertions=0,max_deletions=0)
                for c in ant_matches:
                    ant_find_seq=dna_seq[c.start:c.end]
                    find_seq=revcomp(ant_find_seq)
                    find_string.append(find_seq)
            if find_string==[]:
                find_string="green"
            gg_fuzzysearch.append([seq,find_string])
        for seq in cc_list:
            find_string=[]
            patten_seq=seq[1]
            ant_patten_seq=revcomp(patten_seq)
            for item in sq:
                dna_seq=item[1]
                matches=find_near_matches(ant_patten_seq,dna_seq,max_substitutions=4,max_insertions=0,max_deletions=0)
                for c in matches:
                    find_string.append(dna_seq[c.start:c.end])
                ant_matches=find_near_matches(patten_seq,dna_seq,max_substitutions=4,max_insertions=0,max_deletions=0)
                for c in ant_matches:
                    ant_find_seq=dna_seq[c.start:c.end ]
                    find_seq=revcomp(ant_find_seq)
                    find_string.append(find_seq)
            if find_string==[]:
                find_string="green"
            cc_fuzzysearch.append([seq,[ant_patten_seq,find_string]])
        return(gg_fuzzysearch,cc_fuzzysearch)
    def seq2mismatch(list):
        seq_out=[]
        gg_list=list[0]
        cc_list=list[1]
        for seq in gg_list:
            site_out=[]
            site_query_dna=seq[0]
            query_dna=seq[0][1]
            ref_dnas=seq[1]
            if ref_dnas=="green":
                seq_out.append(seq)
            else:
                for sim_seq in ref_dnas:
                    site=[]
                    ref_dna=sim_seq
                    if len(query_dna) != len(ref_dna):
                        pass
                    else:
                        i = 0
                        while i < len(query_dna):
                            if query_dna[i] != ref_dna[i]:
                                if i < 20:
                                    site.append(i)
                            i = i + 1

                    site_out.append(site)
                seq_out.append([site_query_dna,site_out])
        for seq in cc_list:
            site_out=[]
            site_query_dna=seq[0]
            score_seq=seq[1]
            query_dna=score_seq[0]
            ref_dnas=score_seq[1]
            if ref_dnas=="green":
                seq_out.append([seq[0],ref_dnas])
            else:
                for sim_seq in ref_dnas:
                    site=[]
                    ref_dna=sim_seq
                    if len(query_dna) != len(ref_dna):
                        pass
                    else:
                        i = 0
                        while i < len(query_dna):
                            if query_dna[i] != ref_dna[i]:
                                if i < 20:
                                    site.append(i)
                            i = i + 1
                    site_out.append(site)
                seq_out.append([site_query_dna,site_out])
        return seq_out

    def sites2score(mis_sites):
        seq_score=[]
        for item in mis_sites:
            site_name=item[0]
            all_mis_site=item[1]
            if all_mis_site=="green":
                seq_score.append(item)
            else:
                score=[]
                for site in all_mis_site:
                    mis_site=site
                    weight_all = (0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804,0.685, 0.583)
                    mis_site_num = len(mis_site)
                    score_p3 = 1
                    if mis_site_num >= 1:
                        score_p3 = 1 / (mis_site_num * mis_site_num)
                    score_p2= 1
                    if mis_site_num >= 2:
                        distance = []
                        i = 0
                        while i < mis_site_num - 1:
                            dis = mis_site[i + 1] - mis_site[i]
                            distance.append(dis)
                            i = i + 1
                        sum_dis = 0
                        a = 0
                        while a < len(distance):
                            sum_dis = sum_dis + distance[a]
                            a = a + 1
                        mean_dis = sum_dis / len(distance)
                        score_p2 = 1 / ((((19 - mean_dis) / 19) * 4) + 1)
                    score_p1 = 1
                    b = 0
                    while b < mis_site_num:
                        score_p1 = score_p1 * (1 - weight_all[mis_site[b]])
                        b = b + 1
                    score_hit = score_p1 * score_p2 * score_p3 * 100
                    score.append(score_hit)
                seq_score.append([site_name,score])
        return seq_score

    def end_score(dif_list, start_list):
        sum_score=0
        score_list=[]
        for c in dif_list:
            gRNA=c[0]
            scores=c[1]
            if scores=="green":
                score_list.append(c)
            else:
                for item in scores:
                    sum_score+=item
                grna_score=100/(100+sum_score)
                if grna_score<0.5:
                    grade="red"
                else:
                    grade="green"
                score_list.append([gRNA,grade])
        return score_list


    def readgff(fg):
        with open(fg,"r")as gffin:
            fg = gffin.readlines()
            gene1gff = {}
            d_keys = []
            exons = []
            for line in fg:
                if line.startswith("#"):
                    continue
                lin = line.rstrip().split("\t")
                type = lin[2]
                attribute = lin[8]
                if type == "gene":
                    hit = re.search(r"ID\=([^;]+)", attribute)
                    gene = hit.group(1)
                    d_keys.append(gene)
                if type == "exon":
                    exons.append(lin)
            for i in range(len(d_keys)):
                valuei = []
                for item in exons:
                    if d_keys[i] in item[8]:
                        valuei.append(item)
                gene1gff[d_keys[i]] = valuei
            return gene1gff
    def exon_pos(d,list):
        seq_len=int(23)
        intron_exon = []
        in_exon = []
        exon_intron = []
        out_exon = []
        site_list1=[]
        output_list=[]
        gg_list=list[0]
        cc_list=list[1]
        for item in gg_list:
            site_list1.append(item)
        for item in cc_list:
            site_list1.append(item)
        site_list1=sorted(site_list1,key=lambda d:d[0])
        for key in d.keys():
            exon_starts = []
            exon_ends = []
            for value in d[key]:
                exon_starts.append(int(value[3]))
                exon_ends.append(int(value[4]))
            exon_starts = sorted(exon_starts, reverse=True)
            exon_ends = sorted(exon_ends, reverse=True)
            for i in range(len(site_list1)):
                item=site_list1[i][0]
                gRNA=site_list1[i][1]
                try:
                    for i in range(len(exon_starts)):
                        if item >= exon_starts[i] and item >= exon_ends[i]:
                            out_exon.append(item)
                            output_list.append("%s位点%s位于基因%s的内含子中"%(item,gRNA,key))
                            continue
                        elif item >= exon_starts[i] and item <= exon_ends[i]:
                            if item + seq_len >= exon_ends[i]:
                                exon_intron.append(item)
                                output_list.append("%s位点%s位于基因%s第%d个外显子与其相邻的内含子上" % (item,gRNA, key, (len(exon_starts) - i)))
                                continue
                            else:
                                in_exon.append(item)
                                output_list.append("%s位点%s位于基因%s第%d个外显子上" % (item, gRNA,key, (len(exon_starts) - i)))
                                continue
                        elif item <= exon_starts[i] and item + seq_len >= exon_starts[i]:
                            intron_exon.append(item)
                            output_list.append("%s位点%s位于基因%s内含子与其相邻第%d个外显子上" % (item,gRNA, key, len(exon_starts) - i))
                            continue
                        else:
                            continue
                except:
                    continue
        return output_list
    def write_file(out,score,anno):
        with open(out,"w")as output_file:
            output_file.write("gRNA score ,\"green\" means good and \"red\" means bad:\r\n")
            for c in score:
                c="    ".join(str(x) for x in c[0])+"   "+c[1]
                output_file.write(c+"\r\n")
            output_file.write("gRNA annotation:\r\n")
            for c in anno:
                output_file.write(c+"\r\n")
            output_file.close()
if __name__ == "__main__":#运行
    app = QtWidgets.QApplication(sys.argv)
    window = MainUi()
    window.show()
    sys.exit(app.exec_())