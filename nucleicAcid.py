from vpython import *
import math as m
class RNA:
    def __init__(self, base = ''):
        self.base = base.upper()
        self.pair = {'A' : 'U', 'U' : 'A', 'G' : 'C', 'C' : 'G'}
        self.codon = {'UUU':'Phe', 'UUC':'Phe', 'UUA':'Leu', 'UUG':'Leu', 
                      'CUU':'Leu', 'CUC':'Leu', 'CUA':'Leu', 'CUG':'Leu', 
                      'AUU':'Ile', 'AUC':'Ile', 'AUA':'Ile', 'AUG':'Met', 
                      'GUU':'Val', 'GUC':'Val', 'GUA':'Val', 'GUG':'Val', 
                      'UCU':'Ser', 'UCC':'Ser', 'UCA':'Ser', 'UCG':'Ser', 
                      'CCU':'Pro', 'CCC':'Pro', 'CCA':'Pro', 'CCG':'Pro', 
                      'ACU':'Thr', 'ACC':'Thr', 'ACA':'Thr', 'ACG':'Thr', 
                      'GCU':'Ala', 'GCC':'Ala', 'GCA':'Ala', 'GCG':'Ala', 
                      'UAU':'Tyr', 'UAC':'Tyr', 'UAA':'STOP', 'UAG':'STOP', 
                      'CAU':'His', 'CAC':'His', 'CAA':'Gln', 'CAG':'Gln', 
                      'AAU':'Asn', 'AAC':'Asn', 'AAA':'Lys', 'AAG':'Lys', 
                      'GAU':'Asp', 'GAC':'Asp', 'GAA':'Glu', 'GAG':'Glu', 
                      'UGU':'Cys', 'UGC':'Cys', 'UGA':'STOP', 'UGG':'Trp', 
                      'CGU':'Arg', 'CGC':'Arg', 'CGA':'Arg', 'CGG':'Arg', 
                      'AGU':'Ser', 'AGC':'Ser', 'AGA':'Arg', 'AGC':'Arg', 
                      'GGU':'Gly', 'GGC':'Gly', 'GGA':'Gly', 'GGG':'Gly'}
    def __str__(self):
        return f"RNA 3'-{self.base}-5'"
    def complement(self): #상보적 RNA 객체 생성
        string = list(map(lambda x : self.pair[x], (self.base[::-1])))
        return RNA(''.join(string))
    def transcription(self): #전사, mRNA 객체 생성
        return mRNA(self.complement().base)
    def simpleTranslation(self): #단순 번역, base를 가공없이 바로 단백질로 번역
        def split(self):
            li = list(map(''.join, zip(*[iter(self.base)]*3)))
            if len(li[-1]) < 3:
                li.pop()
            return li
        li = list(map(lambda x : self.codon[x], split(self)))
        return li
    def translation(self): #번역, mRNA로 만든 후, AUG에서 시작, 종결코돈에서 종결
        mrna = self.transcription()
        if not 'AUG' in mrna.base:
            return ''
        else:
            start = RNA(mrna.base[mrna.base.index('AUG'):])
            li = start.simpleTranslation()
        if 'STOP' in li:
            li = li[:li.index('STOP')]
        print('-'.join(li))
    def show(self):
        c = {'U': color.red, 'C':color.green, 'G':color.yellow, 'A':color.blue}
        term = 2
        r = 1
        index = 0
        text(text="3'", pos=vector(-r-1, 0,0), align = 'center')
        for i in self.base:
            sphere(pos=vector(index, 0, 0), color = c[i], radius = r)
            text(text = i, pos = vector(index, r+1, 0), align = 'center')
            index += term
        text(text="5'", pos=vector(index-term+r+1, 0 ,0), align = 'center')
    def ligase(self, string):
        self.base = self.base + string.upper()
        
class mRNA(RNA):
    def __str__(self):
        return "mRNA 3'-" + self.base + "-5'"
    def translation(self):
        if not 'AUG' in self.base:
            return ''
        else:
            start = RNA(self.base[self.base.index('AUG'):])
            li = start.simpleTranslation()
        if 'STOP' in li:
            li = li[:li.index('STOP')]
        print('-'.join(li))
    def transcription(self):
        raise Exception("Can't transcript mRNA")
    def complement(self):
        reverse = lambda x : self.pair[x]
        string = list(map(reverse, (self.base[::-1])))
        return mRNA(''.join(string))

class DNA(RNA):
    def __init__(self, base = ''):
        self.base = base.upper()
        self.pair = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
        self.codon = {'TTT':'Phe', 'TTC':'Phe', 'TTA':'Leu', 'TTG':'Leu', 
                      'CTT':'Leu', 'CTC':'Leu', 'CTA':'Leu', 'CTG':'Leu', 
                      'ATT':'Ile', 'ATC':'Ile', 'ATA':'Ile', 'ATG':'Met', 
                      'GTT':'Val', 'GTC':'Val', 'GTA':'Val', 'GTG':'Val', 
                      'TCT':'Ser', 'TCC':'Ser', 'TCA':'Ser', 'TCG':'Ser', 
                      'CCT':'Pro', 'CCC':'Pro', 'CCA':'Pro', 'CCG':'Pro', 
                      'ACT':'Thr', 'ACC':'Thr', 'ACA':'Thr', 'ACG':'Thr', 
                      'GCT':'Ala', 'GCC':'Ala', 'GCA':'Ala', 'GCG':'Ala', 
                      'TAT':'Tyr', 'TAC':'Tyr', 'TAA':'STOP', 'TAG':'STOP', 
                      'CAT':'His', 'CAC':'His', 'CAA':'Gln', 'CAG':'Gln', 
                      'AAT':'Asn', 'AAC':'Asn', 'AAA':'Lys', 'AAG':'Lys', 
                      'GAT':'Asp', 'GAC':'Asp', 'GAA':'Glu', 'GAG':'Glu', 
                      'TGT':'Cys', 'TGC':'Cys', 'TGA':'STOP', 'TGG':'Trp', 
                      'CGT':'Arg', 'CGC':'Arg', 'CGA':'Arg', 'CGG':'Arg', 
                      'AGT':'Ser', 'AGC':'Ser', 'AGA':'Arg', 'AGC':'Arg', 
                      'GGT':'Gly', 'GGC':'Gly', 'GGA':'Gly', 'GGG':'Gly'}
    def __str__(self):
        return "DNA 3'-" + self.base + "-5'"
    def show(self):
        c = {'T': color.red, 'C':color.green, 'G':color.yellow, 'A':color.blue}
        term = 2 # 염기 간 간격
        r = 1 # 염기반지름
        R = 3 # 회전 반지름
        index = 0 # 초기 위치
        angle = 36 # 1회 회전량
        theta = 0 # 초기 각도
        text(text="3'", pos=vector(-r-1, 0,0), align = 'center')

        before1 = vector(index, R*m.sin(m.radians(theta)), R*m.cos(m.radians(theta)))
        before2 = vector(index, -R*m.sin(m.radians(theta)), -R*m.cos(m.radians(theta)))
        for i in self.base:
            y = R*m.sin(m.radians(theta))
            z = R*m.cos(m.radians(theta))

            sphere(pos=vector(index, y, z), color = c[i], radius = r)
            label(text = i, pos = vector(index, y+2+r, z), align = 'center')

            j = self.pair[i]
            sphere(pos=vector(index, -y, -z), color = c[j], radius = r)
            label(text = j, pos = vector(index, -y-2-r, -z), align = 'center')

            curve(vector(index, y, z), vector(index, -y, -z))
            curve(before1, vector(index, y, z))
            curve(before2, vector(index, -y, -z))

            before1 = vector(index, y, z)
            before2 = vector(index, -y, -z)

            index += term
            theta += angle
        text(text="5'", pos=vector(index-term+r+1, 0 ,0), align = 'center')