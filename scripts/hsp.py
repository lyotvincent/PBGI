class Hsp:

    qStart = 0
    qEnd = 0
    sStart = 0
    sEnd = 0 
    evalue = ""
    score = ""

    def __init__(self):
        self.qStart = 0
        self.qEnd = 0
        self.sStart = 0
        self.sEnd = 0 
        self.evalue = ""
        self.score = ""
    
    def getHsp(self):
        return "q.start:"+str(self.qStart)+"\t"+"q.end:"+str(self.qEnd)+"\t"+"s.start:"+str(self.sStart)+"\t"+"s.end:"+str(self.sEnd)+"\t"+"evalue:"+str(self.evalue)+"\t"+"bit-score:"+self.score