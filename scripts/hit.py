class Hit:
    
    desc = ''
    id = ''
    length = 0
    hsps = []

    def __init__(self):
        self.id = ''
        self.desc = ''
        self.length = 0
        self.hsps = []
    
    def getCoverage(self):
        
        t = [0 for i in range(self.length+1)]

        for h in self.hsps:
            if h.sStart < h.sEnd:
                d = h.sStart
                f = h.sEnd
            else:
                f = h.sStart
                d = h.sEnd
            for i in range(d, f+1):
                t[i] = 1
        tot = 0
        for i in range(len(t)):
            tot += t[i]
        
        return tot/(self.length*100)

    def getHit(self):
        res="*****************************************************************\n"
        res += 'id:'+self.id+'\n'
        res += self.desc+'\n'
        res += 'HSP_length:'+str(self.length)+'\n'
        res += 'HSP_coverage:'+str(self.getCoverage())+'%\n'
        for i in range(len(self.hsps)):
            h = self.hsps[i]
            res += 'HSP_'+str(i+1)+'\t'+h.getHsp()+'\n'
        return res
