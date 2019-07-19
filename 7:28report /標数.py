class Fp():
    def __init__(self, p):
        self.p = p
        self.Set = set(i for i in range(p))
    def sum(self, a, b):
        if (a in self.Set) and (b in self.Set):
            return (a+b)%self.p
    def multplctn(self, a,b):
        if (a in self.Set) and (b in self.Set):
            return (a*b)%self.p
