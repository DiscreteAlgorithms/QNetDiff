
class UnionFind:
    def __init__(self, array) -> None:
        self.parent = dict()
        self.rank = dict()
        for x in array:
            self.parent[x] = x
            self.rank[x] = 0
    
    def find(self, x):
        if self.parent[x] == x:
            return x
        else:
            self.parent[x] = self.find( self.parent[x] )
            return self.parent[x]

    def unite(self, x, y):
        x = self.find(x)
        y = self.find(y)
        if x==y:
            return
        if self.rank[x] < self.rank[y]:
            self.parent[x] = y
        else:
            self.parent[y] = x
            if self.rank[x] == self.rank[y]:
                self.rank[x] += 1

    def isSame(self, x, y):
        return self.find(x) == self.find(y)


            