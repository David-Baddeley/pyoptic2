import sources

# Complete optical system
class System(list) : 

    def __init__(self, iterable) :
        list.__init__(self, iterable)

        # maximum size as determined from the elements
        # FIXME - is this used anywhere?
        self.size = [0,0,0]

        
    def propagate(self, source=None) :
        raytree = []
        
        if source is None:
            source = self[0]
            
        if isinstance(self[0],sources.Source):
            startAt = 1
        else:
            startAt = 0

        # iterate over rays from source
        i = 0 
        ri = iter(source)
        try :                
            while True :
                r = ri.next()
                raybranch = self.propagate_ray(r, startAt=startAt, source=source)
                raytree.append(raybranch)
                i += 1
        except StopIteration :
            pass
            

        return raytree
        
    def propagate_ray(self, r, startAt=1, source=None):
        raybranch = []
        raybranch.append(r)
        for j in range(startAt,len(self)) :
            #print 'System.propagate> element=',j
            if j == 0:
                rp = self[j].propagate(source, raybranch[j-startAt])
            else:
                rp = self[j].propagate(self[j-1],raybranch[j-startAt])
                
            if rp == None:
                break                   
            raybranch.append(rp)
            
        return raybranch
    
    def prepend(self, elements):
        try:
            for e in elements[::-1]:
                self.insert(0, e)
        except AttributeError:
            #single element
            self.insert(0, elements)
            
    def add(self, elements):
        try:
            self.extend(elements)
        except TypeError:
            self.append(elements)
                

                
    def __str__(self) :
        s = ''
        for e in self :
            s +=str(e)
        return s


