
class Query:
    def __init__(self, name, snps):
        self.name = name
        self.snps = snps
        self.lineage = ""
        self.identity = 0
        self.prop_snps_included = 0
        
    def update_lineage(self, lineage):
        self.lineage = lineage[0]
        self.identity = lineage[1]
        self.prop_snps_included = lineage[2]

class Lineage:
    def __init__(self, name, lineages):
        self.name = name
        
        self.snps = []
        self.ancestors = []
        if self.name != "root":
            self.parent = self.assign_parent(lineages)
            self.get_ancestors(lineages)
        self.children = []

    def get_children_names(self):
        print([i.name for i in self.children])
    
    def assign_parent(self, lineages):
        
        name_list = self.name.split(".")
        
        parent_lineage = ""
        if len(name_list) == 1:
            if self.name == "A":
                parent_lineage = "root"
            elif self.name == "B":
                parent_lineage = "A"
        else:
            parent_lineage = '.'.join(self.name.split(".")[:-1])
            
        if parent_lineage in lineages:
            siblings = lineages[parent_lineage].children
            siblings.append(self)
            return lineages[parent_lineage]
            
        else:
            new_parent = Lineage(parent_lineage,lineages)
            siblings = new_parent.children
            siblings.append(self)
            lineages[new_parent.name]=new_parent
            return new_parent
        
    def get_parent_snps(self,snp):
        snps = []
        parent= self.parent
        for snp in parent.snps:
            snps.append(snp)

        if parent.name != "root":
            parent.get_parent_snps(snp)
        return list(set(snps))
    
    def add_snp(self, snp):
        if snp not in self.snps:
            parents_snps = self.get_parent_snps(snp)
            if snp not in parents_snps:
                self.snps.append(snp)
