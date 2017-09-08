class Mnode:
    mutation_name = ""
    #cancer cell fraction
    ccf = 0
    #normal copy number
    cnv_N = 2
    #tumor copy number
    cnv_T = 2
    #genotype -- 0 as Heterozygous ; 1 as homozygous
    genotype = 0
    parent = None
    left_child = None
    right_child = None
    children_number = 0

    def __init__(self):
        self.mutation_name = ""

    def set_parent(self, parent):
        self.parent = parent

    def update_children_number(self):
        self.children_number = 0
        if self.left_child != None:
            self.children_number += 1
        if self.right_child != None:
            self.children_number += 1

    def set_child(self, child, left = True):
        if child != None:
            if left:
                left_child = child
            else:
                right_child = child
            child.set_parent(self)
            self.update_children_number()


