import csv
import random

class Mutation:
    chrome = 1
    loc = 0
    cnv_R = 2
    cnv_RS = 2
    cnv_V = 2
    genotype = 0

    def __init__(self, chrome, loc):
        self.chrome = chrome
        self.loc = loc

    def set_cnv(self, cnv_R, cnv_RS, cnv_V):
        self.cnv_R = cnv_R
        self.cnv_RS = cnv_RS
        self.cnv_V = cnv_V

    def get_infor(self):
        return "chr" + str(self.chrome) + "_" + str(self.loc)

class Node:
    mutations = []
    ccf = 0
    node_name = ""
    parent = None
    left_child = None
    right_child = None
    node_id = "0"

    def __init__(self):
        self.mutations = []

    def set_parent(self, parent):
        self.parent = parent

    def set_child(self, child, left = True):
        if child != None:
            if left:
                self.left_child = child
            else:
                self.right_child = child
            child.set_parent(self)

    def add_mutation(self, a_mut):
        self.mutations.append(a_mut)

    def set_mutations(self, muts):
        for mut in muts:
            self.add_mutation(mut)

    def set_ccf(self, ccf):
        self.ccf = ccf

    def set_name(self, name):
        self.node_name = name

    def get_ccf(self):
        return self.ccf

    def get_mutations(self):
        return self.mutations[:]

    def get_id(self):
        return self.node_id

    def set_id(self, id):
        self.node_id = id

    def get_information(self):
        muts_list = []
        for mut in self.mutations:
            muts_list.append(mut.get_infor())
        return [self.node_id, self.node_name, self.ccf, "+".join(muts_list)]

class THTree:
    root = None
    mutation_pool = []
    name_pool = []
    leaf_cell = []

    def __init__(self):
        self.buildMutations()
        print len(self.mutation_pool)
        self.name_pool = self.buildNamePool()
        self.buildRandomTree()

    def buildNamePool(self):
        char_list = []
        for i in range(65, 91):
            char_list.append(chr(i))
        name_pool = char_list[:]
        for c in char_list:
            for b in char_list:
                name_pool.append(c + b)
        return name_pool

    def getName(self):
        return self.name_pool.pop(0)

    def buildMutations(self):
        chrome_list = []
        for i in range(1,24):
            chrome_list.append(str(i))
        chrome_list.append("X")
        chrome_list.append("Y")
        for chrome in chrome_list:
            loc_sample = random.sample(range(1,66666666), 100)
            for loc in loc_sample:
                tmp_mut = Mutation(chrome, loc)
                self.mutation_pool.append(tmp_mut)


    def buildRandomTree(self):
        mroot = Node()
        mroot.set_id("1")
        mroot.set_ccf(100)
        mroot.set_name(self.getName())
        mroot.set_mutations(self.getMuts())
        self.root = mroot
        self.recRandomBuild(mroot)

    def recRandomBuild(self, node):
        ccp = node.get_ccf()
        # print ccp
        if ccp > 1:
            left = random.randint(1, ccp / 2)
            right = ccp - left
            left_node = Node()
            left_node.set_id(node.get_id() + "0")
            left_node.set_ccf(left)
            left_node.set_name(self.getName())
            left_node.set_mutations(self.getMuts())
            node.set_child(left_node)
            right_node = Node()
            right_node.set_id(node.get_id() + "1")
            right_node.set_ccf(right)
            right_node.set_name(self.getName())
            right_node.set_mutations(self.getMuts())
            node.set_child(right_node, False)
            self.recRandomBuild(left_node)
            self.recRandomBuild(right_node)
        else:
            self.leaf_cell.append(node)

    def getMuts(self):
        number_of_mutations = int(random.gauss(7, 1))
        muts = []
        for i in range(number_of_mutations):
            k = random.randint(0, len(self.mutation_pool))
            tmp_mut = self.mutation_pool.pop(k)
            muts.append(tmp_mut)
        return muts

    def outputCSV(self, fileName):
        title = ["Node_id", "Node_Name", "CCF", "Mutations"]
        output_list = []
        # print self.root.left_child.mutation_name
        self.recout(self.root, output_list)
        # print output_list
        output_list.sort(key=lambda x: int(x[0]))
        with open(fileName, 'wb') as msf:
            spamwriter = csv.writer(msf, dialect="excel")
            spamwriter.writerow(title)
            for line in output_list:
                spamwriter.writerow(line)

    def recout(self, node, output_list):
        output_list.append(node.get_information())
        if node.left_child != None:
            self.recout(node.left_child, output_list)
        if node.right_child != None:
            self.recout(node.right_child, output_list)