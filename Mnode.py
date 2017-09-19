import csv
import random


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
    node_id = "0"

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

    def set_basic_information(self, node_id, mutation_name, ccf):
        self.node_id = node_id
        self.mutation_name = mutation_name
        self.ccf = ccf

    def set_advance_information(self, cnv_N, cnv_T, genotype):
        self.cnv_N = cnv_N
        self.cnv_T = cnv_T
        self.genotype = genotype

    def get_information(self):
        return [self.node_id, self.genotype, self.mutation_name, self.cnv_N, self.cnv_T]

class TEHtree:
    root = None
    full_information = False
    name_pool = []
    def __init__(self, file = ""):
        if len(file) == 0:
            #randomly build a tree
            self.name_pool = self.buildNamePoor()
            self.buildRandomTree()
        else:
            #read the tree structure from a csv
            self.readTree(file)

    def getName(self):
        return self.name_pool.pop(0)

    def buildRandomTree(self):
        mroot = Mnode()
        mroot.set_basic_information("1", self.getName(), 100)
        self.root = mroot
        self.recRandomBuild(mroot)

    def recRandomBuild(self, node):
        ccp = node.ccf
        if ccp > 1:
            left = random.randint(1, ccp/2)
            right = ccp - left
            left_node = Mnode()
            left_node.set_basic_information(node.node_id+"0", self.getname(), left)
            right_node = Mnode()
            right_node.set_basic_information(node.node_id + "1", self.getname(), right)
            self.recRandomBuild(left_node)
            self.recRandomBuild(right_node)

    def readTree(self, file):
        with open(file, 'rb') as msf:
            node_dict = {}
            csvreader = csv.reader(msf)
            Title = csvreader.next()
            if "CNV_N" in Title:
                self.full_information = True
            for line in csvreader:
                node_id = line[0]
                ccf = int(line[-1])
                mutation_name = line[2]
                tmp_node = Mnode()
                tmp_node.set_basic_information(node_id, mutation_name, ccf)
                if self.full_information:
                    tmp_node.set_advance_information(float(line[3]), float(line[4]), int(line[1]))
                node_dict[node_id] = tmp_node
        self.recBuild("1", node_dict)
        self.root = node_dict["1"]

    def recBuild(self, id, node_dict):
        left_id = id+"0"
        right_id = id+"1"
        tmp_node = node_dict.get(left_id, None)
        node_dict[id].set_child(tmp_node)
        if tmp_node != None:
            self.recBuild(left_id, node_dict)
        tmp_node = node_dict.get(right_id, None)
        node_dict[id].set_child(tmp_node, False)
        if tmp_node != None:
            self.recBuild(right_id, node_dict)

    def outputCSV(self, fileName):
        title = ["Node_id", "Genotype", "Mutation_Name", "CNV_N", "CNV_T", "CCF"]
        output_list = []
        self.recOut(self.root, output_list)
        sorted(output_list, key= lambda x: int(x[0]))
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
            self.recount(node.right_child, output_list)