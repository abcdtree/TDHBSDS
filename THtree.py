import csv
import random
import pickle

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    with open(filename, 'rb') as input:
        return pickle.load(input)

def sort_a(x):
    a1, a2 = x.split(":")
    try:
        a1 = int(a1)
    except:
        if a1 == "X":
            a1 = 23
        else:
            a1 = 24
    return (a1, int(a2))


class Mutation:
    chrome = "1"
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
    all_muts = []

    def __init__(self, file=""):
        if len(file) == 0:
            self.buildMutations()
            print len(self.mutation_pool)
            self.name_pool = self.buildNamePool()
            self.buildRandomTree()
        else:
            self.buildFromFile(file)

    def buildFromFile(self, file):
        with open(file, "rb") as csf:
            csvreader = csv.reader(csf)
            title = csvreader.next()
            node_dict = {}
            for line in csvreader:
                node_id = line[0]
                ccf = int(line[2])
                node_name = line[1]
                tmp_node = Node()
                tmp_node.set_id(node_id)
                tmp_node.set_ccf(ccf)
                tmp_node.set_name(node_name)
                muts = []
                m_strs = line[-1].split("+")
                for m_str in m_strs:
                    chrome = m_str.split("_")[0][3:]
                    loc = int(m_str.split("_")[1])
                    tmp_mut = Mutation(chrome, loc)
                    muts.append(tmp_mut)
                    self.all_muts.append(tmp_mut)
                tmp_node.set_mutations(muts)
                node_dict[node_id] = tmp_node
        self.recBuild("1", node_dict)
        self.root = node_dict["1"]

    def recBuild(self, id, node_dict):
        left_id = id + "0"
        right_id = id + "1"
        flag = True
        tmp_node = node_dict.get(left_id, None)
        node_dict[id].set_child(tmp_node)
        if tmp_node != None:
            self.recBuild(left_id, node_dict)
            flag = False
        tmp_node = node_dict.get(right_id, None)
        node_dict[id].set_child(tmp_node, False)
        if tmp_node != None:
            self.recBuild(right_id, node_dict)
            flag = False
        if flag:
            self.leaf_cell.append(node_dict[id])


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
        for i in range(1,23):
            chrome_list.append(str(i))
        chrome_list.append("X")
        chrome_list.append("Y")
        for chrome in chrome_list:
            loc_sample = []
            #loc_sample = random.sample(range(1,66666666), 100)
            for i in range(100):
                loc_sample.append(int(random.random() * 66666666))
            for loc in loc_sample:
                tmp_mut = Mutation(chrome, loc)
                self.mutation_pool.append(tmp_mut)


    def buildRandomTree(self):
        mroot = Node()
        mroot.set_id("1")
        mroot.set_ccf(100)
        mroot.set_name(self.getName())
        mroot.set_mutations(self.getMuts())
        self.all_muts = self.all_muts + mroot.get_mutations()
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
            self.all_muts = self.all_muts + left_node.get_mutations()
            node.set_child(left_node)
            right_node = Node()
            right_node.set_id(node.get_id() + "1")
            right_node.set_ccf(right)
            right_node.set_name(self.getName())
            right_node.set_mutations(self.getMuts())
            self.all_muts = self.all_muts + right_node.get_mutations()
            node.set_child(right_node, False)
            self.recRandomBuild(left_node)
            self.recRandomBuild(right_node)
        else:
            self.leaf_cell.append(node)

    def getMuts(self):
        number_of_mutations = int(random.gauss(7, 1))
        muts = []
        for i in range(number_of_mutations):
            k = random.randint(0, len(self.mutation_pool)-1)
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

    def subRemove(self, list_A, list_B):
        return [x for x in list_A if x not in list_B]

    def random_create_single_cell_data(self, fps, fns, sample_size):
        sample_cells = random.sample(self.leaf_cell, sample_size)
        single_cell_data = []
        for cell in sample_cells:
            positive_list = []
            c_node = cell
            while 1:
                if c_node == None:
                    break
                for mut in c_node.get_mutations():
                    positive_list.append(mut)
                c_node = c_node.parent
            num_of_muts = len(positive_list)
            num_of_fp = int(num_of_muts * fps)
            num_of_fn = int(num_of_muts * fns)
            other_muts = self.subRemove(self.all_muts[:] , positive_list)
            fp_list = random.sample(other_muts, num_of_fp)
            fn_list = random.sample(positive_list, num_of_fn)
            final_list = self.subRemove(positive_list,fn_list) + fp_list
            one_single_cell = []
            for mut in final_list:
                one_single_cell.append(str(mut.chrome) + ":" + str(mut.loc))
            one_single_cell.sort(key=sort_a)
            single_cell_data.append(one_single_cell)
        return single_cell_data

    def outSingleCells(self, fps, fns, sample_size, filename = "sample1.txt"):
        single_cell_data = self.random_create_single_cell_data(fps, fns, sample_size)
        with open(filename, 'w') as outfile:
            i = 1
            for cell in single_cell_data:
                outfile.write("Cell_" + str(i) + "\t" + "\t".join(cell) + "\n")
                i = i + 1



