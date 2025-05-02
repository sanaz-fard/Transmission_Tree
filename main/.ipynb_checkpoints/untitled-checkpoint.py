import csv
import uuid
import numpy as np
import covasim as cv
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from ete3 import Tree, TreeStyle, TextFace, NodeStyle

class TransmissionTree:
    def __init__(self, population=0, type=None, loc=None, to_networkx=False, test_case = None):
        self.population = population
        self.type = type
        self.loc = loc
        self.to_networkx = to_networkx
        self.pars = dict(pop_size=self.population, 
                         pop_type=self.type, 
                         location=self.loc)
        self.test_case = test_case

    def simulations(self):
        sim = cv.Sim(self.pars).run()
        self.trans_tree = sim.make_transtree(to_networkx = self.to_networkx)

    def graph_vis(self):
        nx.draw(self.trans_tree.graph, with_labels = True, pos=nx.spring_layout(self.trans_tree.graph))
        plt.show()

    def transmissions_info(self):
        if self.test_case is None:
            trans_jason = self.trans_tree.to_json()
            trans_log = trans_jason['infection_log']
            self.transmissions = []  ### records (parent, child, time) 
            for i in range(len(trans_log)):
                if trans_log[i]['source'] == None: continue ### this line disregards seed infections as they don't have a parent
                temp = [int(trans_log[i]['source']), int(trans_log[i]['target']), int(trans_log[i]['date'])]
                self.transmissions.append(temp)
        else:
            self.transmissions = self.test_case

    def save_pop(self):
        if self.test_case is None:
            with open('data/'+'transmission_pop_'+str(self.population)+'.csv', 'w') as csvfile:
                spamwriter = csv.writer(csvfile)
                spamwriter.writerow(['Parent', 'Child', 'Origin_time'])
                for i in self.transmissions:
                    spamwriter.writerow([i[0], i[1], i[2]])
        else: pass

    def acyclic(self):                                  ##### if 2 enteries have the same child and origin time,  the first appeared parent
        grouped = defaultdict(list)                     ##### in order to have only 1 parent for each node,                                            
        for item in self.transmissions:                   ##### will be taken and the rest will be removed
            key = tuple(item[-2:])                 
            grouped[key].append(item)                                          ##### Step 1: Group by the last two items
        self.acyclic_transmissions = [group[0] for group in grouped.values()]  ##### Step 2: Keep only one representative from each group
        print(self.acyclic_transmissions)

    def fix_reinfections(self):                       ##### Sort by the last item of each sublist,
        data = self.acyclic_transmissions                       ##### origin_time, so parents who will appear
        sorted_data = sorted(data, key=lambda x: x[-1]) ##### as a child in the future will get a new
        parents = []                                    ##### name based on its node name and the new time of origin
        for i in sorted_data:
            if i[0] not in parents:
                parents.append(i[0])
            if i[1] in parents:
                i[1] = str(i[1])+'_'+str(i[-1])
        self.fixed_reinf = sorted_data
        print(sorted_data)

    def list_to_newick(self):
        # Step 1: Build children map and origin_times, with duplicate handling
        children = defaultdict(list)
        origin_times = {}
        parent_history = set()

        for parent, child, time in self.fixed_reinf:
            if child in parent_history and child != parent:
                # Duplicate if child is already a parent
                child = f"{child}_dup_{uuid.uuid4().hex[:6]}"
            parent_history.add(parent)
            children[parent].append(child)
            origin_times[(parent, child)] = float(time)

        # Step 2: Identify roots
        all_parents = set(children.keys())
        all_children = {c for clist in children.values() for c in clist}
        roots = list(all_parents - all_children)

        # Step 3: Add dummy root if needed
        if len(roots) > 1:
            dummy_root = "ROOT"
            for r in roots:
                children[dummy_root].append(r)
                origin_times[(dummy_root, r)] = 0.0
            root = dummy_root
        else:
            root = roots[0]

        # Step 4: Recursive Newick generator
        def to_newick(node):
            if node not in children:
                return f"{node}"
            subtrees = [to_newick(child) + f":{origin_times[(node, child)]}"
                        for child in children[node]]
            return f"({','.join(subtrees)}){node}"

        return to_newick(root) + ";"

    def simple_tree_vis(self):
        self.newick_str = self.list_to_newick()
        self.phylo_tree = Tree(self.newick_str, format=1)                  ##### Load tree
        self.phylo_tree.show()                                        ##### Show the tree in a GUI window
        
    def style_all_nodes(self):
        tree = self.phylo_tree
        for node in tree.traverse():
            # Show label on all nodes
            name_face = TextFace(node.name, fsize=10)
            node.add_face(name_face, column=0, position="branch-right")

            # Node style
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "black"
            nstyle["size"] = 10
            nstyle["shape"] = "circle"

            # Color code: leaves = blue, internal = red
            if node.is_leaf():
                nstyle["bgcolor"] = "lightblue"
            else:
                nstyle["bgcolor"] = "lightcoral"

            node.set_style(nstyle)

    def get_tree_style(self):
        ts = TreeStyle()
        ts.show_leaf_name = False  # Prevent auto-labeling
        ts.show_branch_length = True
        ts.show_branch_support = False
        ts.scale = 100
        self.ts = ts

    def styled_vis(self):
        self.style_all_nodes()
        self.phylo_tree.show(tree_style=self.ts)

    def all_labeled(self):
        for node in self.phylo_tree.traverse():
            name_face = TextFace(node.name, fsize=10)
            node.add_face(name_face, column=0, position="branch-right")
        self.ts.branch_vertical_margin = 10
        self.phylo_tree.show(tree_style=self.ts)

    def final(self):
        if self.test_case is None:
            self.simulations()
            if self.to_networkx:
                graph_vis()
            self.transmissions_info()
            self.save_pop()
            self.acyclic()
            self.fix_reinfections()
            self.list_to_newick()
            self.simple_tree_vis()
            self.get_tree_style()
            self.styled_vis()
            self.all_labeled()
        else:
            self.transmissions_info()
            self.save_pop()
            self.acyclic()
            self.fix_reinfections()
            self.list_to_newick()
            self.simple_tree_vis()
            self.get_tree_style()
            self.styled_vis()
            self.all_labeled()

if __name__ == __main__:
    test_A = [                                                                   ##### newick: (((X_3)Z)Y)X;
    ['X', 'Y', 1],
    ['Y', 'Z', 2],
    ['Z', 'X', 3],  # X reappears as child at time 3 → becomes X_3
    ]

    test_B = [                                                                 ##### newick: ((((((T)S)Q)P)N, O)M, (M_5, U)R)ROOT;
    ['M', 'N', 1],
    ['M', 'O', 2],
    ['N', 'P', 3],
    ['O', 'P', 3],  # P has two parents at same time: N (kept), O (ignored)
    ['P', 'Q', 4],
    ['R', 'M', 5],  # M becomes a child of R at time 5 → M_5
    ['Q', 'S', 6],
    ['S', 'T', 7],
    ['R', 'U', 8],
    ['U', 'T', 7],] # T has two parents at same time: S (kept), U (ignored)


    test_C = [                                                                 ##### NEWICK: (((((((J)I)G)F)D)B, (E)C)A, ((((O)N)M)L)K)H;
    ['A', 'B', 1],
    ['A', 'C', 1],                                                         ##### (((B(D(F(G(I(J))))),C(E))A_8,K(L(M(N(O)))))H);
    ['B', 'D', 2],
    ['C', 'E', 3],
    ['D', 'F', 4],
    ['E', 'F', 4],  # F has two parents at same time: D (kept), E (ignored)
    ['F', 'G', 5],
    ['H', 'A', 6],  # A becomes child at time 6 → A_6
    ['G', 'I', 7],
    ['I', 'J', 8],
    ['H', 'K', 9],
    ['K', 'L', 10],
    ['J', 'L', 10],  # L has two parents at same time: K (kept), J (ignored)
    ['L', 'M', 11],
    ['M', 'N', 12],
    ['N', 'O', 13],
    ]

    # Create a tree manually
    t_A = Tree()
    root_A = t_A.add_child(name="root")

    # Add children to the root
    d = root_A.add_child(name="X")
    d1 = d.add_child(name="Y")
    d2 = d1.add_child(name="Z")
    d3 = d2.add_child(name="X_3")

    # Create a tree manually
    t_B = Tree()
    root_B = t_B.add_child(name="root")

    # Add children to the root
    a = root_B.add_child(name="M")
    b = root_B.add_child(name="R")

    # Add children to A
    a1 = a.add_child(name="N")
    a2 = a.add_child(name="O")
    a3 = a1.add_child(name="P")
    a4 = a3.add_child(name="Q")
    a5 = a4.add_child(name="S")
    a6 = a5.add_child(name="T")
    a7 = b.add_child(name="M_5")
    a8 = b.add_child(name="U")

    # Create a tree manually
    t_C = Tree()
    root_C = t_C.add_child(name="root")

    # Add children to the root
    c1 = root_C.add_child(name="A")
    c2 = root_C.add_child(name="H")

    c3 = c1.add_child(name="B")
    c4 = c3.add_child(name="D")
    c5 = c4.add_child(name="F")
    c6 = c5.add_child(name="G")
    c7 = c6.add_child(name="I")
    c8 = c7.add_child(name="J")
    c9 = c8.add_child(name="L")
    c10 = c9.add_child(name="M")
    c11 = c10.add_child(name="N")
    c12 = c11.add_child(name="O")

    c13 = c1.add_child(name="C")
    c14 = c13.add_child(name="E")

    c15 = c2.add_child(name="K")
    c16 = c2.add_child(name="A_6")