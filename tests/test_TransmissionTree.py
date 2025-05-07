from codes import TransmissionTree

def test_generated_trees():
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

    sim_obj_A = TransmissionTree(0, None, None, False, test_A)
    sim_obj_A.final()
    # Compare simple structure
    assert sim_obj_A.phylo_tree.write(format=1) == t_A.write(format=1)

    sim_obj_B = TransmissionTree(0, None, None, False, test_B)
    sim_obj_B.final()
    # Compare simple structure
    assert sim_obj_B.phylo_tree.write(format=1) == t_B.write(format=1)

    sim_obj_C = TransmissionTree(0, None, None, False, test_C)
    sim_obj_C.final()
    # Compare simple structure
    assert sim_obj_C.phylo_tree.write(format=1) == t_C.write(format=1)