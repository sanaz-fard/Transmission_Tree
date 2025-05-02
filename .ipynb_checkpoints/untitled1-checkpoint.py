def test_generated_trees():

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