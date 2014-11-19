indel        = -1 # used when creating table
new_indel    = -1
extend_indel = -1

# ATTENTION: To make sure that matrix is symmetric, do NOT change zeroes. They
#            will automatically edited.
similarity_matrix = [  # A   G   C   T
                       [ 1, -1, -1, -1],  # A
                       [ 0,  1, -1, -1],  # G
                       [ 0,  0,  1, -1],  # C
                       [ 0,  0,  0,  1],  # T
                    ]
