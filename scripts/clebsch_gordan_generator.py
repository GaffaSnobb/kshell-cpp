from sympy.physics.quantum.cg import CG
import numpy as np

def main():
    arr = np.zeros((3, 5, 3, 5, 5))
    print(f"{arr.ndim = }")
    print(f"{arr.size = }")
    print(f"{arr.nbytes = }")
    m_vals = {
        3: [-3, -1, 1, 3],
        5: [-5, -3, -1, 1, 3, 5],
        1: [-1, 1],
    }
    counter = 1
    with open("cg_flat_arr.txt", "w") as outfile:
        for j1 in [3, 5, 1]:
            for m1 in m_vals[j1]:

                for j2 in [3, 5, 1]:
                    for j3 in range(abs(j1 - j2), j1 + j2 + 2, 2):
                        for m3 in range(-j3, j3 + 2, 2):
                            for m2 in m_vals[j2]:

                                if (m1 + m2) != m3: continue
                                if counter == 3:
                                    msg = ",\n"
                                    counter = 0
                                else:
                                    msg = ", "
                                
                                cg_creation = CG(
                                    j1 = j1/2,
                                    m1 = m1/2,
                                    j2 = j2/2,
                                    m2 = m2/2,
                                    j3 = j3/2,
                                    m3 = m3/2,
                                )
                                # outfile.write(f"{float(cg_creation.doit()):19.16f}" + msg)
                                counter += 1
                                # arr[j1, m1, j2, m2, j3, m3] = float(cg_creation.doit())
                                # outfile.write(f"({j1:2}, {m1:2}, {j2:2}, {m2:2}, {j3:2}, {m3:2}): {float(cg_creation.doit())},\n")
                                # outfile.write(r"{{" + f"{j1:2}, {m1:2}, {j2:2}, {m2:2}, {j3:2}, {m3:2}" + r"}, " + f"{float(cg_creation.doit())}" + r"}" + ",\n")
                                # outfile.write(r"{" + f"{j1:2}, {m1:2}, {j2:2}, {m2:2}, {j3:2}, {m3:2}" + r"}," + "\n")
                                # print(r"{" + f"{j1:2}, {(m1 + j1)//2:2}, {j2:2}, {(m2 + j2)//2:2}, {j3//2:2}" + r"},")

def main2():
    arr = np.zeros(shape=(3, 6, 3, 6, 6))
    arr_flat = np.zeros(1944)
    print(f"{arr.ndim = }")
    print(f"{arr.size = }")
    print(f"{arr.nbytes = }")
    m_vals = {
        1: [-1, 1],
        3: [-3, -1, 1, 3],
        5: [-5, -3, -1, 1, 3, 5],
    }
    j_vals = [1, 3, 5]
    counter = 1

    flat_idx = 0

    for j1_idx in range(len(j_vals)):

        m1_vals = m_vals[j_vals[j1_idx]]
        # m1_vals = m_vals[5]
        for m1_idx in range(len(m1_vals)):

            # arr[j1_idx, m1_idx] = (j1_idx, m1_idx)

            for j2_idx in range(len(j_vals)):

                m2_vals = m_vals[j_vals[j2_idx]]
                # m2_vals = m_vals[5]
                for m2_idx in range(len(m2_vals)):
                    
                    j3_vals = range( abs(j_vals[j1_idx] - j_vals[j2_idx]), j_vals[j1_idx] + j_vals[j2_idx] + 2, 2 )

                    # j_min = max(
                    #     abs(orbital_idx_to_j_map[creation_orb_idx_0] - orbital_idx_to_j_map[creation_orb_idx_1]),
                    #     std::abs(orbital_idx_to_j_map[annihilation_orb_idx_0] - orbital_idx_to_j_map[annihilation_orb_idx_1])
                    # );
                    # unsigned short j_max = std::min(
                    #     orbital_idx_to_j_map[creation_orb_idx_0] + orbital_idx_to_j_map[creation_orb_idx_1],
                    #     orbital_idx_to_j_map[annihilation_orb_idx_0] + orbital_idx_to_j_map[annihilation_orb_idx_1]
                    # );

                    for j3_idx in range(len(j3_vals)):
                        """
                        NOTE: m3 is not explicitly needed since m1 + m2 = m3 for the CG coeff to be non-zero.
                        """
                        j1 = j_vals[j1_idx]/2
                        m1 = m1_vals[m1_idx]/2
                        j2 = j_vals[j2_idx]/2
                        m2 = m2_vals[m2_idx]/2
                        j3 = j3_vals[j3_idx]/2
                        m3 = m1 + m2
                        
                        if counter == 3:
                            msg = ",\n"
                            counter = 0
                        else:
                            msg = ", "
                        
                        cg_creation = CG(
                            j1 = j1,
                            m1 = m1,
                            j2 = j2,
                            m2 = m2,
                            j3 = j3,
                            m3 = m3,
                        )
                        
                        # outfile.write(f"{float(cg_creation.doit()):19.16f}" + msg)
                        counter += 1


                        # arr[j1_idx, m1_idx, j2_idx, m2_idx, j3_idx] = (j1_idx, m1_idx, j2_idx, m2_idx, j3_idx)
                        # arr[j1_idx, m1_idx, j2_idx, m2_idx, j3_idx] = float(cg_creation.doit())
                        calculated_flat_idx = ((((j1_idx*6 + m1_idx)*3 + j2_idx)*6 + m2_idx)*6 + j3_idx)    # Indexing a 1D array as if it was a 5D array.
                        arr_flat[calculated_flat_idx] = float(cg_creation.doit())

                        if (j1_idx == 1) and (m1_idx == 0) and (j2_idx == 2) and (m2_idx == 4) and (j3_idx == 2):
                            print(f"{j1_idx = }")
                            print(f"{m1_idx = }")
                            print(f"{j2_idx = }")
                            print(f"{m2_idx = }")
                            print(f"{j3_idx = }")
                            print(f"{calculated_flat_idx = }")
                            print(f"{arr_flat[calculated_flat_idx] = }")
                            print(f"{j_vals[j1_idx] = }")
                            print(f"{m1_vals[m1_idx] = }")
                            print(f"{j_vals[j2_idx] = }")
                            print(f"{m2_vals[m2_idx] = }")
                            print(f"{j3_vals[j3_idx] = }")
                            print(f"{j3 = }")
                            print(f"{list(j3_vals) = }")
                            print(f"{[elem//2 for elem in j3_vals] = }")
                        flat_idx += 1
                        # arr[j1_idx, m1_idx, j2_idx, m2_idx, j3_idx] = 99.9
                # break
    
    msg = ", "
    with open("cg_flat_arr.txt", "w") as outfile:
        for i in range(len(arr_flat)):

            if (i + 1)%3 == 0:
                msg = ",\n"
            else:
                msg = ", "

            outfile.write(f"{arr_flat[i]:21.16f}" + msg)

    print(len(arr_flat))
    print(len(arr_flat[arr_flat == 0]))

if __name__ == "__main__":
    main2()