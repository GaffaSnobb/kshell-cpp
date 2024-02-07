from sympy.physics.quantum.cg import CG
import numpy as np

def main():
    arr = np.zeros(shape=(3, 6, 3, 6, 6))
    # arr = np.zeros(shape=(3, 6), dtype=tuple)
    arr_flat = np.zeros(1944)
    print(f"{arr.ndim = }")
    print(f"{arr.size = }")
    print(f"{arr.nbytes = }")
    print(f"{arr_flat.ndim = }")
    print(f"{arr_flat.size = }")
    print(f"{arr_flat.nbytes = }")
    m_vals = {
        1: [-1, 1],
        3: [-3, -1, 1, 3],
        5: [-5, -3, -1, 1, 3, 5],
    }
    j_vals = [1, 3, 5]
    j_max = max(j_vals) # For shifting the m-values.

    for j1 in j_vals:
        for m1 in m_vals[j1]:
            # j1_idx = (j1 + 1)//2 - 1    # [1, 3, 5] -> [0, 1, 2]
            # m1_idx = (m1 + j_max)//2    # [-5, -3, -1, 1, 3, 5] -> [0, 1, 2, 3, 4, 5], [-1, 1] -> [2, 3]
            # arr[j1_idx, m1_idx] = (j1, m1)

            for j2 in j_vals:
                for m2 in m_vals[j2]:
                    m3 = m1 + m2

                    for j3 in range( abs( j1 - j2 ), j1 + j2 + 2, 2 ):

                        j1_idx = (j1 + 1)//2 - 1    # [1, 3, 5] -> [0, 1, 2]
                        m1_idx = (m1 + j_max)//2    # [-5, -3, -1, 1, 3, 5] -> [0, 1, 2, 3, 4, 5], [-1, 1] -> [2, 3]
                        j2_idx = (j2 + 1)//2 - 1
                        m2_idx = (m2 + j_max)//2
                        j3_idx = j3//2              # [0,  2,  4,  6,  8, 10] -> [0, 1, 2, 3, 4, 5], [2, 4, 6, 8] -> [1, 2, 3, 4]

                        cg = CG(
                            j1 = j1/2,
                            m1 = m1/2,
                            j2 = j2/2,
                            m2 = m2/2,
                            j3 = j3/2,
                            m3 = m3/2,
                        )

                        calculated_flat_idx = ((((j1_idx*6 + m1_idx)*3 + j2_idx)*6 + m2_idx)*6 + j3_idx)    # Indexing a 1D array as if it was a 5D array.
                        arr[j1_idx, m1_idx, j2_idx, m2_idx, j3_idx] = float(cg.doit())
                        arr_flat[calculated_flat_idx] = float(cg.doit())

                        # if (j1_idx == 1) and (m1_idx == 0) and (j2_idx == 2) and (m2_idx == 4) and (j3_idx == 2):
                        #     print(f"{j1_idx = }")
                        #     print(f"{m1_idx = }")
                        #     print(f"{j2_idx = }")
                        #     print(f"{m2_idx = }")
                        #     print(f"{j3_idx = }")
                        #     print(f"{calculated_flat_idx = }")
                        #     print(f"{arr_flat[calculated_flat_idx] = }")
                        #     print(f"{j1 = }")
                        #     print(f"{m1 = }")
                        #     print(f"{j2 = }")
                        #     print(f"{m2 = }")
                        #     print(f"{j3 = }")

    msg = ", "
    with open("cg_flat_arr.txt", "w") as outfile:
        for i in range(len(arr_flat)):

            if (i + 1)%3 == 0:
                msg = ",\n"
            else:
                msg = ", "

            outfile.write(f"{arr_flat[i]:21.16f}" + msg)

    print(len(arr_flat[arr_flat == 0]))

if __name__ == "__main__":
    main()