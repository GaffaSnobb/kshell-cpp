import numpy as np

def main():
    dim = 32
    a = np.random.normal(size=(dim, dim))
    a += a.T
    # print(a.ravel())

    for i, e in enumerate(a.ravel()):
        if i%4 == 0: print()
        
        print(f"{e}, ", end="")

if __name__ == "__main__":
    main()