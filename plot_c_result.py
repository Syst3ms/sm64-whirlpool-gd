import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    print("Reading paths")
    memory = []
    Mx = []
    Mz = []
    yaws = []
    with open("path.txt", "r") as f:
        for ln in f:
            if ln != "\n":
                coords = ln.split(",")
                Mx.append(float(coords[0]))
                Mz.append(float(coords[1]))
                if len(coords) >= 3:
                    yaws.append(float(coords[2]))
            else:
                memory.append((Mx, Mz))
                Mx = []
                Mz = []
    
    print("Reading yaws")
    yaw_path = []
    with open("yaws.txt", "r") as f:
        for ln in f:
            vals = ln.split(",")
            yaw_path.append((
                float(vals[0]),
                float(vals[1]), 
                float(vals[2])
            ))

    print("Done reading from files, plotting")

    N = len(memory)
    for i, (oldMx, oldMz) in enumerate(memory):
        plt.plot(oldMx, oldMz, '-', color=(i / N, 0.0, 1.0 - i / N))

    ax = plt.gca()

    ax.set_xlabel("X")
    ax.set_ylabel("Z")
    ax.invert_yaxis()

    rx = np.array([-1326, 1374, -1326, 774, 1274])
    rz = np.array([1198, 948, -1202, 23, -1502])
    txt = [f"Chest {i}" for i in range(4)] + ["Star"]

    plt.plot(rx, rz, 'mo')
    for i in range(len(rx)):
        if i < 4:
            ax.add_patch(plt.Circle((rx[i], rz[i]), 150, fill=False))
        ax.annotate(txt[i], (rx[i], rz[i]))

    plt.plot(Mx, Mz, 'r-')
    n = len(yaw_path)
    for i, (x, z, t) in enumerate(yaw_path):
        plt.arrow(
            x, z,
            np.sin(t) * 300, np.cos(t) * 300,
            width = 10, zorder=100,
            color = (i/n, 1.0, 0.0)
        )
    
    # plt.plot(Mx, Mz, 'ro')
    ax.add_patch(plt.Circle((0,0), 26))
    A = np.linspace(0, 2 * np.pi, 100)
    plt.plot(2000 * np.cos(A), 2000 * np.sin(A), 'k-')
    ax.axis('equal')
    plt.show()
