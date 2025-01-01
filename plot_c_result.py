import matplotlib.pyplot as plt
import numpy as np

def au_to_radians(a):
    return a / 65536 * 2 * np.pi

def whirl(x, z):
    norm = np.sqrt(x ** 2 + z ** 2)
    yaw_offset = - np.pi * 250 / (norm + 1000)
    sin_o = np.sin(yaw_offset)
    cos_o = np.cos(yaw_offset)
    fac = -20 * (1/norm - 1/2000)
    whirl_x = fac * (cos_o * x + sin_o * z)
    whirl_z = fac * (-sin_o * x + cos_o * z)
    return whirl_x, whirl_z

if __name__ == '__main__':
    print("Reading paths")
    memory = []
    Mx = []
    Mz = []
    with open("path.txt", "r") as f:
        for ln in f:
            if ln != "\n":
                coords = ln.split(",")
                Mx.append(float(coords[0]))
                Mz.append(float(coords[1]))
            else:
                memory.append((Mx, Mz))
                Mx = []
                Mz = []
    
    
    print("Reading debug")
    yaw_path = []
    with open("resampled.txt", "r") as f:
        f.readline()
        for ln in f:
            vals = ln.split(",")
            yaw_path.append((float(vals[0]), float(vals[1]), au_to_radians(float(vals[2]))))

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
        if i%3:
            continue
        # angles oriented relative to z+
        dxp, dzp = np.sin(t) * 28, np.cos(t) * 28
        dxw, dzw = whirl(x,z)
        plt.arrow(
            x, z,
            dxw, dzw,
            width = 10, zorder=100,
            color = (0.0, 1.0, 0.0)
        )
        plt.arrow(
            x, z,
            dxp * 5, dzp * 5,
            width = 10, zorder=100,
            color = (1.0, 1.0, 0.0)
        )
        plt.arrow(
            x, z,
            dxp + dxw, dzp + dzw,
            width = 10, zorder=100,
            color = "k"
        )
    
    # plt.plot(Mx, Mz, 'ro')
    ax.add_patch(plt.Circle((0,0), 26))
    A = np.linspace(0, 2 * np.pi, 100)
    plt.plot(2000 * np.cos(A), 2000 * np.sin(A), 'k-')
    ax.axis('equal')
    plt.show()
