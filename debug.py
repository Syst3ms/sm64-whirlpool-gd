import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    before_x = []
    before_z = []
    before_theta = []
    after_x = []
    after_z = []
    after_theta = []

    with open("debug.txt", "r") as f:
        lines = f.readlines()
        for l in lines[:128]:
            coords = tuple(map(float, l.rstrip().split(",")))
            before_x.append(coords[0])
            before_z.append(coords[1])
            before_theta.append(coords[2])
        for l in lines[129:-1]:
            coords = tuple(map(float, l.rstrip().split(",")))
            after_x.append(coords[0])
            after_z.append(coords[1])
            after_theta.append(coords[2])

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

    plt.plot(before_x, before_z, "r-")
    plt.plot(after_x, after_z, "b-")

    for i, t in enumerate(before_theta):
        if i%3:
            continue
        x = before_x[i]
        z = before_z[i]
        dxp, dzp = np.sin(t) * 28, np.cos(t) * 28

        plt.arrow(
            x, z,
            dxp * 5, dzp * 5,
            width = 4, zorder=100,
            color = (1.0, 1.0, 0.0)
        )
    for i, t in enumerate(after_theta):
        if i%3:
            continue
        x = after_x[i]
        z = after_z[i]
        dxp, dzp = np.sin(t) * 28, np.cos(t) * 28

        plt.arrow(
            x, z,
            dxp * 5, dzp * 5,
            width = 4, zorder=100,
            color = (1.0, 1.0, 0.0)
        )

    ax.axis('equal')
    plt.show()