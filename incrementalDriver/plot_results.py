import numpy as np
import matplotlib.pylab as plt


def principal_stresses(stress_tensor):
    """
    Computes the principal stresses from a given 3D stress tensor.

    :param stress_tensor: 3x3 symmetric stress tensor.

    :return: Principal stresses sorted in descending order.
    """
    # Compute eigenvalues of the stress tensor (principal stresses)
    eigenvalues = np.linalg.eigvals(stress_tensor)

    # Sort in descending order (by convention: sigma_1 >= sigma_2 >= sigma_3)
    principal_stresses = np.sort(eigenvalues)[::-1]

    return principal_stresses


def read_file(file_name, cohesion=0, friction_angle=30, psi=0):
    """
    Reads the output file and computes the principal stresses, q, p' and the Mohr-Coulomb yield criterion.

    :param file_name: Name of the output file.
    :param cohesion: Cohesion of the material.
    :param friction_angle: Friction angle of the material.
    :param psi: Dilatancy angle.

    :return: Dictionary containing the computed values.
    """

    with open(file_name, "r") as fo:
        data = fo.read().splitlines()

    friction_angle = friction_angle * np.pi / 180
    psi = psi * np.pi / 180

    data = [list(map(float, i.split())) for i in data[1:]]
    data = np.array(data)


    sigma_1, sigma_2, sigma_3 = np.zeros(len(data)), np.zeros(len(data)), np.zeros(len(data))
    for i in range(len(data)):
        sigma_xx = -data[i][8]
        sigma_yy = -data[i][9]
        sigma_zz = -data[i][10]
        sigma_xy = data[i][11]
        sigma_xz = data[i][12]
        sigma_yz = data[i][13]

        stress_tensor = np.array([[sigma_xx, sigma_xy, sigma_xz],
                                [sigma_xy, sigma_yy, sigma_yz],
                                [sigma_xz, sigma_yz, sigma_zz]])

        # Compute principal stresses
        s1, s2, s3 = principal_stresses(stress_tensor)
        sigma_1[i] = s1
        sigma_2[i] = s2
        sigma_3[i] = s3


    p = (sigma_1 + sigma_2 + sigma_3) / 3
    I2 = sigma_2 * sigma_3 + sigma_3 * sigma_1 + sigma_1 * sigma_2

    J2 = (3 * p) ** 2 / 3 - I2
    q = (3 * J2)

    mc_yield = (sigma_1 - sigma_3)  + (sigma_1 + sigma_3) * np.sin(friction_angle) - 2 * cohesion * np.cos(friction_angle)

    results= {"p": p,
              "q": q,
              "mc_yield": mc_yield,
              "sigma_xx": sigma_xx,
              "sigma_yy": sigma_yy,
              "sigma_zz": sigma_zz,
              "sigma_xy": sigma_xy,
              "sigma_xz": sigma_xz,
              "sigma_yz": sigma_yz,
              "epsilon_xx": data[:, 2],
              "epsilon_yy": data[:, 3],
              "epsilon_zz": data[:, 4],
              "epsilon_xy": data[:, 5],
              "epsilon_xz": data[:, 6],
              "epsilon_yz": data[:, 7],
              "sigma_1": sigma_1,
              "sigma_2": sigma_2,
              "sigma_3": sigma_3
              }

    return results

def make_plots(data):

    fig, ax = plt.subplots(3, 1, figsize=(8, 10))

    ax[0].plot(data["p"], data["q"])
    ax[0].plot(data["p"][0], data["q"][0], "ro")
    ax[0].plot(data["p"][-1], data["q"][-1], "bo")
    ax[0].plot(data["p"], data["mc_yield"], "k--")
    ax[0].set_xlabel("p'")
    ax[0].set_ylabel("q")
    ax[0].grid()

    ax[1].plot(data["epsilon_yy"] , data["q"])
    ax[1].plot(data["epsilon_yy"][0], data["q"][0], "ro")
    ax[1].plot(data["epsilon_yy"][-1], data["q"][-1], "bo")

    ax[1].set_xlabel("Axial strain")
    ax[1].set_ylabel("q")
    ax[1].grid()

    volumetric_strain = data["epsilon_xx"] + data["epsilon_yy"] + data["epsilon_zz"]
    ax[2].plot(data["epsilon_yy"] , volumetric_strain)
    ax[2].plot(data["epsilon_yy"][0], volumetric_strain[0], "ro")
    ax[2].plot(data["epsilon_yy"][-1], volumetric_strain[-1], "bo")

    ax[2].set_xlabel("Axial strain")
    ax[2].set_ylabel("Volumetric strain")
    ax[2].grid()

    plt.show()


if __name__ == "__main__":
    results = read_file("Txt_UnComAni.out", 0, 30, 0)
    make_plots(results)
