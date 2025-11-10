import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
import json

from plot_openmx_band import get_band_from_file


Ha2eV = 27.21138602


def read_file(openmx_file="openmx.out"):
    with open(openmx_file, "r") as out_file:
        lines = out_file.readlines()
    for line_index, line in enumerate(lines):
        if line.find("Chemical Potential (Hatree)") != -1:
            break
    fermi_energy = float(line.split()[-1])
    line_index += 5
    line_index_tmp = line_index + 3
    while len(lines[line_index_tmp].split()) == 2:
        num_band = int(lines[line_index_tmp].split()[0])
        line_index_tmp += 1
    print("[INFO] number of bands: %d" % num_band)

    k_vector_list = []
    eigen_values_list = []
    while lines[line_index].find("kloop=") != -1:
        num_kpoint = int(lines[line_index][lines[line_index].find("kloop=") + 6 :])
        k_vector = np.array(
            [
                float(lines[line_index + 1].split()[1]),
                float(lines[line_index + 1].split()[3]),
                float(lines[line_index + 1].split()[5]),
            ]
        )
        k_vector_list.append(k_vector)
        eigen_values = np.zeros(num_band)
        for band_index in range(num_band):
            band_index_file, eigen_value = lines[line_index + 3 + band_index].split()
            eigen_values[band_index] = eigen_value
            assert band_index == int(band_index_file) - 1
        eigen_values_list.append(eigen_values)
        line_index += num_band + 4
    k_vectors = np.stack(k_vector_list)  # shape: (num_kpoint, 3)
    eigen_values = np.stack(eigen_values_list)  # shape: (num_kpoint, num_band)
    num_kpoint += 1
    print("[INFO] number of k points: %d" % num_kpoint)
    return k_vectors, (eigen_values - fermi_energy) * Ha2eV


def read_file_kpath(openmx_file="openmx.Band"):
    with open(openmx_file, "r") as out_file:
        lines = out_file.readlines()
    fermi_energy = float(lines[0].split()[-1])
    num_band = len(lines[13].split())
    print("[INFO] number of bands: %d" % num_band)
    k_vector_list = []
    eigen_values_list = []
    for line in lines[12:]:
        if len(line.split()) == 4:
            assert int(line.split()[0]) == num_band
            k_vector = np.array(
                [float(line.split()[1]), float(line.split()[2]), float(line.split()[3])]
            )
            k_vector_list.append(k_vector)
        else:
            assert len(line.split()) == num_band
            eigen_values = np.array(
                [float(eigen_value) for eigen_value in line.split()]
            )
            eigen_values_list.append(eigen_values)
    k_vectors = np.stack(k_vector_list)  # shape: (num_kpoint, 3)
    eigen_values = np.stack(eigen_values_list)  # shape: (num_kpoint, num_band)
    print("[INFO] number of k points: %d" % k_vectors.shape[0])

    return k_vectors, (eigen_values - fermi_energy) * Ha2eV


class DensityOfStates:
    def __init__(self, k_vectors, eigen_values, delta):
        self.k_vectors = k_vectors
        # self.eigen_values = eigen_values[:, 20:50]
        self.eigen_values = eigen_values
        self.delta = delta

    def dos_omega(self, omega):
        dos = self.delta(self.eigen_values[..., None] - omega[None, :])
        dos = dos.sum(0)
        return dos.sum(0)

    def jdos_omega(self, omega):
        num_kpoint = self.eigen_values.shape[0]
        num_band = self.eigen_values.shape[1]
        num_omega = omega.shape[0]
        f_n = (self.eigen_values > 0).astype(int)
        jdos = np.zeros([num_kpoint, num_omega])
        for n in range(num_band):
            for m in range(num_band):
                f_nm = f_n[:, n] - f_n[:, m]
                if f_nm.sum() != 0:
                    delta_nm = self.delta(
                        self.eigen_values[:, n, None]
                        - self.eigen_values[:, m, None]
                        - omega[None, :]
                    )
                    jdos += delta_nm * f_nm[:, None]
        jdos = jdos.sum(0)
        return jdos
        # jdos = jdos.sum(0)
        # plt.figure(figsize=[5, 10])
        # plt.plot(ret, omega)
        # plt.show()


def analysis_omega_band(eigen_values, omega_slider, omega_tolerance=0.05):
    color_map_bands = np.zeros_like(eigen_values)
    eigen_values = eigen_values.transpose()
    num_band = eigen_values.shape[1]
    f_n = (eigen_values > 0).astype(int)
    meet_list = []
    for n in range(num_band):
        for m in range(num_band):
            f_nm = f_n[:, n] - f_n[:, m]
            if f_nm.sum() > 0:
                omega_nm = eigen_values[:, n] - eigen_values[:, m] - omega_slider
                epsilon = 0.3
                color_map_bands[n] += (
                    np.exp(-omega_nm * omega_nm / epsilon / epsilon)
                    / np.sqrt(np.pi)
                    / epsilon
                    * 2
                )
                color_map_bands[m] += (
                    np.exp(-omega_nm * omega_nm / epsilon / epsilon)
                    / np.sqrt(np.pi)
                    / epsilon
                )
    return color_map_bands


def plot_band(
    band_axes,
    omega_slider,
    plot_hsk_symbol_list,
    hsk_corrdinate_list,
    kpoints_corrdinate,
    spin_up_bands_data,
):
    min_plot_energy = -1.5
    max_plot_energy = 2.5
    band_axes.xaxis.set_tick_params(labelsize=10)
    band_axes.yaxis.set_tick_params(labelsize=8)
    band_axes.set_xlim(0.0, hsk_corrdinate_list[-1])
    band_axes.set_ylim(min_plot_energy, max_plot_energy)
    band_axes.set_xlabel("", fontsize=10)
    band_axes.set_ylabel("Energy (eV)", fontsize=10)
    band_axes.set_xticks(hsk_corrdinate_list)
    band_axes.set_xticklabels(plot_hsk_symbol_list)
    for hsk_corrdinate in hsk_corrdinate_list:
        band_axes.vlines(
            hsk_corrdinate,
            min_plot_energy,
            max_plot_energy,
            colors="grey",
            linewidth=0.6,
        )
    # Plot the fermi energy surface with a dashed line
    band_axes.hlines(
        0.0,
        0.0,
        hsk_corrdinate_list[-1],
        colors="grey",
        linestyles="dashed",
        linewidth=0.6,
    )
    # Plot the Band Structure
    color_map_bands = analysis_omega_band(spin_up_bands_data, omega_slider)
    lc_list = []
    for index_band, spin_up_band_data in enumerate(spin_up_bands_data):
        points = np.array([kpoints_corrdinate, spin_up_band_data]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        color_map = color_map_bands[index_band]

        norm = plt.Normalize(color_map_bands.min(), color_map_bands.max())
        lc = LineCollection(segments, cmap="coolwarm", norm=norm)
        lc.set_array(color_map)
        lc.set_linewidth(1)
        lc.set_zorder(2)
        band_axes.add_collection(lc)
        lc_list.append(lc)

    # Add colorbar
    if lc_list:
        cbar = plt.colorbar(
            lc_list[0],
            ax=band_axes,
            pad=0.02,
            aspect=30,
        )
        # Set only two ticks: min and max
        vmin = color_map_bands.min()
        vmax = color_map_bands.max()
        cbar.set_ticks([vmin, vmax])
        cbar.set_ticklabels(["0", "+"])
        cbar.ax.tick_params(labelsize=8)

    # os.system("cls")
    # for k_indices, n, m in meet_list:
    #     for k_index in k_indices:
    #         print(n, m, str(k_index) + '/%d' % len(kpoints_corrdinate))
    #     x = np.concatenate([kpoints_corrdinate[k_indices]] * 2)
    #     y = np.concatenate([spin_up_bands_data[n, k_indices].reshape(-1), spin_up_bands_data[m, k_indices].reshape(-1)])
    #     # print(spin_up_bands_data[n, k_indices] - spin_up_bands_data[m, k_indices])
    #     band_axes.scatter(x, y, alpha=1, s=5, color='red')


if __name__ == "__main__":
    plt.rcParams.update(
        {"font.size": 10, "font.family": "Times new roman", "mathtext.fontset": "cm"}
    )
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    fig = plt.figure(figsize=(5.5 / 2.54, 4.8 / 2.54))
    band_axes = fig.add_subplot(111)

    omega_slider = 1.1923076923076923  # 1.1923076923076923   480.4580165287854

    # wosoc
    (
        plot_hsk_symbol_list,
        hsk_corrdinate_list,
        kpoints_corrdinate_list,
        spin_up_bands_data,
    ) = get_band_from_file("wosoc/openmx.Band")
    kpoints_corrdinate = np.array(kpoints_corrdinate_list)
    plot_band(
        band_axes,
        omega_slider,
        plot_hsk_symbol_list,
        hsk_corrdinate_list,
        kpoints_corrdinate,
        spin_up_bands_data,
    )

    # # soc
    # (
    #     plot_hsk_symbol_list,
    #     hsk_corrdinate_list,
    #     kpoints_corrdinate_list,
    #     spin_up_bands_data,
    # ) = get_band_from_file("soc/openmx.Band")
    # kpoints_corrdinate = np.array(kpoints_corrdinate_list)
    # for index_band, spin_up_band_data in enumerate(spin_up_bands_data):
    #     band_axes.plot(
    #         kpoints_corrdinate,
    #         spin_up_band_data,
    #         "k",
    #         linewidth=0.5,
    #         linestyle="--",
    #         zorder=100,
    #     )

    plt.tight_layout()
    plt.savefig("band.png", bbox_inches="tight", dpi=1000)
    plt.savefig("band.svg", bbox_inches="tight", transparent=True)
