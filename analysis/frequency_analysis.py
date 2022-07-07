import os
import numpy as np
from mpi4py import MPI

from read_data import ReadClass

experiment_name = "exRo9N10e01_LES"

figure_directory = "../figure/" + experiment_name
if not os.path.isdir(figure_directory):
    os.makedirs(figure_directory)

reader = ReadClass(experiment_name)
reader.read_geometry(0)
reader.read_budget(0)
time_days = reader.time / 86400
NT = np.size(time_days)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

thin = 32
n1_thinned = reader.n1 // thin
n3_thinned = reader.n3 // thin

freq_vortex = np.sqrt(reader.alpha_plus * reader.alpha_minus)
freq_vortex_days = freq_vortex * 86400
freq_buoyancy = reader.buoyancy_frequency * 86400
freq_inertia = (
    reader.angular_velocity[2]*2 - reader.alpha_minus - reader.alpha_plus
) * 86400
period_vortex = 2 * np.pi / freq_vortex

output_per_day = 86400 / reader.output_interval


def derive_frequency_spectrum(start_day, end_day, var):
    nt_start = int(start_day * output_per_day)
    nt_end = int(end_day * output_per_day)
    nTime = nt_end - nt_start
    delta_freq = 2 * np.pi / (reader.output_interval * nTime) * 86400
    freq_axis = (
        np.linspace(0, nTime/2, nTime//2+1, endpoint=True) * delta_freq)

    window = np.zeros(nTime)
    for it in range(0, nTime):
        a = 0.05
        if it < nTime // 2:
            jt = it
        else:
            jt = nTime - it
        if jt < a * nTime / 2:
            window[it] = 0.5 * (1 - np.cos(2 * np.pi * jt / nTime / a))
        else:
            window[it] = 1
    F_time = np.zeros((n1_thinned, n3_thinned, nTime))
    for it in range(nt_start, nt_end):
        F = reader.get_snapshot(var, rank, it, surface='XZ')
        jt = it - nt_start
        F_time[:, :, jt] = F[::thin, ::thin]
        if rank == 0:
            print(it, var)
    F_freq = np.fft.rfft(F_time * window, axis=2)
    spec = np.mean(np.abs(F_freq)**2, axis=(0, 1)) / nTime**2
    spec[1:-1] = spec[1:-1] * 2
    return freq_axis, spec


def derive_UV_frequency_spectrum(start_day, end_day):
    nt_start = int(start_day * output_per_day)
    nt_end = int(end_day * output_per_day)
    nTime = nt_end - nt_start
    delta_freq = 2 * np.pi / (reader.output_interval * nTime) * 86400
    freq_axis = (np.linspace(0, nTime/2, nTime//2+1, endpoint=True)
                 * delta_freq)

    window = np.zeros(nTime)
    for it in range(0, nTime):
        a = 0.05
        if it < nTime / 2:
            jt = it
        else:
            jt = nTime - it
        if jt < a * nTime / 2:
            window[it] = 0.5 * (1 - np.cos(2 * np.pi * jt / nTime / a))
        else:
            window[it] = 1
    U_time = np.zeros((n1_thinned, n3_thinned, nTime))
    V_time = np.zeros((n1_thinned, n3_thinned, nTime))
    for it in range(nt_start, nt_end):
        time = it * reader.output_interval
        jt = it - nt_start
        U = reader.get_snapshot('U', rank, it, surface='XZ')
        V = reader.get_snapshot('V', rank, it, surface='XZ')
        U_time[:, :, jt] = (U[::thin, ::thin] * np.cos(freq_vortex * time)
                            - V[::thin, ::thin] * np.sin(freq_vortex * time))
        V_time[:, :, jt] = (U[::thin, ::thin] * np.sin(freq_vortex * time)
                            + V[::thin, ::thin] * np.cos(freq_vortex * time))
        if rank == 0:
            print(it)
    U_freq = np.fft.rfft(U_time * window, axis=2)
    V_freq = np.fft.rfft(U_time * window, axis=2)
    spec_U = np.mean(np.abs(U_freq)**2, axis=(0, 1)) / nTime**2
    spec_V = np.mean(np.abs(V_freq)**2, axis=(0, 1)) / nTime**2
    spec_U[1:-1] = spec_U[1:-1] * 2
    spec_V[1:-1] = spec_V[1:-1] * 2
    return freq_axis, spec_U, spec_V


def derive_KE_spectrum(start_day, end_day):
    freq_axis, U_spec, V_spec = derive_UV_frequency_spectrum(
        start_day, end_day)
    freq_axis, W_spec = derive_frequency_spectrum(
        start_day, end_day, 'W')
    KE_spec = (U_spec + V_spec + W_spec) * 0.5
    KE_spec_global = np.zeros_like(KE_spec)
    comm.Reduce(KE_spec, KE_spec_global, op=MPI.SUM, root=0)
    KE_spec_global /= size
    return freq_axis, KE_spec_global


def derive_PE_spectrum(start_day, end_day):
    freq_axis, T_spec = derive_frequency_spectrum(
        start_day, end_day, 'T')
    PE_spec = T_spec * 0.5
    PE_spec_global = np.zeros_like(PE_spec)
    comm.Reduce(PE_spec, PE_spec_global, op=MPI.SUM, root=0)
    PE_spec_global /= size
    return freq_axis, PE_spec_global


freq_axis, KE_spec = derive_KE_spectrum(40, 80)
if (rank == 0):
    np.save('./data/KE_spec_series_1', KE_spec)
    np.save('./data/freq_axis_1', freq_axis)

freq_axis, KE_spec = derive_KE_spectrum(80, 125)
if (rank == 0):
    np.save('./data/KE_spec_series_2', KE_spec)
    np.save('./data/freq_axis_2', freq_axis)

freq_axis, KE_spec = derive_KE_spectrum(160, 220)
if (rank == 0):
    np.save('./data/KE_spec_series_3', KE_spec)
    np.save('./data/freq_axis_3', freq_axis)

freq_axis, PE_spec = derive_PE_spectrum(40, 80)
if (rank == 0):
    np.save('./data/PE_spec_series_1', PE_spec)

freq_axis, PE_spec = derive_PE_spectrum(80, 125)
if (rank == 0):
    np.save('./data/PE_spec_series_2', PE_spec)

freq_axis, PE_spec = derive_PE_spectrum(160, 220)
if (rank == 0):
    np.save('./data/PE_spec_series_3', PE_spec)
