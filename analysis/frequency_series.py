import numpy as np
from mpi4py import MPI

from read_data import ReadClass

experiment_name = "exRo9N10e01_LES"

reader = ReadClass(experiment_name)
reader.read_geometry(0)
reader.read_budget(0)
time_days = reader.time / 86400
NT = np.size(time_days)

num = round(reader.vortex_period / reader.output_interval)
b = np.ones(num) / num

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

thin = 32
n1_thinned = reader.n1 // thin
n3_thinned = reader.n3 // thin

freq_vortex = np.sqrt(reader.alpha_plus * reader.alpha_minus)
period_vortex = 2 * np.pi / freq_vortex
MT = round(period_vortex / reader.output_interval)

nt_start = MT * 20
nt_end = MT * 150

nt_interval = nt_start
nTime = nt_end - nt_start

time_axis = np.arange(nt_start, nt_end) - nt_interval // 2
time_axis = time_axis * reader.output_interval / 86400

delta_freq = 2 * np.pi / (reader.output_interval * nt_interval)
freq_axis = np.arange(0, nt_interval // 2 + 1)
freq_axis = freq_axis * delta_freq * 86400

window = np.zeros(nt_interval)
for it in range(0, nt_interval):
    a = 0.05
    if it < nt_interval // 2:
        jt = it
    else:
        jt = nt_interval - it
    if jt < a * nt_interval / 2:
        window[it] = 0.5 * (1 - np.cos(2 * np.pi * jt / nt_interval / a))
    else:
        window[it] = 1


def derive_spec(var):
    spec_series = np.zeros((nt_interval // 2 + 1, nTime))
    F_time = np.zeros((n1_thinned, n3_thinned, nt_interval))
    for it in range(nt_start - nt_interval, nt_start):
        jt = it - (nt_start - nt_interval)
        F = reader.get_snapshot(var, rank, it, surface='XZ')
        F_time[:, :, jt] = F[::thin, ::thin]
        if rank == 0:
            print(var, it)

    for it in range(nt_start, nt_end):
        jt = it - nt_start
        F = reader.get_snapshot(var, rank, it, surface='XZ')
        F_time[:, :, :-1] = F_time[:, :, 1:]
        F_time[:, :, -1] = F[::thin, ::thin]
        F_freq = np.fft.rfft(F_time * window, axis=2)
        spec_series[:, jt] = (
            np.mean(np.abs(F_freq)**2, axis=(0, 1)) / nt_interval**2)
        if rank == 0:
            print(var, it)
    spec_series[1:-1, :] = spec_series[1:-1, :] * 2
    return spec_series


def derive_U_spec():
    spec_series = np.zeros((int(nt_interval / 2) + 1, nTime))
    U_time = np.zeros((n1_thinned, n3_thinned, nt_interval))
    for it in range(nt_start - nt_interval, nt_start):
        time = it * reader.output_interval
        jt = it - (nt_start - nt_interval)
        U = reader.get_snapshot('U', rank, it, surface='XZ')
        V = reader.get_snapshot('V', rank, it, surface='XZ')
        U_time[:, :, jt] = (U[::thin, ::thin] * np.cos(freq_vortex * time)
                            - V[::thin, ::thin] * np.sin(freq_vortex * time))
        if rank == 0:
            print('U', it)

    for it in range(nt_start, nt_end):
        time = it * reader.output_interval
        U_time[:, :, :-1] = U_time[:, :, 1:]
        jt = it - nt_start
        U = reader.get_snapshot('U', rank, it, surface='XZ')
        V = reader.get_snapshot('V', rank, it, surface='XZ')

        U_time[:, :, -1] = (U[::thin, ::thin] * np.cos(freq_vortex * time)
                            - V[::thin, ::thin] * np.sin(freq_vortex * time))
        U_freq = np.fft.rfft(U_time * window, axis=2)

        spec_series[:, jt] = (
            np.sum(np.abs(U_freq)**2, axis=(0, 1))
            / (n1_thinned * n3_thinned * nt_interval**2))
        if rank == 0:
            print('U', it)
    spec_series[1:-1, :] = spec_series[1:-1, :] * 2
    return spec_series


def derive_V_spec():
    spec_series = np.zeros((int(nt_interval / 2) + 1, nTime))
    V_time = np.zeros((n1_thinned, n3_thinned, nt_interval))
    for it in range(nt_start - nt_interval, nt_start):
        time = it * reader.output_interval
        jt = it - (nt_start - nt_interval)
        U = reader.get_snapshot('U', rank, it, surface='XZ')
        V = reader.get_snapshot('V', rank, it, surface='XZ')
        V_time[:, :, jt] = (U[::thin, ::thin] * np.sin(freq_vortex * time)
                            + V[::thin, ::thin] * np.cos(freq_vortex * time))
        if rank == 0:
            print('V', it)

    for it in range(nt_start, nt_end):
        time = it * reader.output_interval
        V_time[:, :, :-1] = V_time[:, :, 1:]
        jt = it - nt_start
        U = reader.get_snapshot('U', rank, it, surface='XZ')
        V = reader.get_snapshot('V', rank, it, surface='XZ')

        V_time[:, :, -1] = (U[::thin, ::thin] * np.sin(freq_vortex * time)
                            + V[::thin, ::thin] * np.cos(freq_vortex * time))
        V_freq = np.fft.rfft(V_time * window, axis=2)

        spec_series[:, jt] = (
            + np.sum(np.abs(V_freq)**2, axis=(0, 1))
            / (n1_thinned * n3_thinned * nt_interval**2))
        if rank == 0:
            print('V', it)
    spec_series[1:-1, :] = spec_series[1:-1, :] * 2
    return spec_series


U_spec_series = derive_U_spec() * 0.5
V_spec_series = derive_V_spec() * 0.5
W_spec_series = derive_spec(var='W') * 0.5
KE_spec_series = U_spec_series + V_spec_series + W_spec_series
PE_spec_series = derive_spec(var='T') * 0.5

KE_spec_series_global = np.zeros_like(KE_spec_series)
PE_spec_series_global = np.zeros_like(PE_spec_series)
comm.Reduce(KE_spec_series, KE_spec_series_global, op=MPI.SUM, root=0)
comm.Reduce(PE_spec_series, PE_spec_series_global, op=MPI.SUM, root=0)
KE_spec_series_global /= size
PE_spec_series_global /= size

if rank == 0:
    np.save('./data/KE_spec_series_full', KE_spec_series)
    np.save('./data/PE_spec_series_full', PE_spec_series)
