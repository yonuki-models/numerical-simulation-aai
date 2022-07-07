from os import path

import f90nml
import numpy as np

FACTOR_SINGLE_PRECISION = 4


class ReadClass:
    def __init__(self, experiment_name, base_directory='./'):
        self.experiment_name = experiment_name
        self.data_directory = (
            base_directory + '../data/' + experiment_name + '/')
        namelist = f90nml.read(
            self.data_directory + experiment_name + '.out')
        tmp = namelist['PARAMS']

        self.n = tmp['NUMBER_OF_GRID']
        self.division = tmp['DIVISION']
        self.length = tmp['DOMAIN_LENGTH']
        self.time_end = tmp['TIME_END']
        self.output_interval = tmp['OUTPUT_INTERVAL']
        self.output_state_interval = tmp['OUTPUT_STATE_INTERVAL']
        self.buoyancy_frequency = tmp['BUOYANCY_FREQUENCY']
        self.angular_velocity = tmp['ANGULAR_VELOCITY']
        self.alpha_plus = tmp['ALPHA_PLUS']
        self.alpha_minus = tmp['ALPHA_MINUS']
        self.N_Ri = tmp['N_RI']
        self.Ri_min = tmp['RI_MIN']
        self.Ri_max = tmp['RI_MAX']
        self.N_Thorpe = tmp['N_THORPE']

        self.n1 = self.n[0]
        self.n2 = self.n[1]
        self.n3 = self.n[2]

        self.n1_half = self.n1 // 2
        self.n2_half = self.n2 // 2
        self.n3_half = self.n3 // 2

        self.n1_truncated = self.n1 // 3
        self.n2_truncated = self.n2 // 3
        self.n3_truncated = self.n3 // 3

        self.l1 = self.length[0]
        self.l2 = self.length[1]
        self.l3 = self.length[2]

        self.n2C_local = self.n2 // self.division[1]
        self.n1R_local = self.n1 // self.division[0]
        self.r = np.sqrt(self.alpha_plus / self.alpha_minus)
        self.e = 1 - 1.0 / self.r
        self.Ro = (self.alpha_plus + self.alpha_minus) / (
            2 * self.angular_velocity[2])

        self.f = 2 * self.angular_velocity[2]
        self.vortex_frequency = np.sqrt(self.alpha_plus * self.alpha_minus)
        self.vortex_period = 2 * np.pi / self.vortex_frequency

    def _read_data(self, file_name, data_type, shape, record=0):
        size = shape.prod()
        data_bytes = round(int(data_type.name[-2:]) / 8)
        offset = size * data_bytes * record
        file_open = open(self.data_directory + file_name, 'r')
        data_tmp = np.fromfile(file_open, dtype=data_type,
                               count=size, offset=offset)
        data = data_tmp.reshape(shape, order='F')
        return data

    def read_geometry(self, index_finish):
        file_name = 'geometry.out'
        file_size = path.getsize(self.data_directory + file_name)
        number_of_columns = 16
        number_of_record = int(
            file_size / (number_of_columns * FACTOR_SINGLE_PRECISION))
        # print(number_of_record, index_finish)
        if (number_of_record < index_finish):
            print('number_of_record should be greater than index_finish.')
        elif (index_finish != 0):
            number_of_record = index_finish
        shape = np.array([16, number_of_record])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape)
        self.time = data[0, :]
        self.A = np.zeros([number_of_record, 3, 3])
        self.A[:, 0, 0] = data[7, :]
        self.A[:, 0, 1] = data[8, :]
        self.A[:, 0, 2] = data[9, :]
        self.A[:, 1, 0] = data[10, :]
        self.A[:, 1, 1] = data[11, :]
        self.A[:, 1, 2] = data[12, :]
        self.A[:, 2, 0] = data[13, :]
        self.A[:, 2, 1] = data[14, :]
        self.A[:, 2, 2] = data[15, :]
        self.number_of_record = number_of_record

    def read_budget(self, index_finish):
        file_name = 'budget.out'
        file_size = path.getsize(self.data_directory + file_name)
        number_of_columns = 22
        number_of_record = int(
            file_size / (number_of_columns * FACTOR_SINGLE_PRECISION))
        if (number_of_record < index_finish):
            print('number_of_record should be greater than index_finish.')
        elif (index_finish != 0):
            number_of_record = index_finish
        shape = np.array([number_of_columns, number_of_record])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape)
        self.time = data[0, :]
        self.KE = data[1, :]
        self.PE = data[2, :]
        self.E = data[3, :]
        self.WE = data[4, :]
        self.VE = data[5, :]
        self.PWE = data[6, :]
        self.TKE = data[7, :]
        self.TPE = data[8, :]
        self.PKE = data[9, :]
        self.PPE = data[10, :]
        self.PW1 = data[11, :]
        self.PW2 = data[12, :]
        self.PVE = data[13, :]
        self.CKP = data[14, :]
        self.DKE = data[15, :]
        self.DPE = data[16, :]
        self.EN2 = data[17, :]
        self.EN3 = data[18, :]
        self.EN4 = data[19, :]
        self.DCH = data[20, :]
        self.DCV = data[21, :]
        self.number_of_record = number_of_record

    def get_complex_data_mpi(self, it, variable_name, rank):
        n2 = self.n2
        n3 = self.n3
        nt1 = self.n1_truncated
        nt2 = self.n2_truncated
        nt3 = self.n3_truncated

        n1_half = self.n1_half
        n2C_local = self.n2C_local

        C = np.zeros([n1_half+1, n2C_local, n3],
                     dtype='complex64', order='F')
        if (rank * n2C_local > nt2 and (rank + 1) * n2C_local - 1 < n2 - nt2):
            return C
        shape = np.array([nt1+1, n2C_local, 2*nt3+1])
        data_type = np.dtype('<c8')
        str_ip = '{0:04d}'.format(rank)
        file_name = 'state/' + variable_name + str_ip + '.out'
        C_div = self._read_data(file_name, data_type, shape, it)

        C[0:nt1, 0:n2C_local, 0:nt3]  \
            = C_div[0:nt1, 0:n2C_local, 0:nt3]
        C[0:nt1, 0:n2C_local, n3-nt3:n3]  \
            = C_div[0:nt1, 0:n2C_local, nt3+1:2*nt3+1]
        return C

    def get_wavenumber_mpi(self, it, rank):
        n2 = self.n2
        n3 = self.n3

        n1_half = self.n1_half
        n2_half = self.n2_half
        n3_half = self.n3_half
        n2C_local = self.n2C_local

        delta_1 = 2 * np.pi / self.l1
        K1_local = np.zeros(n1_half+1)
        for i1 in range(n1_half+1):
            K1_local[i1] = delta_1 * i1

        delta_2 = 2 * np.pi / self.l2
        K2_local = np.zeros(n2C_local)
        for i2 in range(n2C_local):
            i2g = rank * n2C_local + i2
            if (i2g <= n2_half):
                K2_local[i2] = delta_2 * i2g
            else:
                K2_local[i2] = delta_2 * (i2g - n2)

        delta_3 = 2 * np.pi / self.l3
        K3_local = np.zeros(n3)
        for i3 in range(n3):
            if (i3 <= n3_half):
                K3_local[i3] = delta_3 * i3
            else:
                K3_local[i3] = delta_3 * (i3 - n3)

        K1l, K2l, K3l = (
            np.meshgrid(K1_local, K2_local, K3_local, indexing='ij'))

        K = np.zeros([3, n1_half+1, n2C_local, n3])
        A = self.A[it, :, :]
        K[0, :] = A[0, 0] * K1l + A[0, 1] * K2l + A[0, 2] * K3l
        K[1, :] = A[1, 0] * K1l + A[1, 1] * K2l + A[1, 2] * K3l
        K[2, :] = A[2, 0] * K1l + A[2, 1] * K2l + A[2, 2] * K3l

        return K

    def get_energy_horizontam_mean(self, it):
        NH_max = int(np.sqrt(self.r) * self.n1_truncated + 0.5)
        NV_max = self.n3_truncated
        shape = np.array([NH_max + 1, NV_max + 1])
        data_type = np.dtype('<f4')
        file_name = 'Spec_HV_KE.out'
        spec_HV = self._read_data(file_name, data_type, shape, it)
        kinetic_energy = np.sum(spec_HV[0, :])
        file_name = 'Spec_HV_PE.out'
        spec_HV = self._read_data(file_name, data_type, shape, it)
        potential_energy = np.sum(spec_HV[0, :])

        return kinetic_energy, potential_energy

    def horizontam_mean_energy_series(self, it_start, it_end):
        NH_max = int(np.sqrt(self.r) * self.n1_truncated + 0.5)
        NV_max = self.n3_truncated
        shape = np.array([NH_max + 1, NV_max + 1])
        data_type = np.dtype('<f4')
        KE_file_name = 'Spec_HV_KE.out'
        PE_file_name = 'Spec_HV_PE.out'
        KE_series = np.zeros(it_end - it_start)
        PE_series = np.zeros(it_end - it_start)
        for it in range(it_start, it_end):
            spec_HV_KE = self._read_data(KE_file_name, data_type, shape, it)
            KE_series[it - it_start] = np.sum(spec_HV_KE[0, :])
            spec_HV_PE = self._read_data(PE_file_name, data_type, shape, it)
            PE_series[it - it_start] = np.sum(spec_HV_PE[0, :])
            if it % 200 == 0:
                print(it, '/', it_end,
                      KE_series[it - it_start], np.sum(spec_HV_KE),
                      PE_series[it - it_start], np.sum(spec_HV_PE))

        return KE_series, PE_series

    def get_spectrum_HV(self, it, variable_name):
        file_name = 'Spec_HV_' + variable_name + '.out'
        NH_max = int(np.sqrt(self.r) * self.n1_truncated + 0.5)
        NV_max = self.n3_truncated
        shape = np.array([NH_max + 1, NV_max + 1])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape, it)

        delta_1 = 2 * np.pi / self.l1
        delta_3 = 2 * np.pi / self.l3

        K_axis = np.arange(0, delta_1 * (NH_max + 1), delta_1)
        M_axis = np.arange(0, delta_3 * (NV_max + 1), delta_3)
        return data, K_axis, M_axis

    def get_spectrum_1d(self, it, variable_name):
        file_name = 'Spec_HV_' + variable_name + '.out'
        NH_max = int(np.sqrt(self.r) * self.n1_truncated + 0.5)
        NV_max = self.n3_truncated
        shape = np.array([NH_max + 1, NV_max + 1])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape, it)

        delta_1 = 2 * np.pi / self.l1
        delta_3 = 2 * np.pi / self.l3

        delta = delta_3
        N_max = int(round(
            np.sqrt((NH_max * delta_1)**2 + (NV_max * delta_3)**2) / delta))
        data_out = np.zeros(N_max+1)

        K_axis = np.arange(0, delta_1 * (NH_max + 1), delta_1)
        M_axis = np.arange(0, delta_3 * (NV_max + 1), delta_3)
        M, K = np.meshgrid(M_axis, K_axis)
        index = np.round(np.sqrt(K**2 + M**2) / delta)
        for j in range(N_max+1):
            data_out[j] = np.sum(data[index == j])
        wavenumber = np.arange(0, delta * (N_max + 1), delta)
        data_out = data_out / delta
        return data_out, wavenumber

    def get_spectrum_V(self, it, variable_name):
        file_name = 'Spec_HV_' + variable_name + '.out'
        NH_max = int(np.sqrt(self.r) * self.n1_truncated + 0.5)
        NV_max = self.n3_truncated
        shape = np.array([NH_max + 1, NV_max + 1])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape, it)

        delta_3 = 2 * np.pi / self.l3
        M_axis = np.arange(0, delta_3 * (NV_max + 1), delta_3)
        data_out = np.sum(data, axis=0)
        print('Z spectrum sum:', np.sum(data_out), data_out[0])
        data_out = data_out / delta_3
        return data_out, M_axis

    def get_spectrum_H(self, it, variable_name):
        file_name = 'Spec_HV_' + variable_name + '.out'
        NH_max = int(np.sqrt(self.r) * self.n1_truncated + 0.5)
        NV_max = self.n3_truncated
        shape = np.array([NH_max + 1, NV_max + 1])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape, it)

        delta_1 = 2 * np.pi / self.l1
        K_axis = np.arange(0, delta_1 * (NH_max + 1), delta_1)
        data_out = np.sum(data, axis=1)
        data_out = data_out / delta_1
        return data_out, K_axis

    def get_spectrum_H2(self, it, variable_name):
        file_name = 'Spec_H2_' + variable_name + '.out'
        NH_max = int(np.sqrt(self.r) * self.n1_truncated + 0.5)
        NH_min = int(self.n1_truncated / np.sqrt(self.r) + 0.5)
        shape = np.array([NH_max + 1, NH_min*2 + 1])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape, it)
        delta = 2 * np.pi / self.l1
        K1_axis = np.arange(- delta * NH_min, delta * (NH_min + 1),
                            delta)
        K2_axis = np.arange(0, delta * (NH_max + 1), delta)
        return data, K1_axis, K2_axis

    def get_spectrum_X(self, it, variable_name):
        file_name = 'Spec_H2_' + variable_name + '.out'
        NH_max = int(np.sqrt(self.r) * self.n1_truncated + 0.5)
        NH_min = int(self.n1_truncated / np.sqrt(self.r) + 0.5)
        shape = np.array([NH_max + 1, NH_min*2 + 1])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape, it)
        data[1:, :] = data[1:, :] * 2
        delta = 2 * np.pi / self.l1
        data_out = np.zeros(NH_min+1)
        data_out = np.sum(data[:, NH_min:], axis=0)
        data_out[1:] = data_out[1:] + np.sum(data[:, NH_min-1::-1], axis=0)
        K1_axis = np.arange(0, delta * (NH_min + 1), delta)
        print('X spectrum sum:', np.sum(data_out), data_out[0])
        data_out = data_out / delta
        return data_out, K1_axis

    def get_spectrum_Y(self, it, variable_name):
        file_name = 'Spec_H2_' + variable_name + '.out'
        NH_max = int(np.sqrt(self.r) * self.n1_truncated + 0.5)
        NH_min = int(self.n1_truncated / np.sqrt(self.r) + 0.5)
        shape = np.array([NH_max + 1, NH_min*2 + 1])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape, it)
        data[1:, :] = data[1:, :] * 2
        delta = 2 * np.pi / self.l1
        K2_axis = np.arange(0, delta * (NH_max + 1), delta)
        data_out = np.sum(data, axis=1)
        print('Y spectrum sum:', np.sum(data_out), data_out[0])
        data_out = data_out / delta
        return data_out, K2_axis

    def get_XYZ_G_axis_mpi(self, it, rank, thin):
        LX = self.l1
        LY = self.l2
        LZ = self.l3

        n2R_local = self.n2 // self.division[1]

        X_lin = np.linspace(-LX/2, LX/2, self.n1+1,
                            dtype='float64', endpoint=True)
        Y_lin = np.linspace(-LY/2, LY/2, self.n2+1,
                            dtype='float64', endpoint=True)
        Z_lin = np.linspace(-LZ/2, LZ/2, self.n3+1,
                            dtype='float64', endpoint=True)
        Y_lin_local = Y_lin[rank*n2R_local:(rank+1)*n2R_local+1]
        X, Y, Z = np.meshgrid(
            X_lin[::thin], Y_lin_local[::thin], Z_lin[:],
            indexing='ij')

        it_state = it * self.output_state_interval
        A_inv = np.linalg.inv(self.A[it_state, :, :])
        XG = A_inv[0, 0] * X + A_inv[0, 1] * Y + A_inv[0, 2] * Z
        YG = A_inv[1, 0] * X + A_inv[1, 1] * Y + A_inv[1, 2] * Z
        ZG = A_inv[2, 0] * X + A_inv[2, 1] * Y + A_inv[2, 2] * Z

        XG = np.reshape(XG, np.shape(XG), order='F')
        YG = np.reshape(YG, np.shape(YG), order='F')
        ZG = np.reshape(ZG, np.shape(ZG), order='F')
        return XG, YG, ZG

    def read_Richardson(self):
        shape = np.array([self.N_Ri, self.number_of_record])
        data_type = np.dtype('<f4')
        file_name = 'Richardson_distribution.out'
        data = self._read_data(file_name, data_type, shape)
        return data

    def read_Richardson_2(self, it):
        shape = np.array([self.N_Ri, it])
        data_type = np.dtype('<f4')
        file_name = 'Richardson_distribution.out'
        data = self._read_data(file_name, data_type, shape)
        return data

    def set_Richardson_axis(self):
        Del_Ri = (self.Ri_max - self.Ri_min) / \
            (self.N_Ri - 2)
        Ri_axis = np.linspace(
            self.Ri_min - Del_Ri / 2, self.Ri_max + Del_Ri / 2,
            self.N_Ri, endpoint=True)
        return Ri_axis

    def read_Froude_spectrum(self, it):
        shape = np.array([self.n3_truncated+1])
        data_type = np.dtype('<f4')
        file_name = 'Froude_spectrum.out'
        data = self._read_data(file_name, data_type, shape, it)
        delta_3 = 2 * np.pi / self.l3
        M_axis = np.arange(0, delta_3 * (self.n3_truncated + 1), delta_3)
        return data, M_axis

    def read_Thorpe_displacement(self):
        shape = np.array([self.N_Thorpe, self.number_of_record])
        data_type = np.dtype('<f4')
        file_name = 'Thorpe_displacement.out'
        data = self._read_data(file_name, data_type, shape)
        delta = self.l3 / self.n3
        Thorpe_axis = np.arange(- self.l3 + delta,
                                self.l3, delta)
        return data, Thorpe_axis

    def get_snapshot(self, variable_name, i1, it, surface='XZ'):
        file_name = ('/snapshot/' + surface + '_' + variable_name
                     + '{0:04d}'.format(i1) + '.out')
        shape = np.array([])
        if surface == 'XY':
            shape = np.array([self.n1, self.n2])
        elif surface == 'XZ':
            shape = np.array([self.n1, self.n3])
        data_type = np.dtype('<f4')
        data = self._read_data(file_name, data_type, shape, it)
        return data[:, :]
