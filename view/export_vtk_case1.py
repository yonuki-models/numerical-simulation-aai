# import sys
import os
import numpy as np
from mpi4py import MPI
from mpi4py_fft import PFFT, newDistArray

import vtk

# sys.path.append(os.path.abspath("../analysis/"))
from read_data import ReadClass

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()

setting_name = "exRo10N3e06_LES"
reader = ReadClass(setting_name)
reader.read_geometry(0)

vtk_directory = reader.data_directory + 'vtk/'
if rank == 0 and (not os.path.isdir(vtk_directory)):
    os.makedirs(vtk_directory)

str_rank = '{0:03d}'.format(rank)
thin = 1

n2R_local = reader.n2 // size
n3R_local = reader.n3 // size
global_shape = np.array([reader.n1, reader.n2, reader.n3], dtype=int)
local_shape = np.array([reader.n1, n2R_local, reader.n3], dtype=int)

c2r = PFFT(comm, global_shape, axes=(2, 1, 0), dtype=float, grid=(1, 1, -1))
C_fft = newDistArray(c2r, True)
R_fft = newDistArray(c2r, False)
R_recv = np.zeros([size, reader.n1, n2R_local, n3R_local], order='F')
R = np.zeros(local_shape, order='F')

shape_R_cell = local_shape
shape_R_cell_thinned = shape_R_cell // [thin, thin, 1]
shape_R_point = shape_R_cell + 1
shape_R_point_thinned = shape_R_cell // [thin, thin, 1] + 1

XG = np.zeros(shape_R_point_thinned, order='F')
YG = np.zeros(shape_R_point_thinned, order='F')
ZG = np.zeros(shape_R_point_thinned, order='F')

XA = np.zeros(np.prod(shape_R_point_thinned), order='F')
YA = np.zeros(np.prod(shape_R_point_thinned), order='F')
ZA = np.zeros(np.prod(shape_R_point_thinned), order='F')

D_ghost_send = np.zeros([reader.n1, reader.n3], order='F')
D_ghost_recv = np.zeros([reader.n1, reader.n3], order='F')
D_ext = np.zeros(shape_R_point, order='F')
# D = np.zeros(shape_R_point_thinned, order='F')

UA = np.zeros(np.prod(shape_R_point_thinned), order='F')
VA = np.zeros(np.prod(shape_R_point_thinned), order='F')
WA = np.zeros(np.prod(shape_R_point_thinned), order='F')
TA = np.zeros(np.prod(shape_R_point_thinned), order='F')

# U = np.zeros(shape_R_point_thinned, order='F')
# V = np.zeros(shape_R_point_thinned, order='F')
# W = np.zeros(shape_R_point_thinned, order='F')
# T = np.zeros(shape_R_point_thinned, order='F')


def create_pvts(template, file_name_without_ext):
    nX = shape_R_cell_thinned[0]
    nY = shape_R_cell_thinned[1]
    nZ = shape_R_cell_thinned[2]

    with open(template, 'r') as f:
        text = f.read()

        replace_text = text.replace(
            'WholeExtent=""',
            'WholeExtent="0 ' + str(nX) + ' 0 ' + str(nY*size)
            + ' 0 ' + str(nZ) + '"'
        )
        insert_text = ''
        for j in range(size):
            insert_text += (
                '\n    <Piece Extent="0 ' + str(nX) + ' '
                + str(j*nY) + ' ' + str((j + 1)*nY) + ' 0 ' + str(nZ) + '"'
                + ' Source="' + file_name_without_ext + '_'
                + str(j) + '.vts"/>'
            )
        replace_text = replace_text.replace('TO_BE_REPLACED', insert_text)

    with open(vtk_directory + file_name_without_ext + '.pvts', 'w') as f:
        f.write(replace_text)


def read_scalar_data(variable_name, it):
    C_fft[:, :, :] = reader.get_complex_data_mpi(it, variable_name, rank)
    R_fft[:, :, :] = c2r.backward(C_fft)

    R_send = []
    for j in range(size):
        R_send.append(R_fft[:, j*n2R_local:(j+1)*n2R_local, :])
    R_recv[:, :, :, :] = comm.alltoall(R_send)
    R[:, :, :] = np.reshape(
        R_recv.transpose(1, 2, 3, 0),
        [reader.n1, n2R_local, reader.n3],
        order='F'
    )

    D_ext[:-1, :-1, :-1] = R[:, :, :]
    D_ghost_send[:, :] = R[:, 0, :]
    comm.Sendrecv(
        sendbuf=D_ghost_send, dest=(size+rank-1) % size,
        recvbuf=D_ghost_recv, source=(rank+1) % size)
    D_ext[:-1, -1, :-1] = D_ghost_recv
    D_ext[-1, :, :-1] = D_ext[0, :, :-1]
    D_ext[:-1, :, -1] = D_ext[:-1, :, 0]
    D_ext[-1, :, -1] = D_ext[0, :, 0]

    return D_ext


def output_vtk(it, XG, YG, ZG, thin=thin):
    XG[:, :, :], YG[:, :, :], ZG[:, :, :] = reader.get_XYZ_G_axis_mpi(
        it, rank, thin)
    # ZG[:, :, :] = ZG[:, :, :] * 10
    XA[:] = np.ravel(XG, order='F')
    YA[:] = np.ravel(YG, order='F')
    ZA[:] = np.ravel(ZG, order='F')

    sgrid = vtk.vtkStructuredGrid()
    sgrid.SetDimensions(shape_R_point_thinned)

    points = vtk.vtkPoints()
    points.Allocate(np.product(shape_R_point_thinned))

    for j in range(np.product(shape_R_point_thinned)):
        points.InsertPoint(j, (XA[j], YA[j], ZA[j]))

    UA[:] = np.ravel(read_scalar_data('U', it)[::thin, ::thin, :], order='F')
    VA[:] = np.ravel(read_scalar_data('V', it)[::thin, ::thin, :], order='F')
    WA[:] = np.ravel(read_scalar_data('W', it)[::thin, ::thin, :], order='F')
    vector = vtk.vtkFloatArray()
    vector.SetNumberOfComponents(3)
    vector.SetNumberOfTuples(np.product(shape_R_point_thinned))
    vector.SetName('U')
    for j in range(np.product(shape_R_point_thinned)):
        vector.InsertTuple(j, [UA[j], VA[j], WA[j]])

    TA[:] = np.ravel(read_scalar_data('T', it)[::thin, ::thin, :], order='F')
    scalar = vtk.vtkFloatArray()
    scalar.SetNumberOfComponents(1)
    scalar.SetName('T')
    scalar.SetNumberOfValues(np.product(shape_R_point_thinned))
    for j in range(np.product(shape_R_point_thinned)):
        scalar.InsertValue(j, TA[j])

    sgrid.SetPoints(points)
    sgrid.GetPointData().SetScalars(scalar)
    sgrid.GetPointData().SetVectors(vector)

    extents = [0, shape_R_cell_thinned[0],
               shape_R_cell_thinned[1] * rank,
               shape_R_cell_thinned[1] * (rank + 1),
               0, shape_R_cell_thinned[2]]
    sgrid.SetExtent(extents)

    str_it = '{0:06d}'.format(it)
    file_name = (vtk_directory + 'UVWT' + str_it + '_'
                 + str(rank) + '.vts')

    writer = vtk.vtkXMLDataSetWriter()
    writer.EncodeAppendedDataOff()
    writer.SetFileName(file_name)
    writer.SetInputData(sgrid)
    writer.Write()

    if (rank == 0):
        template = './Template.pvts'
        file_name_without_ext = 'UVWT' + str_it
        create_pvts(template, file_name_without_ext)
        print(it)


for it in range(0, 100, 1):
    output_vtk(it, XG, YG, ZG)

MPI.Finalize()
