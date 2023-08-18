import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

filenames = ['surfaceReactions_O_on_W.nc', 'surfaceReactions_D_on_W.nc', 'surfaceReactions_W_on_W.nc']

models  = ['O on W', 'D on W', 'W on W']
fig, ax = plt.subplots(3, 2, figsize=(10, 10))

for it, filename in enumerate(filenames):
    print(filename)
    data = nc.Dataset("BCA/" + filename, "r")

    rfyld = data['rfyld'][:]  # shape (nA, nE)
    spyld = data['spyld'][:] # (nA, nE)

    E = data['nE'][:]
    A = data['nA'][:]

    print(E.shape, A.shape, rfyld.shape, spyld.shape)
    print(spyld)
    subsampling_factor = 2  # Adjust the subsampling factor as needed

    for it_angle in range(0, len(A), int(len(A) / 4)):
        ax[it][0].semilogx(E[::subsampling_factor], spyld[it_angle, ::subsampling_factor], label=str(round(A[it_angle], 1))+' deg')
        ax[it][1].semilogx(E[::subsampling_factor], rfyld[it_angle, ::subsampling_factor], label=str(round(A[it_angle], 1))+' deg')
        # # ax[it][2].semilogx(E[::subsampling_factor], spyld[it_angle, ::subsampling_factor], label=str(A[it_angle])+' deg')

        ax[it][0].set_xlabel(r'$E$ [eV]')
        ax[it][1].set_xlabel(r'$E$ [eV]')

        ax[0][0].set_ylabel(r'R')
        ax[1][0].set_ylabel(r'R')
        ax[2][0].set_ylabel(r'R')
        ax[0][1].set_ylabel(r'Y')
        ax[1][1].set_ylabel(r'Y')
        ax[2][1].set_ylabel(r'Y')


        # # annotating the plots 
        ax[it][0].set_title(str(models[it])) 
        ax[it][1].set_title(str(models[it]))
        ax[it][0].grid(True, which="both", ls="-")
        ax[it][1].grid(True, which="both", ls="-")

        ax[0][0].legend(loc='best')
        ax[1][0].legend(loc='best')
        ax[2][0].legend(loc='best')

        plt.tight_layout()

plt.savefig('SurfaceReactionsModels.png',dpi=300, bbox_inches = "tight",facecolor='w')

plt.show()
