
# SPM import
import nipype.interfaces.spm as spm

# Other import
import os.path
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import math
import scipy.io as spio # open .mat file
import os
from scipy import signal
import pandas as pd
import scipy.stats as ss

class ASL_processes():

    def __init__(self):
        super(ASL_processes,self).__init__()

    def Rescaling_process(self, params):

        NSlices = params['image'].shape[2]
        Rec_matrix = np.array([params['image'].shape[0],params['image'].shape[1]])
        VoxSasl = np.product(params['image'].header.get_zooms()[0:3])*1e-9
        VoxVasl = VoxSasl
        VoxVm0 =  VoxSasl
        PASLModes = ['FAIR','TILT','STAR','Q2TIPS']
        CASLModes = ['CASL','pCASL','PSEUDO']
        deltaTI = []
        SFasl = None
        SFm0_1 = None
        SFm0_2end = None

        # TR_slice calculation
        if params['scan_mode'] == 3:
            TR_slice = 0
        else:
            TR_slice = 6.4*1e-3 + 2*params['TE']

        # label efficienty changes
        if not params['label_eff']:
            if not params['ASL_mode']:
                params['ASL_mode'] = 'pCASL'
            if params['ASL_mode'] == 'pCASL':
                params['label_eff'] = 0.85
            else:
                params['label_eff'] = 0.98
        if not params['back_supp']:
            if not params['M0_image_file']:
                params['back_supp'] = 0
            else:
                params['back_supp'] = 4
        if params['back_supp'] != 0:
            params['label_eff'] = params['label_eff']*(0.95**params['back_supp'])

        # calculation of scale factors
        if params['method'] == 'whitepaper':
            deltaTI = np.reshape(np.array(range(NSlices)) * TR_slice,[1,1,NSlices])
            SFasl = (VoxVasl/params['lambda_'])*2*params['label_eff']*np.exp(-(params['PLD'] + deltaTI)/params['T1_blood'])
            if params['ASL_mode'] in PASLModes:
                SFasl = SFasl * params['TI1'] * np.exp(-params['TI1'] / params['T1_blood'])
            elif params['ASL_mode'] in CASLModes:
                SFasl = SFasl * params['T1_blood'] * (1 - np.exp(-params['TI1'] / params['T1_blood']))
            SFasl = np.float32(SFasl)/6  # from kg*s to hg*min

            # set scaling parameters for M0 data
            SFm0_1 = VoxVm0
            SFm0_2end = SFm0_1 * (1-np.exp(-params['TRm0']/params['T1t']))
            SFm0_1 = np.float32(SFm0_1)*1e6  # from m^3 to mL
            SFm0_2end = np.float32(SFm0_2end)*1e6

            # TODO noeth

            return SFasl, SFm0_1, SFm0_2end

    def Realigning_process(self, params):
        """
        Realigns ASL images
        :param variable_name: params
        :param value: dictionary of input parameters
        :return: ASLOutFnames, graphicsFname and excludeFname
        """

        filename, file_extension = os.path.splitext(params['ASL_file'])
        path, file = os.path.split(filename)
        R_orig = params['image'].affine #I am not sure if this is getting the data needed
        realign = spm.Realign()
        realign.inputs.jobtype = 'estimate'
        realign.inputs.in_files = params['ASL_file']
        realign.inputs.quality = params['quality']
        realign.inputs.fwhm = params['FWHM']
        realign.inputs.separation = params['separation']
        realign.inputs.register_to_mean = params['register_to_mean']
        realign.inputs.wrap = params['wrapping']
        if params['weighting']:
            realign.inputs.weight_img= params['weighting']
        results = realign.run()

        # get rp data from .txt
        rp = np.loadtxt(results.outputs.realignment_parameters, dtype = float)
        rp = rp.transpose()

        # get ctl and tag data
        ctl_data = rp[:,::2]
        tag_data = rp[:,1::2]

        # estimate average difference between control and tag images
        rp_corr = np.mean(ctl_data, axis = 1) - np.mean(tag_data, axis = 1)
        rp_corr = np.tile(rp_corr, (rp.shape[1], 1)).transpose()/2  # rp_corr = (np.array([rp_corr,]*rp.shape[1]).transpose())/2
        rp_corr[:, ::2] = -rp_corr[:, ::2]
        rp2 = rp + rp_corr

        # first, detect local outliers in motion parameters
        ddrp = np.diff(rp2, 2, axis = 1)/2
        ddrp[3:6,:] = (ddrp[3:6,:]*180)/np.pi
        ddrp = np.pad(np.sqrt(np.sum(np.power(ddrp, 2), axis = 0)), (1,1), 'constant', constant_values=(0, 0)) #HERE THE VALUES BEGIN TO GROWWW
        outliers = ddrp > params['motion_threshold']

        # second, detect large absolute motion wrt the average motion parameters
        mrp2 = np.mean(rp2[:, ~outliers], axis = 1)
        mrp2 = np.expand_dims(mrp2, axis = 1)
        mrp2 = np.repeat(mrp2, rp.shape[1], axis = 1)

        drp = rp2 - mrp2
        drp[3:6,:] = (drp[3:6,:]*180) / np.pi
        drp = np.sqrt(np.sum(np.power(drp, 2), axis = 0))/2
        exclude = np.sqrt(np.power(ddrp, 2) + np.power(drp, 2)) > params['motion_threshold']

        # now, repeat the estimation of difference between control and tag,
        # based on the non-excluded frames
        ctl_ex = exclude[::2]
        tag_ex = exclude[1::2]
        rp_corr = np.mean(ctl_data.transpose()[~ctl_ex, :], axis = 0) - np.mean(tag_data.transpose()[~tag_ex, :], axis = 0)

        # subtract Pc/2 from all ctl frames and add to tag frames
        rp_corr = np.tile(rp_corr, (rp.shape[1], 1)).transpose()/2
        rp_corr[:,::2] = -rp_corr[:,::2]
        rp2 = rp + rp_corr

        # overwrite params
        i = 0
        with open(results.outputs.realignment_parameters,'r+') as f:
            for line in f:
                line = rp2[:,i]
                i += 1

        # calculate excluded frames
        exclude = np.argwhere(exclude)

        # plot
        if params['save_graphics']:
            fig = plt.figure()
            plt.subplot(2,1,1)
            linestyle = ['r', 'g', 'b']
            linestyle2 = ['r:', 'g:', 'b:']
            for i in range(3):
                plt.plot(rp[i,:], linestyle[i], markersize = 1, linewidth = 0.7)
                plt.plot(rp2[i,:], linestyle2[i], markersize = 1, linewidth = 0.7)

            plt.legend(('x-Rotation', 'x-Rotation', 'y-Rotation','y-Rotation', 'z-Rotation', 'z-Rotation'))
            plt.ylabel = 'Translation [mm]'

            bottom, top = plt.ylim()
            i = 0
            for fno in exclude:
                plt.bar(exclude[i][0], top - bottom, 0.8, bottom = bottom, color = 'grey')
                i += 1

            plt.subplot(2,1,2)
            for i in range(3):
                plt.plot(rp[i+3,:]*180/np.pi, linestyle[i], markersize = 1, linewidth = 0.7)
                plt.plot(rp2[i+3,:]*180/np.pi, linestyle2[i], markersize = 1, linewidth = 0.7)

            bottom, top = plt.ylim()
            i = 0
            for fno in exclude:
                plt.bar(exclude[i][0], top - bottom, 0.8, bottom = bottom, color = 'grey')
                i += 1
            plt.ylabel = 'Rotation [deg]'
            plt.savefig(params['graphics_fname'])

        # if params['save_params']:
        #     print('Writing params to')
        #     print(params['params_fname'])
        #
        # else:
        #     # remove or not input??
        #     print('No')

        # get the affines matrices from .mat
        mat_fname = os.path.join(path, file) + '.mat'
        mat = spio.loadmat(mat_fname, squeeze_me=True)
        aff = mat['mat']

        # apply new realignment
        for i in range(rp2.shape[1]):
            Mrea2 = self.spm_matrix(rp2[:,i])
            M2 = np.dot(Mrea2,R_orig)
            aff[:,:,i] = M2
        mat['mat'] = aff
        spio.savemat(mat_fname, mat)

        # Now that we have the correct realignment parameters, write the requested images
        realign = spm.Realign()
        realign.inputs.jobtype = 'write'
        realign.inputs.in_files = params['ASL_file']
        realign.inputs.interp = params['interpolation']
        realign.inputs.wrap = params['wrapping']
        realign.inputs.write_which = params['write_which']
        realign.inputs.out_prefix = params['out_prefix']
        realign.inputs.write_mask = params['write_mask']
        realign.run()

        return exclude

    def Outlier_detection_process(self, params):
        """
        Add a new tag with the outliers, volumes with different mean or STD
        :param variable_name: params
        :param value: dictionary of input parameters
        :return: ASL and M0 images with added tag
        """

        filename, file_extension = os.path.splitext(params['im_file'])
        path, filename = os.path.split(filename)
        im = nib.load(params['im_file'])
        img_data = im.get_fdata()

        NSlices = img_data.shape[3]
        NVoxel = np.int(np.prod(img_data.shape)/NSlices)
        img_data = np.reshape(img_data,(NVoxel, NSlices), order="F")

        if params['mask_file']: # TODO: FIND DATA TO TEST IT
            mask = nib.load(params['mask_file'])
            mask_data = mask.get_fdata()
            mask_data = np.nonzero(mask_data)
            valid_vox = ~np.any(np.isnan(img_data[mask_data,:]), axis=1)
            img_data = img_data[valid_vox,:]  # valid 4D data

            # valid_vox = ~np.any(np.isnan(img_data), axis = 1)
        else:
            valid_vox = ~np.any(np.isnan(img_data), axis = 1)
            img_data = img_data[valid_vox,:]  # valid 4D data

        # get mean and std
        img_mean = np.mean(img_data, axis = 0)
        img_std = np.std(img_data, axis = 0)

        # remove ctl/tag difference
        img_mean[::2] =  img_mean[::2] + np.mean(img_mean[1::2] - img_mean[::2], axis = 0)
        img_std[::2] =  img_std[::2] + np.mean(img_std[1::2] - img_std[::2], axis = 0)

        # remove low-frequency drift that will be ignored in the SPM analysis anyway
        RT = 4
        row = np.arange(1, len(img_mean) + 1, 1)
        cutoff = 48
        img_M_hpf = self.spm_filter(RT, row, cutoff, img_mean)
        img_S_hpf = self.spm_filter(RT, row, cutoff, img_std)

        # detect global outliers with Thompson tau method
        idx_Mgs, lim1 = self.thompson_tau(img_M_hpf, 0.0005, len(img_M_hpf)) # index of mean outliers
        idx_Sgs, lim2 = self.thompson_tau(img_S_hpf, 0.0005, len(img_S_hpf)) # index of std outliers
        idx_global_outlier = np.union1d(idx_Mgs, idx_Sgs)

        # detect global outliers at a more moderate level, to intersect with local outlier detection
        idx_Mg = self.thompson_tau(img_M_hpf, 0.005, len(img_M_hpf)) # index of mean outliers
        idx_Sg = self.thompson_tau(img_S_hpf, 0.005, len(img_S_hpf)) # index of std outliers

        _,idx_Ml = self.hampel(img_mean, 8, 4)
        _,idx_Sl = self.hampel(img_std, 8, 4)

        idx_Ml = list(set(np.nonzero(idx_Ml)[0]) & set(idx_Mg[0]))
        idx_Sl = list(set(np.nonzero(idx_Sl)[0]) & set(idx_Sg[0]))
        idx_local_outlier = np.union1d(idx_Ml,idx_Sl)

        idx_outlier = np.union1d(idx_global_outlier,idx_local_outlier)

        # plot results
        if params['plot_diagnostics']:
            fig = plt.figure()
            plt.subplot(2,1,1)
            plt.title('Outliers')
            plt.plot(img_M_hpf, 'bo', markersize = 1.5)
            if idx_Mgs != []:
                plt.plot(idx_Mgs, img_M_hpf[idx_Mgs], 'rx') # plot global outliers
                plt.plot(idx_Ml, img_M_hpf[idx_Ml], 'rx') # plot local outliers
            plt.hlines(lim1, 0, len(img_mean), label = 'Threshold')

            plt.subplot(2, 1, 2)
            plt.plot(img_S_hpf, 'bo', markersize = 1.5)
            if idx_Sgs != []:
                plt.plot(idx_Sgs, img_S_hpf[idx_Sgs], 'rx') # plot global outliers
                plt.plot(idx_Sl, img_S_hpf[idx_Sl], 'rx') # plot local outliers
            plt.hlines(lim2, 0, len(img_std), label = 'Threshold')

        plt.savefig(os.path.join(path, 'o_plot_' + filename) + '.png')
        return idx_outlier

    def Generate_design_matrix_process(self, params):
        """
        Generates the design matrix including ASL regressor and excluded frames.
        This writes a .txt file containing the matrix 'R', which is then passed to the
        fMRI model specification brick from SPM
        :param variable_name: params
        :param value: dictionary of input parameters
        :return: file name of the generated .txt file
        """
        params['phys_reg_file'] = 'toto'
        if params['phys_reg_file']:
            print('The program is not getting there, so as I dont have this file, i cannot code here')

        # load excluded frames and add a regressor per frame


        return None

    def Generate_TMap_process(self, params):  # possibly handle as a call to Results report ?
        """
        Creates T maps
        :param params: params
        :param value: dictionary of input parameters
        :return: result of the call to spm_run_con
        """
        tmp = spio.loadmat(params['SPM_file'], squeeze_me=True)
        tmp = mat['mat']
        return None

    def Calculate_CBFmap_process(self, params):
        """
        Apply the final scaling to SPM beta maps to obtain quantitative CBF(-change) maps.
        :param params: params
        :param value: dictionary of input parameters
        :return: names of files containing generated CBF(-change) maps   #Files, instead of filenames???
        """
        return None

    #Utility functions

    def fiveD_to_fourD(self,img_data):
        Temp = []
        if len(img_data.shape) == 5:  # there is a 5th dimension
            Temp = np.empty((img_data.shape[0], img_data.shape[1],
                             img_data.shape[2], img_data.shape[3] * img_data.shape[4]))
            for i in range(img_data.shape[4]):
                for j in range(img_data.shape[3]):
                    Temp[:, :, :, i + j * 2] = img_data[:, :, :, j, i]
            img_data = Temp

        return np.array(img_data)

    def Apply_SF(self, img_data, SF):
        """
        Apply the scale factor to the image data. If the scale factor matrix does
        not match the image data matrix, reshapes the scale factor matrix.
        :param im: Image
        :param value: Object
        :param SF: Scale factor matrix
        :param value: Matrix
        :return: Scaled ASL Image
        """
        n_image_dims = len(img_data.shape)
        for dim in range(len(SF.shape), n_image_dims):
            SF = np.expand_dims(SF, axis = dim)  # as it is defined tile, we need to create first a new dimension

        if np.any(np.logical_and(SF.shape != np.array(img_data.shape), SF.shape != np.ones(n_image_dims))):
            # error !!
            print('Illegal shape for scale factor parameter.')
            print('Shape of image data: ', img_data.shape)
            print('Shape of scale factor: ', SF.shape)
            return None

        repetition_vector = np.int16(np.array(img_data.shape) / np.array(SF.shape))
        SF = np.tile(SF, repetition_vector)

        return img_data / SF

    def rang(self, x):
        return np.minimum(np.maximum(x, -1), 1)

    def spm_matrix(self, P):
        '''
        :param P: array with the parameters:
            P[0]  - x translation
            P[1]  - y translation
            P[2]  - z translation
            P[3]  - x rotation about - {pitch} (radians)
            P[4]  - y rotation about - {roll}  (radians)
            P[5]  - z rotation about - {yaw}   (radians)
            P[6]  - x scaling
            P[7]  - y scaling
            P[8]  - z scaling
            P[9]  - x affine
            P[10] - y affine
            P[11] - z affine
        :return: affine transformation matrix
        '''
        #TO TEST IT
        # P = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0])
        # aff_matrix = np.array([[1, 0, 0, 0],
        #                       [0, 1, 0, 0],
        #                       [0, 0, 1, 0],
        #                       [0, 0, 0, 1]])
        q  = [0,0,0,0,0,0,1,1,1,0,0,0]
        if len(P)< 12:
            q[0:len(P)]=P
            P = np.array(q)
        if np.prod(P.shape) == 3:
            A = np.eye(4)
            A[0:3, 3] = P[:]

        T = np.array([[1, 0, 0, P[0]],
                      [0, 1, 0, P[1]],
                      [0, 0, 1, P[2]],
                      [0, 0, 0, 1]])

        R1 = np.array([[1, 0, 0, 0],
                       [0, math.cos(P[3]), math.sin(P[3]), 0],
                       [0, -math.sin(P[3]), math.cos(P[3]), 0],
                       [0, 0, 0, 1]])

        R2 = np.array([[math.cos(P[4]), 0, math.sin(P[4]), 0],
                       [0, 1, 0, 0],
                       [-math.sin(P[4]), 0, math.cos(P[4]), 0],
                       [0, 0, 0, 1]])

        R3 = np.array([[math.cos(P[5]), math.sin(P[5]), 0, 0],
                       [-math.sin(P[5]), math.cos(P[5]), 0, 0],
                       [0, 0, 1, 0],
                       [0, 0, 0, 1]])

        R = np.linalg.multi_dot([R1, R2, R3])

        Z = np.array([[P[6], 0, 0, 0],
                      [0, P[7], 0, 0],
                      [0, 0, P[8], 0],
                      [0, 0, 0, 1]])

        S = np.array([[1, P[9], P[10], 0],
                      [0, 1, P[11], 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])

        return np.linalg.multi_dot([T, R, Z, S])

    def spm_imatrix(self, aff_matrix):
        """
        Return the parameters for creating an affine transformation matrix
        :param aff_matrix: Affine transformation matrix
        :return: An array of parameters:
            P[0]  - x translation
            P[1]  - y translation
            P[2]  - z translation
            P[3]  - x rotation about - {pitch} (radians)
            P[4]  - y rotation about - {roll}  (radians)
            P[5]  - z rotation about - {yaw}   (radians)
            P[6]  - x scaling
            P[7]  - y scaling
            P[8]  - z scaling
            P[9]  - x affine
            P[10] - y affine
            P[11] - z affine
        """
        # TO TEST IT
        # aff_matrix = np.array([[1.0, 0.0002, -0.0003, 0.0184],
        #                       [-0.0002, 1.0, -0.0002, 0.016],
        #                       [0.0003, 0.0002, 1.0, 0.0019],
        #                       [0, 0, 0, 1.0]])
        # ans = np.array([0.0184, 0.016, 0.0019, -0.0002, -0.0003, 0.0002, 1.0, 1.0, 1.0, -0, 0, 0])

        P = np.zeros(12)
        R = aff_matrix[0:3,0:3]
        C = np.linalg.cholesky(np.transpose(R)*R)
        P[0:3] = aff_matrix[0:3,3]
        P[6:9] = np.diag(C)

        if np.linalg.det(R) < 0:
            P[6] = -P[6]

        C = np.linalg.solve(np.diag(np.diag(C)), C)
        P[9] = C[0,1]
        P[10] = C[0,2]
        P[11] = C[1,2]
        R0 = self.spm_matrix(np.array([0, 0, 0, 0, 0, 0, P[6], P[7], P[8], P[9], P[10], P[11]]))
        R0 = R0[0:3,0:3]
        R1 = np.dot(R,np.linalg.inv(R0))
        P[4] = math.asin(rang(R[0,2]))
        if np.power(abs(P[4]) - np.pi/2, 2) < 1e-9:
            P[3] = 0
            P[6] = math.atan2(-rang(R[1,0]), rang(-R1[2,0]/R1[2,0]))
        else:
            c = math.cos(P[4])
            P[3] = math.atan2(rang(R1[1,2]/c), rang(R1[2,2]/c))
            P[5] = math.atan2(rang(R1[0,1]/c), rang(R1[0,0]/c))

        return P

    def maxim(self, data):
        '''
        Calculates de maximum value in a list with nan values
        Inputs:
            - data: array-like.
        Outputs:
            - m: int, the maximun value in the array.
            - ind: the index of the maximum.
        '''
        m = 0
        ind = 0
        for i in range(len(data)):
            if not (np.isnan(data[i])) and data[i] > m:
                m = data[i]
                ind = i
        return m, ind

    def thompson_tau(self, data, alpha, threshold):
        # Information obtained from: https://www.statisticshowto.datasciencecentral.com/modified-thompson-tau-test/
        '''
        Implements the Thompson Tau method and returns a list with the outliers index.
        Inputs:
            - data: an array.
            - alpha: the significance level.
            - Threshold: the number of points tested.
        Outputs:
             - outliers: a list with the indices of the outliers.
        '''
        outliers = []
        n = len(data)
        mean = np.mean(data)
        delta = abs(data - mean)
        std = np.std(data)
        for i in range(threshold):
            d, ind = self.maxim(delta)
            reject = ss.t.ppf(alpha / 2, n - 2)
            tau = (reject * (n - 1)) / (np.sqrt(n) * np.sqrt(n - 2 + np.power(reject, 2)))
            if d > -tau * std or d < tau * std: # MIRAR BIEN SI ESTO ES CORRECTO
                outliers += [ind]
            delta[ind] = None
        thr = np.array([mean - tau * std, mean + tau * std])
        return np.array(outliers), thr

    def hampel(self, data, k=7, t0=3): # https://stackoverflow.com/questions/46819260/filtering-outliers-how-to-make-median-based-hampel-function-faster
        # maybe these are better: https://ocefpaf.github.io/python4oceanographers/blog/2015/03/16/outlier_detection/
        '''
        Inputs:
            vals: pandas series of values from which to remove outliers
            k: size of window (including the sample; 7 is equal to 3 on either side of value)
        Outputs:
            vals: array, the data filtered.
            outlier_idx_list: logical index of the replaced values.
        '''
        outlier_idx_list = []
        # convert to Series to make it more efficient
        vals_orig = pd.Series(data)

        #Make copy so original not edited
        vals=vals_orig.copy()

        #Hampel Filter
        L = 1.4826
        rolling_median = vals.rolling(k).median()
        difference = np.abs(rolling_median-vals)
        median_abs_deviation = difference.rolling(k).median()
        threshold = t0 * L * median_abs_deviation
        outlier_idx = difference > threshold
        outlier_idx_list += [outlier_idx]
        vals[outlier_idx] = np.nan
        vals = vals.values
        outlier_idx_list = np.array(outlier_idx_list)
        return vals, outlier_idx_list

    def spm_dctmtx(self, N, K):
        # N = k K = n
        t = np.arange(0, N, 1)
        C = np.zeros((len(t), K))
        C[:,0] = np.ones((1, len(t)))/np.sqrt(N)
        for i in range(1, K, 1):
            C[:,i] = np.sqrt(2/N)*np.cos(np.pi*(2*t+1)*i/(2*N))
        return C

    def spm_filter(self, RT, row, cutoff, Y): #TODO: make this useful not just for this case
        k = len(row)
        n = int(2*(k*RT)/cutoff + 1)
        X0 = self.spm_dctmtx(k,n)
        K_X0 = X0[:,1:]

        # apply HPF
        Y = Y - np.linalg.multi_dot([K_X0,K_X0.transpose(),Y])
        return Y

    def Resample_image(self, im): # Check if this exists in nibabel! --> resample_from_to,resample_to_output???
        return None


class Anat:
    def __init__(self):  # get the params, the Anat included
        return None

    def Generate_mask_process(self, params):
        """
        GFB_create_masks: Create vascular and anatomic masks using the intersection between
        vascular and anatomical masks with tissue probability maps. The binary masks are
        optionally resampled to the geometry of a reference volume
        :param variable_name: params
        :param value: dictionary of input parameters
        :return: file of the generated mask image (TODO: check this)
        """

        if params["normalised_data"] == False:
            if params["normalisation_parameters"][-4:] == '.mat':
                #TODO: use spm_imatrix
                print('It is not normalised')

        if params["reference_volume"] != '':
            params["tissue_map_data"] = params["tissue_map_data"] > params["tissue_map_threshold"]

        mask_vol = nib.load(params["masks"])

        # get mask_vol for non-normalised data
        if params["normalised_data"] == False: #TODO: Finish
            if params["normalisation_parameters"][-4:] == '.mat':
                print("Finish, I dont know the normParams")
            else:
                print("Finish, I dont know the normParams")

        mask_vol_data = mask_vol.get_fdata()
        mask_vol_data = mask_vol_data > 0.5
        print(mask_vol_data.shape,params["tissue_map_data"].shape)
        ROI = mask_vol_data*params["tissue_map_data"]

    # Utility functions
    def Get_Def(self, im):
        """
        Get the deformation field of the image.
        :param im: Image
        :param value: Object
        :return: Deformation field matrix
        """
        return None
    def Apply_Def(self, def_, im):
        """
        Apply deformation field to the image data.
        :param im: Image
        :param value: Object
        :return: Deformation field matrix applied
        """
        return None



# import pickle
# f1 = open('/home/nietob/Documents/matrix1.txt','wb')
# f2 = open('/home/nietob/Documents/matrix2.txt', 'wb')
# pickle.dump(rp2, f1)
# pickle.dump(R_orig, f2)
# f1.close()
# f2.close()

