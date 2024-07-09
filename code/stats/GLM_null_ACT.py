import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm, pearsonr, zscore
from scipy.spatial.distance import squareform
from statsmodels.stats.multitest import multipletests
import nibabel as nib
import scipy.io
import os
from sklearn import preprocessing
from sklearn import datasets, linear_model
from sklearn.cluster import KMeans
from scipy.interpolate import interp1d
from scipy import stats
import pandas as pd 
from sklearn.metrics import mean_squared_error


file_path = os.path.dirname(os.path.realpath(__file__))
file_path_try = file_path.split('\\')
if len(file_path_try) == 1:
    file_path_try = file_path.split('/')
    root_path = '/'.join(file_path_try[:-2]) + '/'
else:
    root_path = '\\'.join(file_path_try[:-2]) + '\\'

sys.path.append(os.path.join(root_path, 'tools'))
from isc_standalone import circular_timeshift, phase_randomize

parser = argparse.ArgumentParser(
    description="Calculate Pearson Correlation between original time series and behavioral signal")
parser.add_argument("-i", "--input", type=int, nargs=1, help="regional level: 0=node, 1=edge, 2=state", default=0)
parser.add_argument("-s", "--state",  nargs=3, help="state hyperparameters, need when regional level is state; [K, delta, iter_tmp]", default=[4, 0.3, 5])
parser.add_argument("-p", "--iter", type=int, nargs=1, help="iteration to generate null distribution", default=500)
args = parser.parse_args()

def get_statistics(model, train_df=None, target_df=None):
    params = np.append(model.intercept_, model.coef_)
    prediction = model.predict(train_df)
    if len(prediction.shape) == 1:
        prediction = np.expand_dims(prediction, axis=1)
    ones_m=np.ones((train_df.shape[0],1)) #add ones column
    new_trainset=np.hstack([ones_m,train_df]) 
    MSE = mean_squared_error(prediction, target_df)
    variance = MSE * (np.linalg.inv(np.dot(np.transpose(new_trainset), new_trainset)).diagonal())       # MSE = (1, ) & else = (n, ) 가 나와야 함.
    std_error = np.sqrt(variance)
    t_values = params / std_error
    p_values = [2 * (1 - stats.t.cdf(np.abs(i), (new_trainset.shape[0] - new_trainset.shape[1] - 1))) for i in t_values]
    std_error = np.round(std_error, 6)
    t_values = np.round(t_values, 3)
    p_values = np.round(p_values, 4)
    params = np.round(params, 4)
    return params[:], std_error[:], t_values[:], p_values[:]

task_type = 0
task_label = ['whole', 'loss', 'reward']

# load data(X : behavioral data, Y : neural activation or edge connetivity)
neuro_data_whole = scipy.io.loadmat(os.path.join(root_path, "data", "ts_all.mat"))
beh_data_whole = scipy.io.loadmat(os.path.join(root_path, "data", "behav", "signal_data.mat"))

## behavioral data
prec_hrf_data = beh_data_whole['Precision_hrf']
ppe_hrf_data = beh_data_whole['Pos_prediction_error_hrf']
npe_hrf_data = beh_data_whole['Neg_prediction_error_hrf']

as_hrf_data = beh_data_whole['AS_hrf']
cvs_hrf_data = beh_data_whole['CVS_hrf']
fvs_hrf_data = beh_data_whole['FVS_hrf']
fvs_w_hrf_data = beh_data_whole['FVS_w_hrf']
fvs_l_hrf_data = beh_data_whole['FVS_l_hrf']

pw_hrf_data = beh_data_whole['PW_hrf']
nw_hrf_data = beh_data_whole['NW_hrf']
pl_hrf_data = beh_data_whole['PL_hrf']
nl_hrf_data = beh_data_whole['NL_hrf']

PE_certain_hrf_data = beh_data_whole['PE_certain_hrf']
PE_uncertain_hrf_data = beh_data_whole['PE_uncertain_hrf']
FVS_certain_hrf_data = beh_data_whole['FVS_certain_hrf']
FVS_uncertain_hrf_data = beh_data_whole['FVS_uncertain_hrf']
FVS_uncertain1_hrf_data = beh_data_whole['FVS_uncertain1_hrf']
FVS_uncertain2_hrf_data = beh_data_whole['FVS_uncertain2_hrf']

FVS_uncertain2_1_hrf_data = beh_data_whole['FVS_uncertain2_1_hrf']
FVS_uncertain2_2_hrf_data = beh_data_whole['FVS_uncertain2_2_hrf']
FVS_uncertain2_3_hrf_data = beh_data_whole['FVS_uncertain2_3_hrf']
FVS_uncertain2_4_hrf_data = beh_data_whole['FVS_uncertain2_4_hrf']
FVS_uncertain2_5_hrf_data = beh_data_whole['FVS_uncertain2_5_hrf']

CVS_uncertain2_1_hrf_data = beh_data_whole['CVS_uncertain2_1_hrf']
CVS_uncertain2_2_hrf_data = beh_data_whole['CVS_uncertain2_2_hrf']
CVS_uncertain2_3_hrf_data = beh_data_whole['CVS_uncertain2_3_hrf']
CVS_uncertain2_4_hrf_data = beh_data_whole['CVS_uncertain2_4_hrf']
CVS_uncertain2_5_hrf_data = beh_data_whole['CVS_uncertain2_5_hrf']

AS_uncertain2_1_hrf_data = beh_data_whole['AS_uncertain2_1_hrf']
AS_uncertain2_2_hrf_data = beh_data_whole['AS_uncertain2_2_hrf']
AS_uncertain2_3_hrf_data = beh_data_whole['AS_uncertain2_3_hrf']
AS_uncertain2_4_hrf_data = beh_data_whole['AS_uncertain2_4_hrf']
AS_uncertain2_5_hrf_data = beh_data_whole['AS_uncertain2_5_hrf']

## neuro data
Nsubj = 54
ts_all = neuro_data_whole['ts_all_'+str(Nsubj)+'sub']
TR = 1.5
tp = int(np.shape(ts_all)[0]/Nsubj)
Nnodes = np.shape(ts_all)[1]
Nedges = int(Nnodes*(Nnodes-1)/2)
region_level = args.input[0]
region_level_label = dict({0 : 'node', 1: 'edge', 2: 'state', 3: 'system'})
print("GLM between original time series of each " + region_level_label[region_level] +  
      " and behavioral signal")
if region_level == 2:
    K = args.state[0]
    delta = args.state[1]
    iter_tmp = args.state[2]
    print('GLMM hyperparamters for states: ', K, delta, iter_tmp)
num_iter = args.iter[0]

# interpolation : tp(231) -> 3465(tp * 10 *TR)
x = np.linspace(0, tp-1, tp)  # 0 ~ 230
xq = np.linspace(x.min(), x.max(), int(tp * TR * 10))

if region_level == 0:
    neuro_data = np.zeros([Nsubj, tp, Nnodes])
    for si in range(Nsubj):
        tmp_ts = ts_all[0 + si*tp : tp + si*tp]
        neuro_data[si] = tmp_ts
    neuro_data = np.moveaxis(neuro_data, 0, 2)  # [tp x Nnode x Nsubj]
    X = neuro_data
    X_MinMax = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    neuro_data = X_MinMax
    tmp_neuro = interp1d(x, neuro_data, kind = 'quadratic', axis = 0)
    neuro_data_interp = tmp_neuro(xq)  
    Y_values = neuro_data_interp  # (3465, 147, 54)
    Nregions = Nnodes
elif region_level == 1:
    neuro_con_data_whole = scipy.io.loadmat(os.path.join(root_path, 'results', 'leida', 'Leida_eig.mat')) 
    Leading_Eig = neuro_con_data_whole['Leading_Eig']
    print(np.shape(Leading_Eig))
    LE_bold_all = np.reshape(Leading_Eig, [Nsubj, tp, Nnodes])
    del Leading_Eig
    mask = np.triu(np.ones((Nnodes, Nnodes), dtype=bool), k=1)
    i_mask = np.where(mask)
    neuro_con_data = np.zeros([Nsubj, tp, Nedges])
    for si in range(Nsubj):
        data = np.squeeze(LE_bold_all[si,:,:])
        for ti in range(tp):
            FC = np.dot(np.expand_dims(data[ti,:], axis = 1),np.expand_dims(data[ti,:], axis = 0))
            FC = np.nan_to_num(FC, copy=False)
            neuro_con_data[si,ti,:] = FC[i_mask[0],i_mask[1]]
    neuro_con_data = np.moveaxis(neuro_con_data, 0, 2)  # [tp x Nedges x Nsubj]
    tmp_neuro_con = interp1d(x, neuro_con_data, kind = 'quadratic', axis = 0)
    neuro_con_data_interp = tmp_neuro_con(xq)  
    v = squareform(np.linspace(0,Nedges-1,Nedges))
    
    tmp = scipy.io.loadmat(os.path.join(root_path, 'data', 'ROE_dopa.mat'))
    ROE_idx = tmp['mat']
    tmp = np.where(ROE_idx)
    ROE_ir = tmp[0]
    ROE_ic = tmp[1]
    ROE_IDX_dopa = list()
    for ii in range(len(ROE_ir)):
        ROE_IDX_dopa.append(int(v[ROE_ir[ii], ROE_ic[ii]]))
    neuro_con_data_interp_ROE_dopa = neuro_con_data_interp[:,ROE_IDX_dopa,:]

    tmp = scipy.io.loadmat(os.path.join(root_path, 'data', 'ROE_3system.mat'))
    ROE_idx = tmp['mat']
    tmp = np.where(ROE_idx)
    ROE_ir = tmp[0]
    ROE_ic = tmp[1]
    ROE_IDX_3system = list()
    for ii in range(len(ROE_ir)):
        ROE_IDX_3system.append(int(v[ROE_ir[ii], ROE_ic[ii]]))
    neuro_con_data_interp_ROE_3system = neuro_con_data_interp[:,ROE_IDX_3system,:]

    tmp = scipy.io.loadmat(os.path.join(root_path, 'data', 'ROE_FPN_sub_DAN.mat'))
    ROE_idx = tmp['mat']
    tmp = np.where(ROE_idx)
    ROE_ir = tmp[0]
    ROE_ic = tmp[1]
    ROE_IDX_FPN_sub_DAN = list()
    for ii in range(len(ROE_ir)):
        ROE_IDX_FPN_sub_DAN.append(int(v[ROE_ir[ii], ROE_ic[ii]]))
    neuro_con_data_interp_ROE_FPN_sub_DAN = neuro_con_data_interp[:,ROE_IDX_FPN_sub_DAN,:]

elif region_level == 2:
    glmm_data = scipy.io.loadmat(os.path.join(root_path, 'results','glmm','glmm_k'+str(K)+'_'+str(delta)+'.mat'))
    gamma_hats_iter = glmm_data['gamma_hats_iter']
    gamma_hats = np.squeeze(gamma_hats_iter[:,:,iter_tmp-1])
    gamma_hats_tmp = np.reshape(gamma_hats, [Nsubj, tp, K])
    gamma_hats_tmp = np.moveaxis(gamma_hats_tmp, 0, 2)  # [tp, K, Nsubj]
    tmp = interp1d(x, gamma_hats_tmp, kind = 'quadratic', axis = 0)
    gamma_hats_tmp_interp = tmp(xq)  
    Y_values = gamma_hats_tmp_interp  # (3465, 4, 54)
    Nregions = K
elif region_level == 3:
    system_con_data = scipy.io.loadmat(os.path.join(root_path, 'results', 'leida', 'Leida_eig.mat'))['inter_vec']
    Nregions = system_con_data.shape[1]  # 28
    system_con_data = np.reshape(system_con_data, [Nsubj, tp, Nregions])
    system_con_data = np.moveaxis(system_con_data, 0, 2)  # [tp x 28 x Nsubj]
    tmp_system_con = interp1d(x, system_con_data, kind = 'quadratic', axis = 0)
    system_con_data_interp = tmp_system_con(xq)  
    Y_values = system_con_data_interp # (3465, 28, 54)

DS1 = np.stack([cvs_hrf_data, fvs_hrf_data, as_hrf_data, prec_hrf_data], axis = 2)
DS2 = np.stack([cvs_hrf_data, fvs_l_hrf_data, fvs_w_hrf_data, as_hrf_data, prec_hrf_data], axis = 2)
DS3 = np.stack([cvs_hrf_data, ppe_hrf_data, npe_hrf_data, as_hrf_data, prec_hrf_data], axis = 2)

DS4 = np.stack([cvs_hrf_data, pw_hrf_data, pl_hrf_data, nw_hrf_data, nl_hrf_data, as_hrf_data, prec_hrf_data], axis = 2)
DS5 = np.stack([cvs_hrf_data, FVS_certain_hrf_data, FVS_uncertain_hrf_data, as_hrf_data, prec_hrf_data], axis = 2)
DS6 = np.stack([cvs_hrf_data, FVS_certain_hrf_data, FVS_uncertain1_hrf_data, FVS_uncertain2_hrf_data, as_hrf_data, prec_hrf_data], axis = 2)

DS7 = np.stack([cvs_hrf_data, FVS_certain_hrf_data, FVS_uncertain1_hrf_data, FVS_uncertain2_1_hrf_data, FVS_uncertain2_2_hrf_data, FVS_uncertain2_3_hrf_data, FVS_uncertain2_4_hrf_data, FVS_uncertain2_5_hrf_data, as_hrf_data, prec_hrf_data], axis = 2)
DS8 = np.stack([cvs_hrf_data, FVS_certain_hrf_data, FVS_uncertain_hrf_data, CVS_uncertain2_1_hrf_data, CVS_uncertain2_2_hrf_data, CVS_uncertain2_3_hrf_data, CVS_uncertain2_4_hrf_data, CVS_uncertain2_5_hrf_data, as_hrf_data, prec_hrf_data], axis = 2)
DS9 = np.stack([cvs_hrf_data, FVS_certain_hrf_data, FVS_uncertain_hrf_data, AS_uncertain2_1_hrf_data, AS_uncertain2_2_hrf_data, AS_uncertain2_3_hrf_data, AS_uncertain2_4_hrf_data, AS_uncertain2_5_hrf_data, as_hrf_data, prec_hrf_data], axis = 2)

num_m = int(list(globals().keys())[-1][2:])
MAX_RANDOM_SEED = 2**32 - 1

# mass GLM
print(num_iter, 'iteration will be started')
for di in range(1, num_m+1):
    EV_name = 'DS'+str(di)
    if EV_name in locals():
        print('DS:', di)
        if region_level == 1:
            if di == 1:
                Y_values = neuro_con_data_interp_ROE_dopa
            elif np.isin(di, [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]):
                Y_values = neuro_con_data_interp_ROE_3system
            elif np.isin(di, [42]):
                Y_values = neuro_con_data_interp_ROE_FPN_sub_DAN
            Nregions = Y_values.shape[1]
        X = globals()['DS{}'.format(di)]
        params_mass = np.zeros([Nregions, Nsubj, X.shape[2]+1, num_iter])
        std_error_mass = np.zeros([Nregions, Nsubj, X.shape[2]+1, num_iter])
        t_values_mass = np.zeros([Nregions, Nsubj, X.shape[2]+1, num_iter])
        p_values_mass = np.zeros([Nregions, Nsubj, X.shape[2]+1, num_iter])
        random_state=None
        X1 = X[:,114:1715,:]
        Y1 = Y_values[114:1715, :, :]
        X2 = X[:,1814:3415,:]
        Y2 = Y_values[1814:3415, :, :]
        if task_type == 0:
            X = np.hstack([X1, X2]) 
            Y = np.concatenate([Y1, Y2])
        elif task_type == 1:
            X = X1
            Y = Y1
        else:
            X = X2
            Y = Y2  
        for iter in range(num_iter):
            print('iter : ', iter)
            if isinstance(random_state, np.random.RandomState):
                prng=random_state
            else:
                prng = np.random.RandomState(random_state) 
            shifted_Y = phase_randomize(Y, voxelwise=True, random_state=prng) 
            for roi in range(Nregions):
                for si in range(Nsubj):
                    regr = linear_model.LinearRegression() # set up model
                    Xx = X[si,:,:]  # [3150, 6]
                    for vi in range(X.shape[2]):
                        x_tmp = Xx[:,vi]  # (3465,)
                        if np.max(x_tmp) == 0:
                            x_minmax = x_tmp
                        else:
                            x_minmax = (x_tmp - x_tmp.min()) / (x_tmp.max() - x_tmp.min())
                        Xx[:, vi] = x_minmax                           
                    shifted_Y_tmp = np.squeeze(shifted_Y[:,roi,si])  
                    regr.fit(Xx, shifted_Y_tmp) # fit model
                    params_tmp, std_error_tmp, t_values_tmp, p_values_tmp = get_statistics(regr,Xx,shifted_Y_tmp)
                    params_mass[roi, si, :, iter] = params_tmp
                    std_error_mass[roi, si, :, iter] = std_error_tmp
                    t_values_mass[roi, si, :, iter] = t_values_tmp
                    p_values_mass[roi, si, :, iter] = p_values_tmp
            random_state = np.random.RandomState(prng.randint(0, MAX_RANDOM_SEED))
        scipy.io.savemat(os.path.join(root_path,'results','GLM','null_act_'+region_level_label[region_level]+'_DS'+str(format(di, '02'))+'_'+task_label[task_type]+'.mat'), 
                         {'params_mass':params_mass, 'std_error_mass':std_error_mass, 't_values_mass':t_values_mass, 'p_values_mass':p_values_mass})
