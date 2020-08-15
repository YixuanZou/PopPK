"""All the functions.

Author Yixuan Zou
Version Alpha

"""
import os
import re
import operator
import subprocess
import collections
import shutil
import csv
import math
import time
import pandas as pd
import numpy as np
import scipy.stats as stats
from math import log10, floor

def bw_wald(thetas, cov, script, base_script, nobs,
            insert_variables, cwdir, covariates, data,
            quantile, quantile1, options='', version=73, bw_nm=True):
    """Main function for HWAM.
    thetas is the thetas vector from full model
    cov is the covariance matrix from full model
    script is NONMEM script with all covariates
    base_script is NONMEM script without any covariates
    nobs is the number of obserations from dataset
    insert_variables is all the covarate-parameter relationship
    cwdir is the current working directory
    covariates is all the covariates
    data is the NONMEM dataset
    quantile is the significance level for BE using approximation.
    quantile1 is the significance level for BE using NONMEM
    options is the option used for NONMEM run in the conmmand line1
    version is the NONMEM version used
    bw_nm is true means we need to run second confirmation in NONMEM
    """
    base_theta_num = count_theta(base_script)
    theta_num = count_theta(script)
    omega_num = count_omega(script)
    p = theta_num + omega_num - count_fixed(script)
    n = nobs
    non_select_list = []
    select_list = range(0, theta_num-base_theta_num)
    bw_list = []
    best_lrt = 0
    flag = True

    temp_final_select_list = []
    i = 0
    while i < (theta_num-base_theta_num):
        temp_sbc_list = []
        for x in select_list:
            temp_non_select_list = non_select_list + [x]
            q = len(temp_non_select_list)
            _, the_lrt = sbc_wam(thetas, cov, temp_non_select_list, p, q, n)
            temp_sbc_list.append((temp_non_select_list, x, the_lrt))
        temp_set, new_point, item_min = \
            min(temp_sbc_list, key=operator.itemgetter(2))

        if not flag:
            temp_final_select_list.append(new_point)

        if flag:
            if item_min - best_lrt > stats.chi2.ppf(quantile1, 1):
                flag = False
                temp_final_select_list.append(new_point)
        if flag:
            bw_list.append(([x+base_theta_num+1 for x in temp_set], item_min))

        best_lrt = item_min
        non_select_list = temp_set
        select_list = [x for x in select_list if x not in non_select_list]
        i += 1

    best_set = temp_final_select_list[:]
    update_list = [x+base_theta_num+1 for x in best_set]
    print 'Best set is ', update_list
    back_select_list = update_list
    with open('wam_interim_result.txt', 'r') as target:
        target.write(str(update_list))
    if bw_nm:
        _, _, back_select_list = covariate_model_bw(base_script, base_theta_num,
                                                    cwdir, data, nobs,
                                                    covariates,
                                                    quantile,
                                                    select_list=update_list,
                                                    options=options,
                                                    version=73)
    return bw_list, back_select_list

def covariate_model_bw(script, base_theta_num,
                       cwdir, data, nobs, covariates,
                       quantile, select_list=[], options='',
                       version=74):
    """Backward method in PK-PD. Do backward elimination.

    script is the NONMEM script without any covariates
    base_theta_num is the number of thetas in base model without any covariates
    cwdir is the current working directory
    nobs is the number of observation in the dataset
    data is the dataset
    nobs is the number of observation
    covariates are the potential covariates.
    quantile is the p value used to backward selection.
    select_list is the list that contains covariates pre-selected

    It returns the best model script and result.
    """
    # create template for full model
    base_script = 'base_model.ctl'
    full_script = 'full_script.ctl'
    _, _, _, insert_variables = \
        full_model(base_script, full_script, data, covariates)

    # initialize the best script and result file
    best_script = 'base_bw.ctl'
    best_result_file = 'base_bw.ext'
    best_result_path = 'base_bw.out'

    # calculate the number of thetas
    # create a full set of list
    num_theta = count_theta(full_script)
    full_set = range(1, num_theta-base_theta_num+1)

    # create update list for those covariates already selected
    update_list = [x+base_theta_num-1 for x in full_set
                   if x+base_theta_num not in select_list]

    # update the template according to the covariates selected
    model_update(full_script, best_script, update_list,
                 pattern='0 FIXED\n', var='THETA')

    # run the model to get base OFV
    run_model(best_script, best_result_path, options, version=version)

    # only take the last part of results
    best_result_file = sep_ext(best_result_file)
    select = [x-base_theta_num for x in select_list]

    # backward elimination begins
    step = 1
    while True:
        print 'Backward elimination step %d.' % step
        result_file_list = []
        result_path_list = []
        script_list = []
        for i, index in enumerate(select):
            print 'Deleting Theta', index + base_theta_num
            result_path_i = 'model_covariate_del_%d_%d.out' % (step, i)
            result_file_i = 'model_covariate_del_%d_%d.ext' % (step, i)
            control_file_i = 'model_covariate_del_%d_%d.ctl'\
                % (step, i)

            model_update(best_script, control_file_i,
                         [index+base_theta_num-1],
                         pattern='0 FIXED\n', var='THETA')

            run_model(control_file_i, result_path_i, options, version=version)
            result_file_i = sep_ext(result_file_i)
            script_list.append(control_file_i)
            result_path_list.append(result_path_i)
            result_file_list.append(result_file_i)
            obj = read_obj(result_file_i)
            if obj == 999999999:
                print 'Deleting Theta', index + base_theta_num
                shutil.copy(control_file_i, 'crash_'+control_file_i)
            else:
                print 'Objective function of step %d is %.4f' % (step, obj)

        worst_index = best_model_covariate(result_file_list)
        temp_best_result_file = 'final_model_covariate_del_%d_%d.ext' % \
            (step, worst_index)
        temp_best_script = 'model_covariate_del_%d_%d.ctl' % \
            (step, worst_index)
        dof = degree_free(temp_best_script, best_script)
        b_ss = model_compare(temp_best_result_file, best_result_file,
                             quantile, dof)
        if b_ss:
            print 'There is nothing to delete for step %s.' % step
            break

        best_script = temp_best_script
        best_result_file = temp_best_result_file
        obj = read_obj(best_result_file)
        select = [x for x in select if x != select[worst_index]]

        if len(select) == 0:
            best_script = base_script
            best_result_file = base_script.replace(".ctl", ".ext")
            print 'The best model is the base model.'
            break

        model_update_initial(best_script, best_result_file, base_theta_num)
        step += 1

    shutil.copy(best_script, 'best_bw_' + best_script)
    best_select = range(1, 1+base_theta_num) + [x+base_theta_num for x in select]
    with open('bw_set.txt', 'w') as target:
        target.write(str(best_select))

    return best_script, best_result_file, best_select

def model_update_initial(script, old_result_file, base_theta_num, nsig=4):
    """Update the model by using the estimates from previous model."""
    old_out = read_result_file(old_result_file)
    estimates = old_out[old_out.ITERATION == -1000000000].iloc[:, 1:-1]

    thetas = collections.OrderedDict()

    for name in estimates:
        if 'THETA' in name:
            thetas[name] = estimates.loc[:, name].values[0]
    # pull out the estimates from the data

    theta_num = len(thetas)
    theta_value = thetas.values()

    initial_new_model(script, theta_value, 'THETA',
                      theta_num, '(1E-10, %.{0}f)\n'.format(nsig),
                      '(-10, %.{0}f, 10)\n'.format(nsig), base_theta_num)

def model_compare(current_result_file, old_result_file, quantile, dof):
    """Choose the model by chi-square test. True means current is better."""
    current_obj = read_obj(current_result_file)
    old_obj = read_obj(old_result_file)
    crit = stats.chi2.ppf(quantile, abs(dof))
    if crit < abs(old_obj - current_obj):
        return True
    else:
        return False


def degree_free(current_script, old_script):
    """Find the degree of freedom of hypothesis testing."""
    current_npar = count_par(current_script)
    old_npar = count_par(old_script)
    dof = current_npar - old_npar
    return dof

def count_par(script):
    """Count how many parameters in the model."""
    npar = 0
    theta_num = 0
    in_theta = False
    eta_num = 0
    in_eta = False
    with open(script, 'r') as target:
        for line in target:
            if line.startswith('CLP'):
                in_theta = True
            if in_theta and line.startswith('CL='):
                in_theta = False
            if in_theta and 'THETA' in line:
                theta_num += 1
            if line.startswith('CL='):
                in_eta = True
            if in_eta and line.startswith('S1='):
                in_eta = False
            if in_eta and '*EXP' in line:
                eta_num += 1
    npar = eta_num + theta_num - count_fixed(script) + 2
    return npar

def best_model_covariate(result_file_list):
    """Select the best model when adding one covariate.

    You should input the a list of results.
    It returns to covariate and index of the best model.
    """
    objs = []
    for result_file in result_file_list:
        mod_index = int(re.findall(r'\d+', result_file)[-1])
        obj = read_obj(result_file)
        objs.append((mod_index, obj))

    best_index = min(objs, key=operator.itemgetter(1))[0]
    return best_index

def sbc_wam(thetas, cov, the_set, p, q, n):
    """Calculate sbc and lrt of the test."""
    thetas_redu = thetas[:, the_set]
    cov_redu = cov[np.ix_(the_set, the_set)]
    if len(the_set) > 1:
        lrt = thetas_redu.dot(np.linalg.inv(cov_redu)
                              ).dot(thetas_redu.transpose())[0, 0]
    elif len(the_set) == 1:
        lrt = (thetas_redu * (1/cov_redu) * thetas_redu)[0][0]
    else:
        lrt = 0
    sbc = -lrt - (p-q)*math.log(n)
    return sbc, lrt

def count_fixed(script):
    """Count how many parameter fixed in control file."""
    count = 0
    with open(script, 'r') as target:
        for line in target:
            if 'FIXED' in line:
                count += 1
    return count

def read_obj(result_file, sig=10):
    """Read the obj of the result file."""
    try:
        out = read_result_file(result_file)
        if 'ITERATION' in out.columns:
            obj = out[out.ITERATION == -1000000000].loc[:, 'OBJ'].values[0]
        else:
            obj = 999999999
        return obj
    except IOError:
        print 'This model does not exist.'

def count_theta(script):
    """Count theta for general NONMEM script."""
    theta_num = 0
    with open(script, 'r') as target:
        for line in target:
            if 'THETA(' in line:
                theta_num += 1
    return theta_num


def count_tv(script):
    """Count typical value for NONMEM script."""
    tv_num = 0
    with open(script, 'r') as target:
        for line in target:
            if line.startswith('TV'):
                tv_num += 1
    return tv_num


def sep_ext(result_file):
    """Seperate extension file for foce and important sampling."""
    output_lines = []
    with open(result_file, 'r') as target:
        table_num = 0
        for line in target:
            if line.startswith('TABLE'):
                table_num += 1

    with open(result_file, 'r') as target:
        i = 0
        for line in target:
            if line.startswith('TABLE'):
                i += 1
            if i == table_num:
                output_lines.append(line)
    ext = 'final_' + result_file
    with open(ext, 'w') as target:
        for line in output_lines:
            target.write(line)
    return ext


def read_theta(result_file, tv_num, theta_num):
    """Read the theta estimates from the result_file."""
    data_frame = read_result_file(result_file)
    estimates = data_frame[data_frame.ITERATION == -1000000000]
    thetas = estimates.iloc[:, (tv_num+1):(theta_num+1)]
    return thetas.values


def read_cov(cov_file, tv_num, theta_num):
    """Read the covariance matrix from .cov file."""
    data_frame = pd.read_csv(cov_file,
                             sep=' ',
                             skipinitialspace=True,
                             skiprows=1)
    del data_frame['NAME']
    cov = data_frame.iloc[tv_num:theta_num, tv_num:theta_num]
    return cov.values

def run_model(control_file, result_path_file,
              options='', version=73, nobuild=False,
              intel=False):
    """Run the model in the shell."""
    if not nobuild:
        cmd_script =\
            'nmfe%d.bat %s %s ' % (version, control_file, result_path_file)
        cmd_script = cmd_script + options
        redirect = r'CALL "C:\Program Files (x86)'
        redirect = redirect +\
            r'\IntelSWTools\compilers_and_libraries_2018.1.156'
        redirect = redirect +\
            r'\windows\bin\ipsxe-comp-vars.bat" intel64 vs2017'
        batch_file = 'run_nonmem.bat'
        with open(batch_file, 'w') as target:
            target.write(redirect)
            target.write('\n')
            target.write(cmd_script)
        try:
            if intel:
                subprocess.call(batch_file)
            else:
                subprocess.call(cmd_script)
        except Exception:
            print 'There is an error in your script.'
    else:
        cmd_script =\
         'nmfe%d.bat %s %s -nobuild' %\
         (version, control_file, result_path_file)
        cmd_script = cmd_script + options
        redirect = r'CALL "C:\Program Files (x86)'
        redirect = redirect +\
            r'\IntelSWTools\compilers_and_libraries_2018.1.156'
        redirect = redirect +\
            r'\windows\bin\ipsxe-comp-vars.bat" intel64 vs2017'
        batch_file = 'run_nonmem.bat'
        with open(batch_file, 'w') as target:
            target.write(redirect)
            target.write('\n')
            target.write(cmd_script)
        try:
            if intel:
                subprocess.call(batch_file)
            else:
                subprocess.call(cmd_script)
        except Exception:
            print 'There is an error in your script.'

def full_model(base_script, script, data, covariates, IMPMAP=False):
    """Construct the full model for WAM on a generalized base."""
    # create related files for full model script
    result_path = script.replace('ctl', 'out')
    result_file = script.replace('ctl', 'ext')
    cov_file = script.replace('ctl', 'cov')

    # find all the parameter with non fixed omega in the model
    tvs_all = []
    with open(base_script, 'r') as target:
        for line in target:
            if line.startswith('TV'):
                tvs_all.append(line.split('=')[0].replace('TV', ''))
    tvs = [tvs_all[i] for i in omega_nonzero(base_script)]

    # truncate all covariates
    trun_covariates = truncate_cov(covariates)

    # number of omega
    num_omega = len(tvs)

    # construct a dictionary for the parameter covariate relationship
    insert_variables = dict()
    num = 1
    for tv in tvs:
        for trun_covariate in trun_covariates:
            insert_variables[num] = \
                (tv, trun_covariate,
                 '%s%s' % (tv, trun_covariate),
                 trun_covariates[trun_covariate],
                 covariates[trun_covariates[trun_covariate]])
# pk_parameter, truncated covariate, insert_variable, covariate. variable_type
            num += 1

    # copy base_script to script
    shutil.copy(base_script, script)

    output_lines = []
    with open(base_script, 'r') as target:
        for line in target:
            if line.startswith('$ERROR'):
                if IMPMAP:
                    for i in range(1, num_omega+1):
                        output_lines.append('MU_%d=LOG(TV%s)\n' %
                                            (i, tvs[i-1]))

            if line.startswith('$ESTIMATION'):
                line1 = '$ESTIMATION METHOD=1 INTERACTION' + \
                    ' MAXEVAL=9999 PRINT=10 NOABORT\n'
                output_lines.append(line1)
                if IMPMAP:
                    line2 = '$ESTIMATION METHOD=IMPMAP INTERACTION' +\
                        ' EONLY=1 ISAMPLE=10000 NITER=5 SEED=12345 PRINT=1\n'
                    output_lines.append(line2)
                continue
            output_lines.append(line)

    with open(script, 'w') as target:
        for line in output_lines:
            target.write(line)

    # extract the individual demographic information
    demo_df = data.groupby('CID').first()

    # add covariate parameter one by one
    for index in range(1, len(insert_variables)+1):
        variable_type = insert_variables[index][4]
        covariate = insert_variables[index][3]
        tv = insert_variables[index][0]
        print 'add', insert_variables[index][2]
        if variable_type != 'factor':
            covariate_median = np.median(demo_df[covariate])
            add_num_covariate(index, insert_variables,
                              covariate_median, script, script)
        else:
            level = len(demo_df[covariate].unique())
            add_factor_covariate(index, insert_variables, level,
                                 script, script)
    return result_path, cov_file, result_file, insert_variables

def add_num_covariate(index, insert_variables,
                      median, best_script, new_script):
    """Add numerical covariate to the best structural model."""
    tv = insert_variables[index][0]
    insert_variable = insert_variables[index][2]
    covariate = insert_variables[index][3]
    theta_num = 0
    with open(best_script, 'r') as script:
        for line in script:
            if 'THETA(' in line:
                theta_num += 1

    # add continuous covariates
    output_lines = []
    with open(best_script, 'r') as script:
        for line in script:
            if line.startswith('$PRED') or  line.startswith('$PK'):
                output_lines.append(line)
                output_lines.append('%s=(%s/%.8f)**THETA(%d)\n' %
                                    (insert_variable,
                                     covariate,
                                     median,
                                     theta_num+1))
                continue
            if line.startswith('TV' + tv + '='):
                output_lines.append(line.strip('\n') +
                                    '*%s\n' % insert_variable)
                continue
            if line.startswith('$OMEGA'):
                output_lines.append('(-10, 0.1, 10)\n')
            output_lines.append(line)

    with open(new_script, 'w') as script:
        for line in output_lines:
            script.write(line)


def add_factor_covariate(index, insert_variables,
                         level, best_script, new_script):
    """Add factor covariate to the best structural model."""
    tv = insert_variables[index][0]
    insert_variable = insert_variables[index][2]
    covariate = insert_variables[index][3]
    theta_num = 0
    # count theta number
    with open(best_script, 'r') as script:
        for line in script:
            if 'THETA(' in line:
                theta_num += 1

    # add factor covariates
    output_lines = []
    with open(best_script, 'r') as script:
        for line in script:
            if line.startswith('$PRED') or  line.startswith('$PK'):
                output_lines.append(line)
                output_lines.append('%s=1\n' % insert_variable)
                for j in range(1, level):
                    output_lines.append('IF(%s==%d) %s=EXP(THETA(%d))\n' %
                                        (covariate, j,
                                         insert_variable, theta_num+j))
                continue
            if line.startswith('TV' + tv + '='):
                output_lines.append(line.strip('\n') +
                                    '*%s\n' % insert_variable)
                continue
            if line.startswith('$OMEGA'):
                for j in range(1, level):
                    output_lines.append('(-10, 0.1, 10)\n')
            output_lines.append(line)

    with open(new_script, 'w') as target:
        for line in output_lines:
            target.write(line)

def read_result_file(result_file):
    """Read the result file into python as data frame."""
    try:
        data_frame = pd.read_csv(result_file,
                                 sep=' ',
                                 skipinitialspace=True,
                                 skiprows=1)
    except IOError:
        data_frame = pd.DataFrame()
        print '%s does not exist.' % result_file
    return data_frame

def round_sig(x, sig=2):
    """Round the value"""
    return round(x, sig-int(floor(log10(abs(x))))-1)

def read_obj(result_file, sig=10):
    """Read the obj of the result file."""
    try:
        out = read_result_file(result_file)
        if 'ITERATION' in out.columns:
            obj = out[out.ITERATION == -1000000000].loc[:, 'OBJ'].values[0]
        else:
            obj = 999999999
        return obj
    except IOError:
        print 'This model does not exist.'


def omega_nonzero(list_file):
    """Find the non-zero omegas and return them to a list."""
    with open(list_file, 'r') as target:
        omega_list = []
        count = 0
        in_omega = False
        for line in target:
            if line.startswith('$SIGMA'):
                in_omega = False
            if in_omega and 'FIXED' not in line and line[0].isdigit():
                omega_list.append(count)
            if in_omega:
                count += 1
            if line.startswith('$OMEGA'):
                in_omega = True
    return omega_list

def truncate_cov(covariates):
    """Truncate the covariates to decrease length of code.

    covariates is dictionary containing variable name and type.
    We return a dictionary of trun_covariate and its full name.
    """
    trun_covariates = dict()
    for covariate in covariates:
        i = 4
        flag = True
        while flag:
            flag = (covariate[0:i] in trun_covariates)
            if flag is False:
                trun_covariates[covariate[0:i]] = covariate
            i += 1
    return trun_covariates

def create_folder(path, folder_name):
    """Create folder for the scripts under current directory."""
    path = path + '\\' + folder_name
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    return path

def num_obs(data):
    """Find how many patients in the data."""
    count = sum(x == 0 for x in data.EVID.values)
    return count

def count_omega(script):
    """Count number of OMEGA."""
    omega_num = 0
    with open(script, 'r') as target:
        for line in target:
            if 'ETA(' in target:
                omega_num += 1
    return omega_num

def model_update(old_file, new_file, update_list,
                 pattern='0 FIXED\n', var='OMEGA'):
    """Update model according to certain rule."""
    output_lines = []
    in_para = False
    para_num = 0
    with open(old_file, 'r') as target:
        for line in target:
            if in_para and line.startswith('$'):
                in_para = False

            if in_para:
                if para_num in update_list:
                    output_lines.append(pattern)
                else:
                    output_lines.append(line)
                para_num += 1
                continue

            if line.startswith('$'+var):
                in_para = True
            output_lines.append(line)

    with open(new_file, 'w') as target:
        for line in output_lines:
            target.write(line)

def move_files(source, pattern, dest):
    """Move the scripts to the folder."""
    for the_file in os.listdir(dest):
        file_path = os.path.join(dest, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as err:
            print err

    for the_file in os.listdir(source):
        if re.search(pattern, the_file):
            shutil.move(the_file, dest)

def initial_new_model(script, values, parameter,
                      theta_num, pattern1='%.8f\n',
                      pattern2='%.8f\n', base_theta_num=0):
    """Insert old values into the new model."""
    output_lines = []
    in_para = False
    num = 0
    with open(script, 'r') as target:
        for line in target:
            if in_para and line.startswith('$'):
                in_para = False
            if in_para:
                if parameter == 'THETA':
                    if num in range(0, theta_num):
                        if num <= (base_theta_num- 1):
                            if values[num] == 0:
                                output_lines.append('0 FIXED \n')
                            elif values[num] == 1:
                                output_lines.append('1 FIXED \n')
                            else:
                                output_lines.append(pattern1 % values[num])
                        else:
                            if values[num] == 0:
                                output_lines.append('0 FIXED \n')
                            elif values[num] == 1:
                                output_lines.append('1 FIXED \n')
                            else:
                                output_lines.append(pattern2 % values[num])
                    else:
                        output_lines.append(line)
                    num += 1
                    continue
                else:
                    if num in range(0, theta_num):
                        if values[num] == 0:
                            output_lines.append('0 FIXED \n')
                        elif values[num] == 1:
                            output_lines.append('1 FIXED \n')
                        else:
                            output_lines.append(pattern1 % values[num])
                    else:
                        output_lines.append(line)
                    num += 1
                    continue

            if line.startswith('$' + parameter):
                in_para = True

            output_lines.append(line)

    with open(script, 'w') as target:
        for line in output_lines:
            target.write(line)
