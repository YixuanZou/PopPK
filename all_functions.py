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
from random import sample

def num_obs(data):
    """Find how many patients in the data."""
    count = sum(x == 0 for x in data.EVID.values)
    return count


def select_structural_model(cwdir, data_path, nobs, variable_name,
                            the_input, quantile):
    """Select the best structural model."""
    best = 1
    control_file = 'model_%d.ctl'
    result_path_file = 'model_%d.out'
    result_file = 'model_%d.ext'

    for i in range(1, 4):
        structural_model(i, variable_name, data_path, the_input)

    i = 1
    j = 1
    ss_list = []
    while i < 4:
        step = 1
        if j == i:
            print 'Running %d compartment model.' % i
        else:
            print 'Running updated %d compartment model.' % i

        result_path_file_i = result_path_file % i
        result_file_i = result_file % i
        control_file_i = control_file % i
        run_model(control_file_i, result_path_file_i)
        obj = read_obj(result_file_i)
        if obj == 999999999:
            print '%d compartment model crashes.' % i
            shutil.copy(control_file_i, 'crash_'+control_file_i)
            break

        converge, cv_error, min_omega, _ = model_converge(result_file_i)
        if i > 1:
            old_result_file = result_file % (i - 1)
            old_control_file = control_file % (i - 1)
            dof = degree_free(control_file_i, old_control_file)
            ss = model_compare(result_file_i, old_result_file, quantile, dof)
            if ss:
                if converge:
                    ss_list.append('YES')
                    best = i
                    step = 1
                else:
                    if cv_error:
                        print cv_error, \
                            ': coefficient variation is' +\
                            ' bigger than 50%, fix it to 0 and run it.'
                        _, name_list = find_index(cv_error)
                        for index, var_name in name_list:
                            model_update(control_file_i, control_file_i,
                                         [index],
                                         pattern='0 FIXED\n', var=var_name)
                        step = 0
                    elif min_omega:
                        print min_omega,\
                              'is the minimal omega.' + \
                              ' When not converging, fix it to 0 and run it.'
                        update_list, _ = find_index(min_omega)
                        model_update(control_file_i, control_file_i,
                                     update_list)
                        step = 0
                    else:
                        print '%i compartment model does not converge.' % i
                        shutil.copy(control_file_i,
                                    'unconverge_'+control_file_i)
                        break
            else:
                ss_list.append('NO')
                print '%d compartment model is not improved.' % i
                break
        else:
            ss_list.append('NONE')
            if not converge:
                if cv_error:
                    print cv_error, \
                        ': coefficient variation is' +\
                        ' bigger than 50%, fix it to 0 and run it.'
                    _, update_list = find_index(cv_error)
                    for index, var_name in update_list:
                        model_update(control_file_i, control_file_i,
                                     [index],
                                     pattern='0 FIXED\n', var=var_name)
                    step = 0
                elif min_omega:
                    print min_omega,\
                          'is the minimal omega.' + \
                          ' When not converging, fix it to 0 and run it.'
                    update_list, _ = find_index(min_omega)
                    model_update(control_file_i, control_file_i, update_list)
                    step = 0
                else:
                    print '%i compartment model does not converge.' % i
                    shutil.copy(control_file_i, 'unconverge_'+control_file_i)
                    break

        if (step == 1) & (i < 3):
            model_update_initial(control_file % (i+1), result_file_i, i)
        i += step
        if step == 1:
            j = i
            print 'The objective function of %d compartment:' % (i-1), obj
        else:
            j += 1

    if not converge and i == 1:
        print 'Even one compartment model crashes.'
    else:
        best_script = control_file % best
        best_result = result_file % best
        obj = read_obj(best_result)
        check_sigma(best_script, best_result, quantile, best)
        print '%d compartment model is the best base model.' % best
        print 'The objective function of %d compartment:' % best, obj
        base_report(best, ss_list, cwdir, nobs, quantile, best)
        return best, 'model_%d.ctl' % best, 'model_%d.ext' % best


def check_sigma(best_script, best_result, quantile, ncomp):
    """Check if we need to remove one sigma."""
    control_file_list = []
    result_file_list = []
    for i in range(0, 2):
        control_file = 'model_err%d.ctl' % i
        model_update(best_script, control_file, [i], var='SIGMA')
        control_file_list.append(control_file)
        result_file = 'model_err%d.ext' % i
        result_file_list.append(result_file)
        result_file_path = 'model_err%d.out' % i
        run_model(control_file, result_file_path)

    index = best_model_covariate(result_file_list)
    current_result = 'model_err%d.ext' % index
    current_script = 'model_err%d.ctl' % index
    if model_compare(current_result, best_result, quantile, 1):
        converge, _, _, _ = model_converge(result_file_list[index])
        if converge:
            shutil.copy(current_script, best_script)


def base_report(best, ss_list, cwdir, nobs, quantile, ncomp):
    """Generate report after running base model."""
    control_file = 'model_%d.ctl'
    result_file = 'model_%d.ext'
    result_table = []
    title = ['Base model development']
    colnames = ['MODEL', 'NPAR', 'OBJ', 'AIC', 'BIC',
                'MINIMIZATION', 'COVARIANCE', 'SELECT', 'SIGNIFICANT']
    conclusion = ['After hypertheis testing with p=%f,' % (1-quantile) +
                  ' the best base model is %d compartment model.' % best]
    contents = []
    for i in range(1, best+1):
        control_file_i = control_file % i
        result_file_i = result_file % i
        xml_file = 'model_%d.xml' % i
        if control_file_i in os.listdir(cwdir):
            num_fixed = count_fixed(control_file_i)
            npar = 2 + 4*i - num_fixed
            obj = read_obj(result_file_i)
            aic = obj + 2*npar
            bic = obj + npar*math.log(nobs)
            minimize = min_status(xml_file)
            converge, _, _, _ = model_converge(result_file_i)
            select = (i == best)
            ss = ss_list[i-1]
            obs = []
            obs.extend([i, npar, obj, aic, bic,
                        minimize, converge, select, ss])
            contents.append(obs)

    result_table.append(title)
    result_table.append(colnames)
    result_table.extend(contents)
    result_table.append(conclusion)
    with open('final_base_result.csv', 'wb') as csvfile:
        writer = csv.writer(csvfile)
        [writer.writerow(r) for r in result_table]


def min_status(xml_file):
    """Check if obv reaches its minimam. 0 means minimazation successful."""
    status = False
    with open(xml_file, 'r') as target:
        for line in target:
            if 'termination_status' in line:
                status_ind = int(re.findall('\d+', line)[0])
                if status_ind == 0:
                    status = True
    return status


def count_fixed(script):
    """Count how many parameter fixed in control file."""
    count = 0
    with open(script, 'r') as target:
        for line in target:
            if 'FIXED' in line:
                count += 1
    return count


def run_model(control_file, result_path_file,
              options='', version=74, nobuild=False,
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
        try:
            if intel:
                with open(batch_file, 'w') as target:
                    target.write(redirect)
                    target.write('\n')
                    target.write(cmd_script)
                subprocess.call(batch_file)
            else:
                subprocess.call(cmd_script, shell=True)
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
        try:
            if intel:
                with open(batch_file, 'w') as target:
                    target.write(redirect)
                    target.write('\n')
                    target.write(cmd_script)
                subprocess.call(batch_file)
            else:
                subprocess.call(cmd_script, shell=True)
        except Exception:
            print 'There is an error in your script.'


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


def mcmc_model_update(script, sample_size):
    """Update NONMEM script for Monte Carlo Method."""
    with open(script, 'r') as target:
        for line in target:
            if line.startswith('$EST'):
                index = int(re.findall('\d+', line).pop())
    output_lines = []
    with open(script, 'r') as target:
        if index == 0:
            skip = False
            for line in target:
                if line.startswith('CTYPE'):
                    continue
                if line.startswith('$THETA'):
                    skip = True
                if line.startswith('$EST'):
                    skip = False
                if skip:
                    continue
                if line.startswith('$EST'):
                    new_line = '$MSFI ' + 'M' + str(index) + '.msf' + '\n'
                    output_lines.append(new_line)
                    line = line.replace('1000', str(sample_size))
                    line = line.replace('50', '15')
                    line = line.replace('M' + str(index), 'M' + str(index + 1))
                    line = line + 'CTYPE=2 CITER=5 CALPHA=0.05\n'
                output_lines.append(line)
        else:
            skip = False
            for line in target:
                if line.startswith('CTYPE'):
                    continue
                if line.startswith('$THETA'):
                    skip = True
                if line.startswith('$MSFI'):
                    skip = False
                if skip:
                    continue
                if line.startswith('$MSFI'):
                    continue
                if line.startswith('$EST'):
                    new_line = '$MSFI ' + 'M' + str(index) + '.msf' + '\n'
                    output_lines.append(new_line)
                    old_sample_size = re.findall('\d+', line)[0]
                    line = line.replace(old_sample_size, str(sample_size))
                    line = line.replace('M' + str(index), 'M' + str(index + 1))
                    line = line + 'CTYPE=2 CITER=5 CALPHA=0.05\n'
                output_lines.append(line)
    with open(script, 'w') as target:
        for line in output_lines:
            target.writelines(line)


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


def find_index(str_list):
    """Extract the index for those parameters make not converge."""
    index_list = []
    name_list = []
    for str in str_list:
        point = re.findall(r'\d+', str)
        index_list.append(int(point[0])-1)
        name = re.findall('[A-Z]+', str)
        name_list.append((int(point[0])-1, name[0]))
    return list(set(index_list)), name_list


def structural_model(i, variable_name, data_path, the_input):
    """Write structural model script."""
    file_name = 'model_%d' % i
    model = '%s.ctl' % file_name
    with open(model, 'w') as target:
        target.truncate()
        target.write('$PROBLEM %d-COMPARTMENT MODEL\n' % i)
        target.write(the_input)
        target.write('$DATA' + ' ' + data_path + ' IGNORE=C\n')
        if i == 1:  # one compartment model
            target.write('$SUBROUTINES  ADVAN1 TRANS2\n')
            target.write('$PK\n')
            target.write('CLP=THETA(1)\n')
            target.write('VP=THETA(2)\n')
            target.write('CLT=CLP\n')
            target.write('VT=VP\n')
            target.write('CL=CLT*EXP(ETA(1))\n')
            target.write('V=VT*EXP(ETA(2))\n')
            target.write('S1=V\n')

        if i == 2:
            target.write('$SUBROUTINES  ADVAN3 TRANS4\n')
            target.write('$PK\n')
            target.write('CLP=THETA(1)\n')
            target.write('V1P=THETA(2)\n')
            target.write('QP=THETA(3)\n')
            target.write('V2P=THETA(4)\n')
            target.write('CLT=CLP\n')
            target.write('V1T=V1P\n')
            target.write('QT=QP\n')
            target.write('V2T=V2P\n')
            target.write('CL=CLT*EXP(ETA(1))\n')
            target.write('V1=V1T*EXP(ETA(2))\n')
            target.write('Q=QT*EXP(ETA(3))\n')
            target.write('V2=V2T*EXP(ETA(4))\n')
            target.write('S1=V1\n')

        if i == 3:
            target.write('$SUBROUTINES ADVAN11 TRANS4\n')
            target.write('$PK\n')
            target.write('CLP=THETA(1)\n')
            target.write('V1P=THETA(2)\n')
            target.write('Q2P=THETA(3)\n')
            target.write('V2P=THETA(4)\n')
            target.write('Q3P=THETA(5)\n')
            target.write('V3P=THETA(6)\n')
            target.write('CLT=CLP\n')
            target.write('V1T=V1P\n')
            target.write('Q2T=Q2P\n')
            target.write('V2T=V2P\n')
            target.write('Q3T=Q3P\n')
            target.write('V3T=V3P\n')
            target.write('CL=CLT*EXP(ETA(1))\n')
            target.write('V1=V1T*EXP(ETA(2))\n')
            target.write('Q2=Q2T*EXP(ETA(3))\n')
            target.write('V2=V2T*EXP(ETA(4))\n')
            target.write('Q3=Q3T*EXP(ETA(5))\n')
            target.write('V3=V3T*EXP(ETA(6))\n')
            target.write('S1=V1\n')

        target.write('$ERROR\n')
        target.write('IPRED=F\n')
        target.write('Y=F*(1+ERR(1))+ERR(2)\n')
        target.write('$THETA\n' + '(0.001, 1)\n' * 2 * i)
        target.write('$OMEGA\n' + '0.1\n' * 2 * i)
        target.write('$SIGMA\n0.1\n0.1\n')
        target.write('$ESTIMATION METHOD=1 INTERACTION ' +
                     'MAXEVAL=9999 PRINT=5 NOABORT\n')
        target.write('$COV UNCONDITIONAL\n')
        target.write('$TABLE' + ' ' + variable_name +
                     ' ' + 'IPRED PRED CWRES\n')
        target.write('FILE=%s.CSV  NOAPPEND NOPRINT' % file_name.upper())


def model_converge(result_file):
    """Test whether the model converges."""
    cv_error = []
    min_omega = []
    converge = True
    out = read_result_file(result_file)
    xml_file = result_file.replace('ext', 'xml')

    if 'ITERATION' in out.columns:
        if -1000000001 in out.ITERATION.values:
            if -1000000007 in out.ITERATION.values:
                code = out[out.ITERATION == -1000000007].iloc[:, 1].values[0]
                for parameter in list(out.columns.values)[1:-1]:
                    std = out[out.ITERATION == -1000000001].loc[
                        :, parameter].values[0]
                    estimate = out[out.ITERATION == -1000000000].loc[
                        :, parameter].values[0]
                    cv = std / estimate  # coefficient of variation
                    if estimate not in [0, 1]:
                        index = int(re.findall(r'\d+', parameter)[0])
                        if ('OMEGA' in parameter or 'SIGMA' in parameter or
                           ('THETA' in parameter)):
                            if cv > 0.5 or estimate < 0.0001:
                                converge = False
                                cv_error.append(parameter)
                                break
                        else:
                            if cv > 0.5:
                                converge = False

            else:
                status = min_status(xml_file)
                if status:
                    for parameter in list(out.columns.values)[1:-1]:
                        std = out[out.ITERATION == -1000000001].loc[
                            :, parameter].values[0]
                        estimate = out[out.ITERATION == -1000000000].loc[
                            :, parameter].values[0]
                        cv = std / estimate  # coefficient of variation
                        if estimate not in [0, 1]:
                            index = int(re.findall(r'\d+', parameter)[0])
                            if ('OMEGA' in parameter or 'SIGMA' in parameter or
                               ('THETA' in parameter)):
                                    if cv > 0.5 or estimate < 0.0001:
                                        converge = False
                                        cv_error.append(parameter)
                                        break
                            else:
                                if cv > 0.5:
                                    converge = False
                else:
                    converge = False
        elif -1000000000 in out.ITERATION.values:
            if -1000000007 in out.ITERATION.values:
                code = out[out.ITERATION == -1000000007].iloc[:, 1].values[0]
            converge = False
            estimates = out[out.ITERATION == -1000000000]
            omegas = dict()
            for name in estimates:
                if 'OMEGA' in name:
                    omegas[name] = estimates.loc[:, name].values[0]
            omega_Diagnogal = dict()
            for omega in omegas:
                if omegas[omega] != 0:
                    omega_Diagnogal[omega] = omegas[omega]
            min_omega = []
            min_omega.append(min(omega_Diagnogal.iteritems(),
                                 key=operator.itemgetter(1))[0])
        else:
            converge = False
    else:
        converge = False

    if -1000000007 in out.ITERATION.values:
        return converge, cv_error, min_omega, code
    else:
        return converge, cv_error, min_omega, None


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


def count_theta(script):
    """Count theta for general NONMEM script."""
    theta_num = 0
    with open(script, 'r') as target:
        for line in target:
            if 'THETA(' in line:
                theta_num += 1
    return theta_num


def round_sig(x, sig=2):
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


def worst_model_covariate(result_file_list):
    """Select the worst model when adding one covariate.

    You should input the a list of results.
    It returns to covariate and index of the worst model.
    """
    objs = []
    for result_file in result_file_list:
        mod_index = int(re.findall(r'\d+', result_file)[-1])
        obj = read_obj(result_file)
        objs.append((mod_index, obj))

    worst_index = max(objs, key=operator.itemgetter(1))[0]
    return worst_index


def pk_para(ncomp):
    """Generate PK parameter in NONMEM."""
    temp_list = []
    if ncomp == 1:
        temp_list = ['CL', 'V']
    elif ncomp == 2:
        temp_list = ['CL', 'V1', 'Q', 'V2']
    else:
        temp_list = ['CL', 'V1', 'Q2', 'V2', 'Q3', 'V3']
    return temp_list


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
                output_lines.append(line.rstrip() +
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
                output_lines.append(line.rstrip() +
                                    '*%s\n' % insert_variable)
                continue
            if line.startswith('$OMEGA'):
                for j in range(1, level):
                    output_lines.append('(-10, 0.1, 10)\n')
            output_lines.append(line)

    with open(new_script, 'w') as target:
        for line in output_lines:
            target.write(line)


def create_folder(path, folder_name):
    """Create folder for the scripts under current directory."""
    path = path + '\\' + folder_name
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    return path


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

    # copy base_script to script running IMPMAP
    shutil.copy(base_script, script)

    # create a template script running foce
    foce_script = 'foce_script.ctl'
    shutil.copy(base_script, foce_script)

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
            add_num_covariate(index, insert_variables,
                              covariate_median, foce_script, foce_script)
        else:
            level = len(demo_df[covariate].unique())
            add_factor_covariate(index, insert_variables, level,
                                 script, script)
            add_factor_covariate(index, insert_variables, level,
                                 foce_script, foce_script)
    return result_path, cov_file, result_file, insert_variables


def run_full(base_script, ncomp, insert_variables,
             data, covariates, update_list=None):
    """Run full covariance model.

    Run full model first with FOCE, if it converges and pull out the
    estimates and covariance matrix. (Since the FOCE results are not good,
    only EM will be used)
    Otherwise run full model with mu referrence, since EM method will
    always give you covariance matrix, this function will always produce
    thetas and covariance matrix based on the best structural model.
    If the objective function value fo EM method is within 5% error of the
    objective function value from FOCE. Then we consider we can continue to
    use the covariance matrix from EM method.
    """
    if update_list is None:
        update_list = []
    num_del = len(update_list)
    foce_script = 'foce_script_%s.ctl' % num_del
    foce_result_path, foce_cov_file, foce_result_file, _ = \
        full_model(base_script, foce_script, ncomp,
                   data, covariates, IMPMAP=False)

    model_update(foce_script, foce_script, update_list,
                 pattern='0 FIXED\n', var='THETA')

    run_model(foce_script, foce_result_path)
    # out = read_result_file(foce_result_file)
    # converge = (-1000000001 in out.ITERATION.values)
    theta_num = count_theta(foce_script)
    foce = False
    # if converge:
    #    thetas = read_theta(foce_result_file, ncomp, theta_num)
    #    cov = read_cov(foce_cov_file, ncomp, theta_num)
    #    foce = True
    #    return thetas, cov, foce_script, foce
    # else:
    mu_script = 'mu_script_%s.ctl' % num_del
    mu_result_path, mu_cov_file, mu_result_file, _ = \
        full_model(base_script, mu_script,
                   ncomp, data, covariates, IMPMAP=True)

    model_update(mu_script, mu_script, update_list,
                 pattern='0 FIXED\n', var='THETA')

    run_model(mu_script, mu_result_path)
    imp_result_file = sep_ext(mu_result_file)[0]
    thetas = read_theta(imp_result_file, ncomp, theta_num)
    cov = read_cov(mu_cov_file, ncomp, theta_num)
    foce_obj = read_obj(foce_result_file)
    mu_obj = read_obj(sep_ext(mu_result_file)[1])
    if abs(mu_obj - foce_obj) / foce_obj < .05:
        return thetas, cov, mu_script, foce
    else:
        return None, None, None, None


def run_full_FOCE(base_script, ncomp, insert_variables,
                  data, covariates, update_list=None):
    """Run full covariance model.

    Run full model first with FOCE, if it converges and pull out the
    estimates and covariance matrix. (Since the FOCE results are not good,
    only EM will be used)
    Otherwise run full model with mu referrence, since EM method will
    always give you covariance matrix, this function will always produce
    thetas and covariance matrix based on the best structural model.
    If the objective function value fo EM method is within 5% error of the
    objective function value from FOCE. Then we consider we can continue to
    use the covariance matrix from EM method.
    """
    if update_list is None:
        update_list = []
    num_del = len(update_list)
    foce_script = 'foce_script_%s.ctl' % num_del
    foce_result_path, foce_cov_file, foce_result_file, _ = \
        full_model(base_script, foce_script, ncomp,
                   data, covariates, IMPMAP=False)

    model_update(foce_script, foce_script, update_list,
                 pattern='0 FIXED\n', var='THETA')

    run_model(foce_script, foce_result_path)
    out = read_result_file(foce_result_file)
    converge = (-1000000001 in out.ITERATION.values)
    theta_num = count_theta(foce_script)
    foce = False
    if converge:
        thetas = read_theta(foce_result_file, ncomp, theta_num)
        cov = read_cov(foce_cov_file, ncomp, theta_num)
        foce = True
        return thetas, cov, foce_script, foce
    else:
        mu_script = 'mu_script_%s.ctl' % num_del
        mu_result_path, mu_cov_file, mu_result_file, _ = \
            full_model(base_script, mu_script,
                       ncomp, data, covariates, IMPMAP=True)

        model_update(mu_script, mu_script, update_list,
                     pattern='0 FIXED\n', var='THETA')

        run_model(mu_script, mu_result_path)
        imp_result_file = sep_ext(mu_result_file)[0]
        thetas = read_theta(imp_result_file, ncomp, theta_num)
        cov = read_cov(mu_cov_file, ncomp, theta_num)
        foce_obj = read_obj(foce_result_file)
        mu_obj = read_obj(sep_ext(mu_result_file)[1])
        if abs(mu_obj - foce_obj) / foce_obj < .05:
            return thetas, cov, mu_script, foce
        else:
            return None, None, None, None


def count_tv(script):
    """Count typical value for NONMEM script."""
    tv_num = 0
    with open(script, 'r') as target:
        for line in target:
            if line.startswith('TV'):
                tv_num += 1
    return tv_num


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


def read_coi(coi_file):
    """Read the covariance matrix from .cov file."""
    # the order is theta, sigma, omega
    data_frame = pd.read_csv(coi_file,
                             sep=' ',
                             skipinitialspace=True,
                             skiprows=1)
    del data_frame['NAME']
    mat = data_frame.values
    index_list = pd.unique(np.where(data_frame != 0)[0])
    mat = mat[np.ix_(index_list, index_list)]
    return mat


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


def modify_foce(script):
    """Deactivate important sampling option."""
    num = 0
    output_lines = []
    with open(script, 'r') as target:
        for line in target:
            if num == 1:
                output_lines.append(';' + line)
                num = 0
                continue
            if line.startswith('$ESTIMATION'):
                num += 1
            output_lines.append(line)

    with open(script, 'w') as target:
        for line in output_lines:
            target.write(line)


def sim_base_script_old(script, insert_variables, ncomp, sim_script):
    """Create simulation base script for score test."""
    shutil.copy(script, sim_script)
    for i, index in enumerate(insert_variables):
        variable_type = insert_variables[index][4]
        if variable_type == 'factor':
            model_update(sim_script, sim_script, [i+2*ncomp],
                         '0 FIXED\n', 'THETA')
        else:
            model_update(sim_script, sim_script, [i+2*ncomp],
                         '0 FIXED\n', 'THETA')
    output_lines = []
    with open(sim_script, 'r') as target:
        for line in target:
            if line.startswith('('):
                est = re.findall('\d+.\d+', line)[-1]
                output_lines.append(est + ' FIXED\n')
            elif re.match('^\d+.\d+', line):
                output_lines.append(line.rstrip() + ' FIXED\n')
            elif 'MAXEVAL=' in line:
                output_lines.append('$ESTIMATION METHOD=1 INTERACTION' +
                                    ' MAXEVAL=0 PRINT=5 NOABORT\n')
            elif 'METHOD=IMPMAP' in line or '$COV' in line:
                continue
            else:
                output_lines.append(line)
    with open(sim_script, 'w') as target:
        for line in output_lines:
            target.write(line)
    return sim_script


def sim_base_script(script, ncomp, sim_script):
    """Create simulation base script."""
    shutil.copy(script, sim_script)
    output_lines = []
    with open(sim_script, 'r') as target:
        for line in target:
            if line.startswith('('):
                est = re.findall('\d+.\d+', line)[-1]
                output_lines.append(est + ' FIXED\n')
            elif re.match('^\d+.\d+', line):
                output_lines.append(line.rstrip() + ' FIXED\n')
            elif 'MAXEVAL=' in line:
                output_lines.append('$ESTIMATION METHOD=1 INTERACTION' +
                                    ' MAXEVAL=0 PRINT=5 NOABORT\n')
            elif 'METHOD=IMPMAP' in line or '$COV' in line:
                continue
            else:
                output_lines.append(line)
    with open(sim_script, 'w') as target:
        for line in output_lines:
            target.write(line)
    return sim_script


def find_est(script, index, var_name):
    """Find the estimate and fix it to do simulation."""
    with open(script, 'r') as target:
        num = 0
        in_var = False
        for line in target:
            if in_var and line.startswith('$'):
                in_var = False
            if in_var:
                if num == index:
                    if re.match('\d+.\d+', line):
                        est = re.findall('\d+.\d+', line)[0]
                    else:
                        est = re.findall('\d+', line)[0]
                num += 1
            if line.startswith('$'+var_name):
                in_var = True
    return est


def unfix_list(script, var_name):
    """Find the list of index for a variable that has not been fixed."""
    with open(script, 'r') as target:
        num = 0
        in_var = False
        index_list = []
        for line in target:
            if in_var and line.startswith('$'):
                in_var = False
            if in_var:
                if 'FIXED' not in line:
                    index_list.append(num)
                num += 1
            if line.startswith('$'+var_name):
                in_var = True
    return index_list


def powerset(seq):
    """Generate a power set."""
    r = [[]]
    for e in seq:
        r += [x+[e] for x in r]
    return r


def wam(thetas, cov, base_theta_num, script, full_result_file,
        insert_variables, nobs, num_max=10):
    """Do variable selection by Wald method."""
    theta_num = count_theta(script)
    p = theta_num + base_theta_num + 2 - count_fixed(script)
    n = nobs
    time0 = time.time()
    non_select_list = powerset(range(0, theta_num-base_theta_num))
    sbc_list = []
    time1 = time.time()
    with open('generate_powerset.txt', 'w') as target:
        target.write(str(time1-time0))
    for the_set in non_select_list:
        q = len(the_set)
        the_sbc, the_lrt = sbc_wam(thetas, cov, the_set, p, q, n)
        sbc_list.append((the_set, the_sbc, the_lrt))
    time2 = time.time()
    with open('sbc_caclulation.txt', 'w') as target:
        target.write(str(time2-time1))
    max_sbc_list = []
    for i in range(0, num_max):
        the_max = max(sbc_list, key=operator.itemgetter(1))
        max_sbc_list.append(the_max)
        sbc_list.remove(the_max)

    wam_list = []
    real_sbc_list = []
    full_obj = read_obj(full_result_file)
    for i, the_tuple in enumerate(max_sbc_list):
        update_list = [x+base_theta_num for x in the_tuple[0]]
        the_script = 'wald_script%d.ctl' % i
        result_path = 'wald_script%d.out' % i
        result_file = 'wald_script%d.ext' % i

        model_update(script, the_script, update_list,
                     pattern='0 FIXED\n', var='THETA')

        run_model(the_script, result_path)
        reduce_obj = read_obj(result_file)
        lrt = reduce_obj - full_obj
        q = len(update_list)
        sbc = -lrt - (p-q)*math.log(n)
        real_sbc_list.append(sbc)
        index_list = [x+1 for x in update_list]
        select_index_list = range(1, theta_num+1)
        for x in index_list:
            select_index_list.remove(x)
        the_tuple = (i+1, select_index_list) + the_tuple[1:3] + (sbc, lrt)
        wam_list.append(the_tuple)

    real_rank = np.array(real_sbc_list).argsort().argsort()
    for i in range(0, num_max):
        the_rank = (num_max - real_rank[i])
        wam_list[i] += (the_rank,)

    return wam_list


def new_score_test(sim_script, paras, ncomp, base_result_file,
                   quantile=0.90):
    """New score test method.

    Use score function and observed Fisher information matrix,
    Use forward selection to do variable screening.
    """
    # number of thetas
    paras = np.array(ncomp*2*[0] + list(paras))
    size = len(paras)
    # create the fisher information matrix
    obs_fim = np.zeros((size, size))
    # create the selected covariate list
    select_list = []
    # the full theta list
    full_list = range(0, size)
    full_list = [x for x in full_list if x not in range(2*ncomp)]
    # the non_select_list
    non_select_list = [x for x in full_list if x not in select_list]
    # create a list to store the lrt_list
    select_lrt_list = []
    # set the var name to be theta
    var_name_list = dict()
    for i in range(0, size):
        var_name_list[i] = (i, 'THETA')

    # best lrt, initially it is 0 compared with base model_
    best_lrt = 0

    # forward selection begins
    step = 1
    while True:
        lrt_list = []
        if len(select_list) == 0:
            for index in range(2*ncomp, size):
                # the index is the covariate-parameter index
                obs_fim[index, index] =\
                    -second_derive(sim_script, base_result_file,
                                   index, var_name_list)
                lrt = paras[index] * paras[index] / obs_fim[index, index]
                lrt_list.append((lrt, index))
            print lrt_list
            max_lrt, max_index = max(lrt_list, key=operator.itemgetter(0))
            if max_lrt - best_lrt > stats.chi2.ppf(quantile, 1):
                select_list.append(max_index)
                non_select_list =\
                    [x for x in non_select_list if x not in select_list]
                best_lrt = max_lrt
                select_lrt_list.append(max_lrt)
            else:
                print "There is no significant covariate in screening step."
                break
        elif len(non_select_list) == 0:
            print "There is no insignificant covariate in screening step."
            select_list = [x+1 for x in select_list]
            print 'The selected list is ', select_list
            break
        else:
            for index_non_select in non_select_list:
                last_select_index = max_index
                obs_fim[index_non_select, last_select_index] =\
                    -mix_second_derive(sim_script,
                                       last_select_index,
                                       index_non_select,
                                       var_name_list)
                obs_fim[last_select_index, index_non_select] =\
                    obs_fim[index_non_select, last_select_index]
                temp_select_list = select_list + [index_non_select]
                paras_redu = paras[temp_select_list]
                fim_redu = obs_fim[np.ix_(temp_select_list, temp_select_list)]
                inv_fim = np.linalg.inv(fim_redu)
                lrt = paras_redu.dot(inv_fim).dot(paras_redu.transpose())
                lrt_list.append((lrt, index_non_select))
            max_lrt, max_index = max(lrt_list, key=operator.itemgetter(0))
            print lrt_list
            if max_lrt - best_lrt > stats.chi2.ppf(quantile, 1):
                select_list.append(max_index)
                non_select_list =\
                    [x for x in full_list if x not in select_list]
                best_lrt = max_lrt
                select_lrt_list.append(max_lrt)
            else:
                print 'The forward screening ends at step %d.' % step
                select_list = [x+1 for x in select_list]
                print 'The selected list is,', select_list
                with open('score_sw_set.txt', 'w') as target:
                    target.write(str(select_list))
                break
        step += 1
    return select_list, select_lrt_list


def true_new_score_test(sim_script, paras, ncomp, base_result_file,
                        quantile=0.90, num_omega=2, num_error=2):
    """New score test method.

    Use score function and observed Fisher information matrix,
    Use forward selection to do variable screening.
    """
    # number of thetas
    paras = np.array(ncomp*2*[0] + list(paras))
    theta_num = count_theta(sim_script)
    num_para = theta_num + num_omega + num_error
    theta_range = range(0, theta_num)
    sigma_range = range(theta_num, theta_num+num_error)
    var_name_list = dict()
    for i in range(0, num_para):
        if i in theta_range:
            var_name_list[i] = (i, 'THETA')
        elif i in sigma_range:
            var_name_list[i] = (i-theta_num, 'SIGMA')
        else:
            var_name_list[i] = (i-theta_num-num_error, 'OMEGA')

    # base model fim
    base_obs_fim = read_coi(base_result_file.replace('ext', 'coi'))
    # create the selected covariate list
    select_list = []
    # the full theta list
    full_list = range(0, theta_num)
    full_list = [x for x in full_list if x not in range(2*ncomp)]
    # create a list to store the lrt_list
    select_lrt_list = []

    # forward selection begins
    lrt_list = []
    for index in range(2*ncomp, theta_num):
        # the index is the covariate-parameter index
        obs_fim_value =\
            -second_derive(sim_script, base_result_file,
                           index, var_name_list)
        b = []
        b_range =\
            range(2*ncomp) + range(theta_num, theta_num+num_omega+num_error)
        for b_index in b_range:
            value = -mix_second_derive(sim_script, index,
                                       b_index, var_name_list)
            b.append(value)
        b = np.array(b)
        partial_obs_fim_value = obs_fim_value -\
            b.dot(np.linalg.inv(base_obs_fim)).dot(b.transpose())
        lrt = paras[index] * paras[index] / partial_obs_fim_value
        print paras[index] * paras[index] / obs_fim_value
        lrt_list.append((lrt, index))
        print lrt
        if lrt > stats.chi2.ppf(quantile, 1):
            select_list.append(index)
            select_lrt_list.append(lrt)
    select_list = [x+1 for x in select_list]
    print 'The selected list is,', select_list
    with open('score_sw_set.txt', 'w') as target:
        target.write(str(select_list))
    with open('score_chi_square.txt', 'w') as target:
        target.write(str(select_lrt_list))
    return select_list, select_lrt_list


def covariate_model_step(foce_script, base_script,
                         cwdir, insert_variables, data, nobs, covariates,
                         in_quantile, out_quantile, select_list=[]):
    """Stepwise method in PK-PD.

    Add covariate to the structural model by forward selection. Then
    do backward elimination.

    'script' is the structural model script.
    'base_script' is NONMEME script without any covariates.
    'cwdir' is the current working directory.
    'nobs' is the number of observation in the dataset.
    'data' is the dataset.
    'nobs' is the number of observation
    'covariates' are the potential covariates.
    'in_quantile' is the p value used to do forward selection.
    'out_quantile' is the p value used to backward selection.
    'select_list' is the list that contains covariates already selected.

    It returns the best model script and result.
    """
    run_count = 0
    base_theta_num = count_theta(base_script)
    theta_num = count_theta(foce_script)
    omega_num = count_omega(foce_script)

    # initialize the best script and result file
    best_script = 'base_sw.ctl'
    best_result_file = 'base_sw.ext'
    best_result_path = 'base_sw.out'

    # calculate the number of thetas
    # create a full set of list
    full_set = range(1, theta_num-base_theta_num+1)

    # create update list for those covariates already selected
    update_list = [x+base_theta_num-1 for x in full_set
                   if x+base_theta_num not in select_list]

    # update the template according to the covariates selected
    model_update(foce_script, best_script, update_list,
                 pattern='0 FIXED\n', var='THETA')

    # run the model to get base OFV
    run_model(best_script, best_result_path)

    select = [x-base_theta_num for x in select_list]

    # forward selection begins
    step = 1
    while True:
        result_file_list = []
        result_path_list = []
        script_list = []
        print "Forward selection step %d begins!" % step

        no_select = [x for x in full_set
                     if x not in select]
        if len(no_select) == 0:
            print "There is no covariate to be added."
            break
        for i, index in enumerate(no_select):
            pk_parameter = insert_variables[index][0]
            insert_variable = insert_variables[index][2]
            covariate = insert_variables[index][3]
            variable_type = insert_variables[index][4]
            print 'Adding', variable_type, covariate, 'on', pk_parameter
            result_path_i = 'model_covariate_add_%d_%d.out' % (step, i)
            result_file_i = 'model_covariate_add_%d_%d.ext' % (step, i)
            control_file_i = 'model_covariate_add_%d_%d.ctl' % (step, i)

            model_update(best_script, control_file_i, [index+base_theta_num-1],
                         pattern='(-10, 0.01, 10)\n', var='THETA')

            run_model(control_file_i, result_path_i)
            run_count += 1
            script_list.append(control_file_i)
            result_path_list.append(result_path_i)
            result_file_list.append(result_file_i)
            obj = read_obj(result_file_i)
            if obj == 999999999:
                print 'Step %d model adding %s %s on %s crashes.' %\
                    (step, variable_type, covariate, pk_parameter)
                shutil.copy(control_file_i, 'crash_'+control_file_i)
            else:
                print 'Objective function of step ' + \
                    '%d model adding %s %s on %s is %.4f' \
                    % (step, variable_type, covariate, pk_parameter, obj)

        # to avoid infinite loop
        run_time = 1

        converge = False
        while result_file_list is not None and converge is False:
            print 'The number of results in the list is', len(result_file_list)
            best_index = best_model_covariate(result_file_list)
            temp_best_result_file = 'model_covariate_add_%d_%d.ext' % \
                (step, best_index)
            temp_best_script = 'model_covariate_add_%d_%d.ctl' % \
                (step, best_index)
            temp_best_result_path = 'model_covariate_add_%d_%d.out' % \
                (step, best_index)
            dof = degree_free(temp_best_script, best_script)
            f_ss = model_compare(temp_best_result_file, best_result_file,
                                 in_quantile, dof)
            if f_ss is False:
                print 'There is no improvement for step %s.' % step
                break

            converge, cv_error, _, code =\
                model_converge(temp_best_result_file)
            # force every model to be converged
            converge = True
            if converge:
                print 'The model converges.'
                best_script = temp_best_script
                best_result_file = temp_best_result_file
                obj = read_obj(best_result_file)
                select += [no_select[best_index]]
            else:
                print 'The model does not converge.'
                if cv_error:
                    _, name_list = find_index(cv_error)
                    for index, var_name in name_list:
                        if 'THETA' in var_name:
                            print var_name, index+1, \
                                ': coefficient variation is' +\
                                ' bigger than 50%, fix it to 0 and run it.'
                            model_update(temp_best_script, temp_best_script,
                                         [index],
                                         pattern='0 FIXED\n', var=var_name)
                            run_model(temp_best_script, temp_best_result_path)
                            run_count += 1
                        else:
                            print 'Something wrong with omega or sigma.'

                elif code == 134:
                    update_list = unfix_list(temp_best_result_file, 'OMEGA')
                    print 'Best model does not converge because it is ' +\
                        'near the boundary, set OMEGA to 0.01 ' +\
                        'and rerun the model.'
                    model_update(temp_best_script, temp_best_script,
                                 update_list, '0.01\n', 'OMEGA')
                    run_model(temp_best_script, temp_best_result_path)
                    run_count += 1
                else:
                    result_file_list.remove(temp_best_result_file)
                    shutil.copy(temp_best_script,
                                'unconverge_'+temp_best_script)
            run_time += 1
            if run_time > 5:
                print "The process is in infinite loop becasue of convergence."
                print "Make it converge to continue the process."
                best_script = temp_best_script
                best_result_file = temp_best_result_file
                obj = read_obj(best_result_file)
                select += [no_select[best_index]]
                break

        if f_ss is False:
            print 'Step %d is not improved.' % step
            break
        model_update_initial(best_script, best_result_file, base_theta_num)
        step += 1
    shutil.copy(best_script, 'best_fw_' + best_script)

    # backward elimination begins
    step = 1
    while True:
        print 'Backward elimination step %d.' % step
        if len(select) == 0:
            print 'There is nothing to eliminate.'
            break
        d_temp_covariate_list = []
        result_file_list = []
        result_path_list = []
        script_list = []
        for i, index in enumerate(select):
            pk_parameter = insert_variables[index][0]
            insert_variable = insert_variables[index][2]
            covariate = insert_variables[index][3]
            variable_type = insert_variables[index][4]
            d_temp_covariate_list.append(insert_variable)
            print 'Deleting', variable_type, covariate, \
                'on', pk_parameter
            result_path_i = 'model_covariate_del_%d_%d.out' % (step, i)
            result_file_i = 'model_covariate_del_%d_%d.ext' % (step, i)
            control_file_i = 'model_covariate_del_%d_%d.ctl'\
                % (step, i)

            model_update(best_script, control_file_i,
                         [index+base_theta_num-1],
                         pattern='0 FIXED\n', var='THETA')

            run_model(control_file_i, result_path_i)
            run_count += 1
            script_list.append(control_file_i)
            result_path_list.append(result_path_i)
            result_file_list.append(result_file_i)
            obj = read_obj(result_file_i)
            if obj == 999999999:
                print 'Step %d model deleting %s %s on %s crashes.' %\
                    (step, variable_type, covariate, pk_parameter)
                shutil.copy(control_file_i, 'crash_'+control_file_i)
            else:
                print 'Objective function of step ' + \
                    '%d model deletes %s %s on %s is %.4f' \
                    % (step, variable_type, covariate,
                        pk_parameter, obj)

        worst_index = best_model_covariate(result_file_list)
        temp_best_result_file = 'model_covariate_del_%d_%d.ext' % \
            (step, worst_index)
        temp_best_script = 'model_covariate_del_%d_%d.ctl' % \
            (step, worst_index)
        temp_best_result_path = 'model_covariate_del_%d_%d.out' % \
            (step, worst_index)
        dof = degree_free(temp_best_script, best_script)
        b_ss = model_compare(temp_best_result_file, best_result_file,
                             out_quantile, dof)
        if b_ss:
            print 'There is nothing to delete for step %s.' % step
            break

        best_script = temp_best_script
        best_result_file = temp_best_result_file
        obj = read_obj(best_result_file)
        select = [x for x in select if x != select[worst_index]]

        model_update_initial(best_script, best_result_file, base_theta_num)
        step += 1
    shutil.copy(best_script, 'best_sw_' + best_script)
    best_select = range(1, 1+base_theta_num) +\
     [x+base_theta_num for x in select]
    with open('sw_set.txt', 'w') as target:
        target.write(str(best_select))

    with open('model_run_count.txt', 'w') as target:
        target.write(str(run_count))
    return best_script, best_result_file


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

def count_omega(script):
    """Count number of OMEGA."""
    omega_num = 0
    with open(script, 'r') as target:
        for line in target:
            if 'ETA(' in target:
                omega_num += 1
    return omega_num


def bw_wald(thetas, cov, script, base_script, nobs,
            insert_variables, cwdir, covariates, data,
            quantile, quantile1, options='', version=74, bw_nm=True):
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
    with open('wam_interim_result.txt', 'w') as target:
        target.write(str(update_list))
    if bw_nm:
        _, _, back_select_list = covariate_model_bw(base_script, base_theta_num,
                                                    cwdir, data, nobs,
                                                    covariates,
                                                    quantile,
                                                    select_list=update_list,
                                                    options=options,
                                                    version=74)
    return bw_list, back_select_list


def sbc_score(paras, obs_fim, the_set, theta_num, ncomp, p, q, n,
              num_error=2, omega_fixed=2):
    """Calculate sbc and lrt of the test."""
    paras_redu = paras[[x-2*ncomp for x in the_set]]
    fim_index = range(2*ncomp) + the_set +\
        range(theta_num, theta_num+2*ncomp-omega_fixed+num_error)
    sub_obs_fim = obs_fim[np.ix_(fim_index, fim_index)]
    sub_fim_index = range(2*ncomp, 2*ncomp+len(the_set))
    inv_fim = np.linalg.inv(sub_obs_fim)[np.ix_(sub_fim_index, sub_fim_index)]
    if len(the_set) > 0:
        lrt = paras_redu.dot(inv_fim).dot(paras_redu.transpose())
    else:
        lrt = 0
    sbc = lrt - q*math.log(n)
    return sbc, lrt


def sbc_score_old(paras, obs_fim, the_set, theta_num, ncomp, p, q, n):
    """Calculate sbc and lrt of the test."""
    the_set_1 = [x-2*ncomp for x in the_set]
    paras_redu = paras[the_set_1]
    inv_fim = np.linalg.inv(obs_fim[np.ix_(the_set, the_set)])
    if len(the_set) > 0:
        lrt = paras_redu.dot(inv_fim).dot(paras_redu.transpose())
    else:
        lrt = 0
    sbc = (lrt - (q)*math.log(n))
    return sbc, lrt


def grd_all(ncomp, sim_script, insert_variables):
    """Calculate the gradient of all parameters."""
    # grd_list = []
    theta_num = count_theta(sim_script)
    grd_theta = first_derive(sim_script, theta_num, 'THETA', ncomp)
    # grd_omgega = first_derive(sim_script, 2*ncomp, 'OMEGA')
    # grd_simga = first_derive(sim_script, 2, 'SIGMA')
    # grd_list = grd_theta + grd_omgega + grd_simga
    grd_array = np.asarray(grd_theta)
    return grd_array


def first_derive(script, num, var_name, ncomp=0, h=0.0001, order=1):
    """Calculate gradient for each parameter."""
    if order == 1:
        grd_list = []
        for i in range(ncomp*2, num):
            script_list = []
            result_path_list = []
            result_file_list = []
            for j in range(0, 2):
                temp_script = 'sim_first_%d.ctl' % j
                temp_result_path = 'sim_first_%d.out' % j
                temp_result_file = 'sim_first_%d.ext' % j
                script_list.append(temp_script)
                result_path_list.append(temp_result_path)
                result_file_list.append(temp_result_file)
            est = float(find_est(script, i, var_name))
            model_update(script, script_list[0], [i],
                         '%.8f FIXED\n' % (est+h), var_name)
            model_update(script, script_list[1], [i],
                         '%.8f FIXED\n' % (est-h), var_name)
            obj_list = []
            for j in range(0, 2):
                run_model(script_list[j], result_path_list[j])
                obj = read_obj(result_file_list[j])
                obj_list.append(obj)
            # print obj_list
            obj_diff = obj_list[0] - obj_list[1]
            grd = obj_diff / (-2*2*h)
            print 'derivative of', var_name, i+1, 'is', grd
            grd_list.append(grd)
    if order == 2:
        grd_list = []
        for i in range(ncomp*2, num):
            script_list = []
            result_path_list = []
            result_file_list = []
            for j in range(0, 4):
                temp_script = 'sim_first_%d.ctl' % j
                temp_result_path = 'sim_first_%d.out' % j
                temp_result_file = 'sim_first_%d.ext' % j
                script_list.append(temp_script)
                result_path_list.append(temp_result_path)
                result_file_list.append(temp_result_file)
            est = float(find_est(script, i, var_name))
            model_update(script, script_list[0], [i],
                         '%.8f FIXED\n' % (est+2*h), var_name)
            model_update(script, script_list[1], [i],
                         '%.8f FIXED\n' % (est+h), var_name)
            model_update(script, script_list[2], [i],
                         '%.8f FIXED\n' % (est-h), var_name)
            model_update(script, script_list[3], [i],
                         '%.8f FIXED\n' % (est-2*h), var_name)
            obj_list = []
            for j in range(0, 4):
                run_model(script_list[j], result_path_list[j])
                obj = read_obj(result_file_list[j])
                obj_list.append(obj)
            # print obj_list
            obj_diff =\
                (-obj_list[0] + 8*obj_list[1] - 8*obj_list[2] + obj_list[3])
            grd = obj_diff / (-2*12*h)
            print 'derivative of', var_name, i+1, 'is', grd
            grd_list.append(grd)
    return grd_list


def second_derive(script, base_result_file, index,
                  var_name_list, h=0.0001, order=1, nobuild=False):
    """Calculate second order derivative.

    script is the simulation base model script.
    """
    if order == 1:
        var_name = var_name_list[index][1]
        real_index = var_name_list[index][0]
        base_script = script
        base_result_path = base_script.replace('ctl', 'out')
        run_model(base_script, base_result_path)
        base_result_file = base_script.replace('ctl', 'ext')
        base_obj = read_obj(base_result_file)
        print 'Second derivative for %s %d.' % (var_name, real_index+1)

        script_list = []
        result_path_list = []
        result_file_list = []
        for i in range(0, 2):
            temp_script = 'sim_sec_%d.ctl' % i
            temp_result_path = 'sim_sec_%d.out' % i
            temp_result_file = 'sim_sec_%d.ext' % i
            script_list.append(temp_script)
            result_path_list.append(temp_result_path)
            result_file_list.append(temp_result_file)

        est = float(find_est(script, real_index, var_name))
        model_update(script, script_list[0], [real_index],
                     '%.8f FIXED\n' % (est+h), var_name)
        model_update(script, script_list[1], [real_index],
                     '%.8f FIXED\n' % (est-h), var_name)
        obj_list = []
        for i in range(0, 2):
            run_model(script_list[i], result_path_list[i], nobuild=nobuild)
            obj = read_obj(result_file_list[i])
            obj_list.append(obj)
        sec_der = (obj_list[0] + obj_list[1] - 2*base_obj) / (-2*h*h)
        print sec_der
    if order == 2:
        var_name = var_name_list[index][1]
        real_index = var_name_list[index][0]
        base_script = script
        base_result_path = base_script.replace('ctl', 'out')
        run_model(base_script, base_result_path)
        base_result_file = base_script.replace('ctl', 'ext')
        base_obj = read_obj(base_result_file)
        print 'Second derivative for %s %d.' % (var_name, real_index+1)

        script_list = []
        result_path_list = []
        result_file_list = []
        for i in range(0, 4):
            temp_script = 'sim_sec_%d.ctl' % i
            temp_result_path = 'sim_sec_%d.out' % i
            temp_result_file = 'sim_sec_%d.ext' % i
            script_list.append(temp_script)
            result_path_list.append(temp_result_path)
            result_file_list.append(temp_result_file)

        est = float(find_est(script, real_index, var_name))
        model_update(script, script_list[0], [real_index],
                     '%.8f FIXED\n' % (est+2*h), var_name)
        model_update(script, script_list[1], [real_index],
                     '%.8f FIXED\n' % (est+h), var_name)
        model_update(script, script_list[2], [real_index],
                     '%.8f FIXED\n' % (est-h), var_name)
        model_update(script, script_list[3], [real_index],
                     '%.8f FIXED\n' % (est-2*h), var_name)
        obj_list = []
        for i in range(0, 4):
            run_model(script_list[i], result_path_list[i], nobuild=nobuild)
            obj = read_obj(result_file_list[i])
            obj_list.append(obj)
        sec_der = (-obj_list[0] + 16*obj_list[1] - 30*base_obj +
                   16*obj_list[2] - obj_list[3]) / (-2*12*h*h)
        print sec_der
    return sec_der


def mix_second_derive(script, ind1, ind2, var_name_list, h=0.0001,
                      nobuild=False):
    """Calculate mix second order derivative."""
    var_name1 = var_name_list[ind1][1]
    var_name2 = var_name_list[ind2][1]
    index1 = var_name_list[ind1][0]
    index2 = var_name_list[ind2][0]
    print 'Calculate second partial derivative for %s %d and %s %d.' \
        % (var_name1, index1 + 1, var_name2, index2 + 1)

    script_list = []
    result_path_list = []
    result_file_list = []
    for i in range(0, 4):
        temp_script = 'sim_mix_sec_%d.ctl' % i
        temp_result_path = 'sim_mix_sec_%d.out' % i
        temp_result_file = 'sim_mix_sec_%d.ext' % i
        script_list.append(temp_script)
        result_path_list.append(temp_result_path)
        result_file_list.append(temp_result_file)

    est1 = float(find_est(script, index1, var_name1))
    est2 = float(find_est(script, index2, var_name2))
    model_update(script, script_list[0], [index1],
                 '%.8f FIXED\n' % (est1+h), var_name1)
    model_update(script_list[0], script_list[0], [index2],
                 '%.8f FIXED\n' % (est2+h), var_name2)
    model_update(script, script_list[1], [index1],
                 '%.8f FIXED\n' % (est1+h), var_name1)
    model_update(script_list[1], script_list[1], [index2],
                 '%.8f FIXED\n' % (est2-h), var_name2)
    model_update(script, script_list[2], [index1],
                 '%.8f FIXED\n' % (est1-h), var_name1)
    model_update(script_list[2], script_list[2], [index2],
                 '%.8f FIXED\n' % (est2+h), var_name2)
    model_update(script, script_list[3], [index1],
                 '%.8f FIXED\n' % (est1-h), var_name1)
    model_update(script_list[3], script_list[3], [index2],
                 '%.8f FIXED\n' % (est2-h), var_name2)

    obj_list = []
    for i in range(0, 4):
        run_model(script_list[i], result_path_list[i], nobuild=nobuild)
        obj = read_obj(result_file_list[i])
        obj_list.append(obj)
    # print obj_list
    ser_der = (obj_list[0] - obj_list[1] - obj_list[2] + obj_list[3])\
        / (-2*4*h*h)
    return ser_der


def hessian(script, base_result_file, ncomp, num_error=2):
    """Approximate hessian matrix by finite difference."""
    theta_num = count_theta(script)
    size = theta_num + 2*ncomp + num_error
    # size = theta_num
    theta_range = range(0, theta_num)
    omega_range = range(theta_num, theta_num+2*ncomp)
    var_name_list = dict()
    for i in range(0, size):
        if i in theta_range:
            var_name_list[i] = (i, 'THETA')
        elif i in omega_range:
            var_name_list[i] = (i-theta_num, 'OMEGA')
        else:
            var_name_list[i] = (i-theta_num-2*ncomp, 'SIGMA')

    mat = np.zeros((size, size))
    for i in range(0, size):
        for j in range(0, size):
            if i == j:
                if (i > 1):
                    mat[i, i] = second_derive(script, base_result_file, i,
                                              var_name_list, nobuild=True)/2.0
                else:
                    mat[i, i] = second_derive(script, base_result_file, i,
                                              var_name_list)/2.0
                print mat[i, j]
            elif i < j:
                if (i == 1 and j == 2):
                    mat[i, j] = mix_second_derive(script, i, j, var_name_list)
                else:
                    mat[i, j] = mix_second_derive(script, i, j, var_name_list,
                                                  nobuild=True)
                print mat[i, j]
    mat = mat.T + mat
    return mat


def score_test(paras, obs_fim, ncomp, script, base_result_file,
               insert_variables, nobs, obs_matrix=True):
    """Do the score test.

    paras is gradient vector.
    obs_fim is the observed fisher information matrix.
    ncomp is the number of compartment for the best base model.
    base_result_file is the result file for the best base model.
    script is the full model script of foce.
    insert_variables is the covariate that you need to insert.

    It will return all the results in score_test.
    """
    theta_num = count_theta(script)
    p = theta_num + 2*ncomp + 2 - count_fixed(script)
    # p is the number of all parameters in the model
    n = nobs
    select_list = powerset(range(0, theta_num-2*ncomp))
    sbc_list = []
    for select_set in select_list:
        q = len(select_set)
        if q > 0:
            # q is the number of covariates fixed to 0
            update_list = [x+2*ncomp for x in select_set]
            the_sbc, the_lrt = sbc_score(paras, obs_fim, update_list,
                                         theta_num,
                                         ncomp, p, q, n)
            sbc_list.append((select_set, the_sbc, the_lrt))

    num_max = 10
    max_sbc_list = []
    for i in range(0, num_max):
        if i <= len(sbc_list):
            the_max = max(sbc_list, key=operator.itemgetter(1))
            max_sbc_list.append(the_max)
            sbc_list.remove(the_max)

    print max_sbc_list
    score_list = []
    real_sbc_list = []
    base_obj = read_obj(base_result_file)
    for i, the_tuple in enumerate(max_sbc_list):
        update_list = [x+2*ncomp for x in the_tuple[0]]
        fixed_zero_list =\
            [x for x in range(2*ncomp, theta_num) if x not in update_list]
        the_script = 'score_script%d.ctl' % i
        result_path = 'score_script%d.out' % i
        result_file = 'score_script%d.ext' % i

        model_update(script, the_script, fixed_zero_list,
                     pattern='0 FIXED\n', var='THETA')

        run_model(the_script, result_path)
        reduce_obj = read_obj(result_file)
        lrt = base_obj - reduce_obj
        q = len(update_list)
        sbc = lrt - (q)*math.log(n)
        real_sbc_list.append(sbc)
        # index_list = [x+1 for x in update_list]
        # select_index_list = range(1, theta_num+1)
        # for x in index_list:
        #     select_index_list.remove(x)
        select_index_list = range(1, 2*ncomp+1) + [x+1 for x in update_list]
        the_tuple = (i+1, select_index_list) + the_tuple[1:3] + (sbc, lrt)
        score_list.append(the_tuple)

    real_rank = np.array(real_sbc_list).argsort().argsort()
    for i in range(0, len(max_sbc_list)):
        the_rank = (len(max_sbc_list) - real_rank[i])
        score_list[i] += (the_rank,)

    return score_list


def score_test_inv(paras, obs_fim, ncomp, script, base_result_file,
                   insert_variables, nobs, obs_matrix=True):
    """Do the score test.

    paras is gradient vector.
    obs_fim is the observed fisher information matrix.
    ncomp is the number of compartment for the best base model.
    base_result_file is the result file for the best base model.
    script is the full model script of foce.
    insert_variables is the covariate that you need to insert.

    It will return all the results in score_test.
    """
    theta_num = count_theta(script)
    p = theta_num + 2*ncomp + 2 - count_fixed(script)
    # p is the number of all parameters in the model
    n = nobs
    select_list = powerset(range(0, theta_num-2*ncomp))
    sbc_list = []
    for select_set in select_list:
        q = len(select_set)
        if q > 0:
            # q is the number of covariates fixed to 0
            update_list = [x+2*ncomp for x in select_set]
            the_sbc, the_lrt = sbc_score(paras, obs_fim, update_list,
                                         theta_num,
                                         ncomp, p, q, n)
            sbc_list.append((select_set, the_sbc, the_lrt))
    #
    # with open('sbc_result.txt', 'w') as target:
    #     for line in sbc_list:
    #         target.write(line)

    num_min = 10
    min_sbc_list = []
    for i in range(0, num_min):
        if i <= len(sbc_list):
            the_min = min(sbc_list, key=operator.itemgetter(1))
            min_sbc_list.append(the_min)
            sbc_list.remove(the_min)

    print min_sbc_list
    score_list = []
    real_sbc_list = []
    base_obj = read_obj(base_result_file)
    for i, the_tuple in enumerate(min_sbc_list):
        update_list = [x+2*ncomp for x in the_tuple[0]]
        fixed_zero_list =\
            [x for x in range(2*ncomp, theta_num) if x not in update_list]
        the_script = 'score_script%d.ctl' % i
        result_path = 'score_script%d.out' % i
        result_file = 'score_script%d.ext' % i

        model_update(script, the_script, fixed_zero_list,
                     pattern='0 FIXED\n', var='THETA')

        run_model(the_script, result_path)
        reduce_obj = read_obj(result_file)
        lrt = base_obj - reduce_obj
        q = len(update_list)
        sbc = lrt - (q)*math.log(n)
        real_sbc_list.append(sbc)
        # index_list = [x+1 for x in update_list]
        # select_index_list = range(1, theta_num+1)
        # for x in index_list:
        #     select_index_list.remove(x)
        select_index_list = range(1, 2*ncomp+1) + [x+1 for x in update_list]
        the_tuple = (i+1, select_index_list) + the_tuple[1:3] + (sbc, lrt)
        score_list.append(the_tuple)

    real_rank = np.array(real_sbc_list).argsort().argsort()
    for i in range(0, len(min_sbc_list)):
        the_rank = real_rank[i]
        score_list[i] += (the_rank,)

    return score_list


def modify_positive_matrix(matrix, method=1):
    if method == 1:
        u, V = np.linalg.eig(matrix)
        print 'The original eigenvalues are ', u
        u_pos = u[u > 0]
        u_neg = u[u <= 0]
        s = sum(u_neg) / 2.0
        t = s ** 2 * 100 + 1
        for i in range(0, len(u)):
            if u[i] <= 0:
                u[i] = min(u_pos) * (s - u[i]) ** 2 / t
        print 'The modified eigenvalues are ', u
        matrix = np.dot(V, np.dot(np.diag(u), np.transpose(V)))
    elif method == 2:
        u, V = np.linalg.eig(matrix)
        print 'The original eigenvalues are ', u
        for i in range(0, len(u)):
            if u[i] <= 0:
                u[i] = 0.01
        print 'The modified eigenvalues are ', u
        matrix = np.dot(V, np.dot(np.diag(u), np.transpose(V)))
    return matrix


def RMSE_Value(error_vector):
    """Calculate RMSE."""
    return ((error_vector.dot(error_vector)) / len(error_vector))**(1.0/2)

def create_NONMEM_data_1(sample_size, sparsity, cov_data):
    """Create NONMEM simulation data for one compartment.

    sample_size is the sample size.
    sparsity is the sparsity of the simulation design.
    cov_data is the covariate data from NHANES dataset.
    """
    cov_df = pd.read_csv(cov_data)
    sample_size = pd.to_numeric(sample_size)
    sub_cov_data  = cov_df.sample(n=sample_size, replace=False)
    col_names = ['CID', 'TIME', 'EVID', 'AMT', 'RATE', 'DV', 'AGE',\
        'WT', 'SEX', 'ALB', 'ALP', 'ALT']
    if sparsity.startswith("Intensive"):
        t = [0, 2, 4, 8, 12, 24]
    else:
        t = [0] + sample(range(1, 13), 1) + sample(range(13, 25), 1)
    n = len(t)
    ID = [0] * (n*sample_size)
    TIME = [0] * (n*sample_size)
    EVID = [0] * (n*sample_size)
    AMT = [0] * (n*sample_size)
    RATE = [0] * (n*sample_size)
    DV = [0] * (n*sample_size)
    AGE = [0] * (n*sample_size)
    WT = [0] * (n*sample_size)
    SEX = [0] * (n*sample_size)
    ALB = [0] * (n*sample_size)
    ALP = [0] * (n*sample_size)
    ALT = [0] * (n*sample_size)
    for i in range(1, sample_size+1):
        ID[(i-1)*n:i*n] = [i] * n
        TIME[(i-1)*n:i*n] = t
        EVID[(i-1)*n:i*n] = [0] * n
        EVID[(i-1)*n] = 1
        AMT[(i-1)*n:i*n] = [0] * n
        AMT[(i-1)*n] = 1000
        RATE[(i-1)*n:i*n] = [0] * n
        DV[(i-1)*n:i*n] = [0] * n
    list_data = list(zip(ID, TIME, EVID, AMT, RATE, DV,
                         AGE, WT, SEX, ALB, ALP, ALT))
    df = pd.DataFrame(data=list_data, columns = col_names)
    cov_index = ['AGE', 'WT', 'SEX', 'ALB', 'ALP', 'ALT']
    for i in range(1, sample_size+1):
        for cov_index_j in cov_index:
            df.loc[df['CID'] == i, cov_index_j] =\
                cov_df.loc[i, cov_index_j]
    df.to_csv('full_sim_data.csv', index=False)

def create_NONMEM_data_2(cov_data):
    """Create NONMEM simulation data for Rituximab.

    cov_data is the covariate data from NHANES dataset.
    """
    rituximab_data = pd.read_csv('rituximab.csv')
    cov_df = pd.read_csv(cov_data)
    sample_size = len(np.unique(rituximab_data['CID']))
    sub_cov_data  = cov_df.sample(n=sample_size, replace=False)
    col_names = ['CID', 'TIME', 'EVID', 'AMT', 'RATE', 'DV', 'AGE',
        'SEX', 'BSA', 'ALB', 'ALP', 'ALT', 'AST', 'HCT',
        'LDH', 'PC', 'SCR', 'TB']
    t1 = [0, 3.0/24, 6.0/24, 48.0/24]
    t2 = [14.0, 14+3.0/24, 14+6.0/24, 14+48.0/24, 14+336.0/24,
            14+1088.0/24, 14+2352.0/24, 14+3696.0/24]
    t = t1 + t2
    n = len(t)
    ID = [0] * (n*sample_size)
    TIME = [0] * (n*sample_size)
    EVID = [0] * (n*sample_size)
    AMT = [0] * (n*sample_size)
    RATE = [0] * (n*sample_size)
    DV = [0] * (n*sample_size)
    AGE = [0] * (n*sample_size)
    SEX = [0] * (n*sample_size)
    BSA = [0] * (n*sample_size)
    ALB = [0] * (n*sample_size)
    ALP = [0] * (n*sample_size)
    ALT = [0] * (n*sample_size)
    AST = [0] * (n*sample_size)
    HCT = [0] * (n*sample_size)
    LDH = [0] * (n*sample_size)
    PC = [0] * (n*sample_size)
    SCR = [0] * (n*sample_size)
    TB = [0] * (n*sample_size)
    for i in range(1, sample_size+1):
        ID[(i-1)*n:i*n] = [i] * n
        TIME[(i-1)*n:i*n] = t
        EVID[(i-1)*n:i*n] = [0] * n
        EVID[(i-1)*n] = 1
        EVID[(i-1)*n+len(t1)] = 1
        AMT[(i-1)*n:i*n] = [0] * n
        AMT[(i-1)*n] = 1000
        AMT[(i-1)*n+len(t1)] = 1000
        RATE[(i-1)*n:i*n] = [0] * n
        RATE[(i-1)*n] = round(1000.0/(255.0/(60.0*24.0)))
        RATE[(i-1)*n+len(t1)] = round(1000.0/(195.0/(60.0*24.0)))
        DV[(i-1)*n:i*n] = [0] * n
    list_data = list(zip(ID, TIME, EVID, AMT, RATE, DV,
                         AGE, SEX, BSA, ALB, ALP, ALT, AST,
                         HCT, LDH, PC, SCR, TB))
    df = pd.DataFrame(data=list_data, columns = col_names)
    cov_index = ['AGE','SEX', 'BSA', 'ALB', 'ALP', 'ALT', 'AST', 'HCT',
        'LDH', 'PC', 'SCR', 'TB']
    for i in range(1, sample_size+1):
        for cov_index_j in cov_index:
            df.loc[df['CID'] == i, cov_index_j] =\
                cov_df.loc[i, cov_index_j]
    df.to_csv('full_sim_data.csv', index=False)

def run_one_comp(sample_size, sparsity, num_sim, home_dir):
    covariates = {'AGE': 'continuous variable',
                  'WT': 'continuous variable',
                  'SEX': 'factor',
                  'ALB': 'continuous variable',
                  'ALP': 'continuous variable',
                  'ALT': 'continuous variable'}
    # select_0 = [5, 6, 12]
    WAM_RMSE = []
    WAMBE_RMSE = []
    SW_RMSE = []
    # specify the simulation times
    i = 0
    while i < num_sim:
        try:
            # change to current work directory
            print 'Begin Simulation %d!.' % i
            project_path = home_dir + '\\One_Compartment'
            os.chdir(project_path)
            current_path = create_folder(project_path, 'Sample_Size_%s_%s' % (sample_size, sparsity))
            file_list = ['nmfe74.bat', 'data_sim_0.ctl', 'data_sim_0_Corr.ctl',
                         'base_model.ctl', 'base_model_Corr.ctl', 'pat_sim_data.csv']
            for file_ in file_list:
                shutil.copy(os.path.join(project_path, file_), current_path)

            os.chdir(current_path)

            # change the median in simulation control file
            data_demo = pd.read_csv('pat_sim_data.csv')
            wt_median = round(np.median(data_demo.WT), 3)
            outlines = []
            sim_file_0 = 'data_sim_0.ctl'
            if sparsity.endswith('Corr'):
                sim_file_0 = 'data_sim_0_Corr.ctl'
            with open(sim_file_0, 'r') as target:
                for line in target:
                    if line.startswith('CLWT') or line.startswith('VWT'):
                        old_median = float(re.findall(r'\d+\.\d+', line)[0])
                        line = line.replace(str(old_median), str(wt_median))
                    outlines.append(line)
            with open('data_sim.ctl', 'w') as target:
                for line in outlines:
                    target.write(line)

            # create NONMEM dataset
            create_NONMEM_data_1(sample_size, sparsity, 'pat_sim_data.csv')

            # run simulation in NONMEM to simulate data_sim.csv
            # using the previous template rituxan.csv
            # reset the random seed in NONMEM simulation file
            outlines = []
            with open('data_sim.ctl', 'r') as target:
                for line in target:
                    if line.startswith('$SIMULATION'):
                        line = '$SIMULATION (123456%d) ONLYSIMULATION SUBPROBLEM=1\n' % i
                    if line.startswith('$TABLE'):
                        line = '$TABLE ID TIME EVID AMT RATE DV' +\
                            ' AGE WT SEX ALB ALP ALT FILE=data_sim.csv\n'
                    outlines.append(line)

            with open('data_sim.ctl', 'w') as target:
                for line in outlines:
                    target.write(line)

            # run simulation in NONMEM
            run_model('data_sim.ctl', 'data_sim.out')

            # clean the data for future use
            data_name = 'data_sim.csv'
            data = pd.read_csv(data_name, delim_whitespace=True, skiprows=1)
            temp_data = data.copy()
            temp_data = temp_data.rename(columns={'ID': 'CID'})
            temp_data.to_csv('data_sim.csv', index=False)
###############################################################################
            os.chdir(current_path)
            data = pd.read_csv('data_sim.csv')
            nobs = num_obs(data)
            id_num = max(data.CID)
            print 'Sample size is', id_num

            # run the base model for the simulated dataset
            if sparsity.endswith('Corr'):
                shutil.copy('base_model_Corr.ctl', 'base_model.ctl')
            base_script = 'base_model.ctl'
            base_out_file = 'base_model.out'
            base_result_file = 'base_model.ext'
            run_model(base_script, base_out_file)
            out = read_result_file(base_result_file)
            base_converge = (-1000000001 in out.ITERATION.values)
            if not base_converge:
                continue

            # create a full model with all covariates script for NONMEM
            full_script = 'full_script_0.ctl'
            foce_script = "foce_script.ctl"
            full_result_path, full_cov_file, full_result_file, insert_variables =\
                full_model(base_script, full_script, data, covariates, IMPMAP=True)

            # run the full model and extract the result
            full_mod_time_0 = time.time()
            run_model(full_script, full_result_path)
            final_result_file = sep_ext(full_result_file)
            theta_num = count_theta(full_script)
            base_theta_num = count_theta(base_script)
            thetas = read_theta(final_result_file, base_theta_num, theta_num)
            cov = read_cov(full_cov_file, base_theta_num, theta_num)
            full_obj = read_obj(final_result_file)
            full_mod_time = time.time() - full_mod_time_0
###############################################################################
            os.chdir(current_path)
            cwdir = create_folder(current_path, 'test_wam_%d' % i)
            file_list = ['nmfe74.bat', 'data_sim.ctl', 'data_sim.csv',
                         'base_model.ctl', 'base_model.ext',
                         'foce_script.ctl', 'full_script_0.ctl',
                         'full_script_0.ext',
                         'final_full_script_0.ext', 'full_sim_data.csv']
            for file_ in file_list:
                shutil.copy(os.path.join(current_path, file_), cwdir)
            os.chdir(cwdir)

            wam_time_0 = time.time()
            wam_result = wam(thetas, cov, base_theta_num,
                             foce_script, final_result_file,
                             insert_variables, nobs)
            wam_time = time.time() - wam_time_0 + full_mod_time

            wam_title = [('Rank_Sim', 'Theta Selected', 'SBC_Sim', 'LAMDA_Sim',
                          'SBC_Real', 'LAMDA_Real', 'Rank_Real')]
            wam_table = wam_title + wam_result
            with open('wald_result.csv', 'wb') as csvfile:
                writer = csv.writer(csvfile)
                [writer.writerow(r) for r in wam_table]
            with open('full_model_time.txt', 'w') as target:
                target.write(str(full_mod_time))
            with open('wam_time.txt', 'w') as target:
                target.write('The script take %.4f s.' % wam_time)

            # find the best model
            wam_result = pd.read_csv('wald_result.csv')
            best_index = wam_result[wam_result.Rank_Real == 1].loc[:, 'Rank_Sim']
            best_script = 'best_' + 'wald_script%d.ctl' % (int(best_index)-1)
            shutil.copy('wald_script%d.ctl' % (int(best_index)-1),
                        best_script)
            best_result_out = best_script.replace('ctl', 'out')

            # change the result name
            outlines = []
            with open(best_script, 'r') as target:
                for line in target:
                    if line.startswith('FILE='):
                        line = 'FILE=wald_script.csv NOAPPEND NOPRINT\n'
                    outlines.append(line)
            with open(best_script, 'w') as target:
                for line in outlines:
                    target.write(line)

            # rerun the model and get the IPRED
            run_model(best_script, best_result_out)
            data_name = 'wald_script.csv'
            result_csv = pd.read_csv(data_name, delim_whitespace=True,
                                     skiprows=1)
            # result_csv = result_csv[result_csv.ID != 'ID']
            result_csv.apply(pd.to_numeric, errors='coerce')
            result_csv = result_csv.convert_objects(convert_numeric=True)
            result_csv_obs = result_csv[result_csv.EVID < 1]
            PRED = result_csv_obs['PRED']

            # change the simulation file
            outlines = []
            with open('data_sim.ctl', 'r') as target:
                for line in target:
                    if line.startswith('$SIMULATION'):
                        line = '$SIMULATION (12345) ONLYSIMULATION SUBPROBLEM=100\n'
                    if line.startswith('$TABLE'):
                        line = line.replace('data_sim', 'data_sim_RMSE')
                    outlines.append(line)
            with open('data_sim_RMSE.ctl', 'w') as target:
                for line in outlines:
                    target.write(line)
            # run the simulation file and clean the results
            run_model('data_sim_RMSE.ctl', 'data_sim_RMSE.out')
            data_name = 'data_sim_RMSE.csv'
            result_csv = pd.read_csv(data_name, delim_whitespace=True, skiprows=1)
            result_csv = result_csv[result_csv.ID != 'ID']
            result_csv = result_csv[result_csv.ID != 'TABLE']
            result_csv.apply(pd.to_numeric, errors='coerce')
            result_csv = result_csv.convert_objects(convert_numeric=True)
            result_csv_obs = result_csv[result_csv.EVID < 1]
            DV = result_csv_obs['DV']
            DIFF = np.tile(PRED, 100) - DV
            rmse = RMSE_Value(DIFF)
            WAM_RMSE.append(rmse)
###############################################################################
            quantile_dict = {1: (0.05, 0.05),
                             2: (0.05, 0.01),
                             3: (0.01, 0.05),
                             4: (0.01, 0.01),}
            for index in quantile_dict:
                os.chdir(current_path)
                p_value_1, p_value_2 = quantile_dict[index]
                # inner elimination significance level
                quantile_1 = 1 - p_value_1
                # elimination significance level in NONMEM
                quantile_2 = 1 - p_value_2
                cwdir = create_folder(current_path,
                                      'test_wald_%d_i' % i + str(index))
                for file_ in file_list:
                    shutil.copy(os.path.join(current_path, file_), cwdir)
                os.chdir(cwdir)

                # run wald method with backward elimintation
                bw_time_0 = time.time()
                base_script = 'base_model.ctl'
                bw_list, back_select_set = bw_wald(thetas, cov,
                                                   full_script, base_script,
                                                   nobs,
                                                   insert_variables,
                                                   cwdir,
                                                   covariates,
                                                   data,
                                                   quantile_2,
                                                   quantile_1)
                back_select_set = [("The best set is " + str(back_select_set), )]
                bw_time = time.time() - bw_time_0 + full_mod_time

                bw_title = [('Variable Deleted', 'LRT')]
                bw_table = bw_title + bw_list + back_select_set
                with open('backward_result.csv', 'wb') as csvfile:
                    writer = csv.writer(csvfile)
                    [writer.writerow(r) for r in bw_table]
                with open('wam_be_time.txt', 'w') as target:
                    target.write('The script take %.4f s.' % bw_time)

                # find the best model
                for file_ in os.listdir(cwdir):
                    if file_.startswith('best_bw_') and file_.endswith('.ctl'):
                        best_script = file_
                        best_result_out = best_script.replace('ctl', 'out')

                # change the result name
                outlines = []
                with open(best_script, 'r') as target:
                    for line in target:
                        if line.startswith('FILE='):
                            line = 'FILE=wambe.csv NOAPPEND NOPRINT\n'
                        outlines.append(line)
                with open(best_script, 'w') as target:
                    for line in outlines:
                        target.write(line)

                # rerun the model and get the IPRED
                run_model(best_script, best_result_out)
                data_name = 'wambe.csv'
                result_csv = pd.read_csv(data_name, delim_whitespace=True,
                                         skiprows=1)
                # result_csv = result_csv[result_csv.ID != 'ID']
                result_csv.apply(pd.to_numeric, errors='coerce')
                result_csv = result_csv.convert_objects(convert_numeric=True)
                result_csv_obs = result_csv[result_csv.EVID < 1]
                PRED = result_csv_obs['PRED']
                DIFF = np.tile(PRED, 100) - DV
                rmse = RMSE_Value(DIFF)
                WAMBE_RMSE.append(rmse)
###############################################################################
            quantile_dict = {0.95: '005',
                             0.99: '001'}
            for quantile in quantile_dict:
                p_value = quantile_dict[quantile]
                # create folder for stepwise method
                os.chdir(current_path)
                cwdir = create_folder(current_path, 'test_sw_%d_' % i + p_value)
                for file_ in file_list:
                    shutil.copy(os.path.join(current_path, file_), cwdir)
                os.chdir(cwdir)

                start_time = time.time()
                final_script, final_result_file = \
                    covariate_model_step(foce_script, base_script,
                                         cwdir, insert_variables,
                                         data, nobs,
                                         covariates,
                                         quantile, out_quantile=0.99)
                end_time = time.time()
                run_time = end_time - start_time
                with open('sw_set.txt', 'a') as target:
                    target.write('The script take %.4f s.' % run_time)

                # find the best model
                for file_ in os.listdir(cwdir):
                    if file_.startswith('best_sw_') and file_.endswith('.ctl'):
                        best_script = file_
                        best_result_out = best_script.replace('ctl', 'out')

                # change the result name
                outlines = []
                with open(best_script, 'r') as target:
                    for line in target:
                        if line.startswith('FILE='):
                            line = 'FILE=sw.csv NOAPPEND NOPRINT\n'
                        outlines.append(line)
                with open(best_script, 'w') as target:
                    for line in outlines:
                        target.write(line)

                # rerun the model and get the IPRED
                run_model(best_script, best_result_out)
                data_name = 'sw.csv'
                result_csv = pd.read_csv(data_name, delim_whitespace=True,
                                         skiprows=1)
                # result_csv = result_csv[result_csv.ID != 'ID']
                result_csv.apply(pd.to_numeric, errors='coerce')
                result_csv = result_csv.convert_objects(convert_numeric=True)
                result_csv_obs = result_csv[result_csv.EVID < 1]
                PRED = result_csv_obs['PRED']
                DIFF = np.tile(PRED, 100) - DV
                rmse = RMSE_Value(DIFF)
                SW_RMSE.append(rmse)
###############################################################################
            i += 1
        except:
            print 'There is an error.'
            # continue

    WAMBE1_RMSE = [WAMBE_RMSE[index]
        for index in range(len(WAMBE_RMSE)) if index % 4 == 0]
    WAMBE2_RMSE = [WAMBE_RMSE[index]
        for index in range(len(WAMBE_RMSE)) if index % 4 == 1]
    WAMBE3_RMSE = [WAMBE_RMSE[index]
        for index in range(len(WAMBE_RMSE)) if index % 4 == 2]
    WAMBE4_RMSE = [WAMBE_RMSE[index]
        for index in range(len(WAMBE_RMSE)) if index % 4 == 3]
    SW_1_RMSE = [SW_RMSE[index]
        for index in range(len(SW_RMSE)) if index % 2 == 0]
    SW_2_RMSE = [SW_RMSE[index]
        for index in range(len(SW_RMSE)) if index % 2 == 1]

    df = pd.DataFrame(data = [SW_1_RMSE, SW_2_RMSE, WAM_RMSE,
                              WAMBE1_RMSE, WAMBE2_RMSE, WAMBE3_RMSE,
                              WAMBE4_RMSE])
    df = df.transpose()
    df.columns = ['SCM1', 'SCM2', 'WAM', 'M1', 'M2', 'M3', 'M4']
    os.chdir(current_path)
    df.to_csv(current_path + '\\one_comp_RMSE_%s.csv' % sample_size, index=False)

def run_rituximab(num_sim, home_dir):
    covariates = {'AGE': 'continuous variable',
              'SEX': 'factor',
              'BSA': 'continuous variable',
              'ALT': 'continuous variable',
              'ALP': 'continuous variable',
              'AST': 'continuous variable',
              'ALB': 'continuous variable',
              'TB': 'continuous variable',
              'HCT': 'continuous variable',
              'PC': 'continuous variable',
              'LDH':'continuous variable',
              'SCR': 'continuous variable'}
    # select_0 = [11, 16, 24, 29]
    WAM_RMSE = []
    WAMBE_RMSE = []
    SW_RMSE = []
    # specify the simulation times
    i = 0
    while i < num_sim:
        try:
            # create a new folder for wald method
            # change to current work directory
            print 'Begin Simulation %d!.' % i
            project_path = home_dir + '\\Rituximab'
            os.chdir(project_path)

            current_path = create_folder(project_path, 'Simulation')
            file_list = ['nmfe74.bat', 'data_sim_0.ctl',
                         'base_model.ctl', 'pat_sim_data.csv',
                         'rituximab.csv']
            for file_ in file_list:
                shutil.copy(os.path.join(project_path, file_), current_path)
            os.chdir(current_path)

            # create NONMEM dataset
            create_NONMEM_data_2('pat_sim_data.csv')

            # run simulation in NONMEM to simulate data_sim.csv
            # using the previous template rituxan.csv
            # reset the random seed in NONMEM simulation file
            outlines = []
            with open('data_sim_0.ctl', 'r') as target:
                for line in target:
                    if line.startswith('$SIMULATION'):
                        line = '$SIMULATION (123456%d) ONLYSIMULATION SUBPROBLEM=1\n' % i
                    outlines.append(line)

            with open('data_sim.ctl', 'w') as target:
                for line in outlines:
                    target.write(line)

            run_model('data_sim.ctl', 'data_sim.out')

            # clean data_sim.csv
            data_name = 'data_sim.csv'
            data = pd.read_csv(data_name, delim_whitespace=True, skiprows=1)
            temp_data = data.copy()
            temp_data = temp_data.rename(columns={'ID': 'CID'})
            temp_data.to_csv('data_sim.csv', index=False)
            os.chdir(current_path)
            data = pd.read_csv('data_sim.csv')
            nobs = num_obs(data)

            # run the base model for the simulated dataset
            base_script = 'base_model.ctl'
            base_out_file = 'base_model.out'
            base_result_file = 'base_model.ext'
            run_model(base_script, base_out_file)
            out = read_result_file(base_result_file)
            base_converge = (-1000000001 in out.ITERATION.values)
            if not base_converge:
                continue

            # create a full model with all covariates script for NONMEM
            full_script = 'full_script_0.ctl'
            foce_script = "foce_script.ctl"
            full_result_path, full_cov_file, full_result_file, insert_variables =\
                full_model(base_script, full_script, data, covariates, IMPMAP=True)

            # run the full model and extract the result
            full_mod_time_0 = time.time()
            run_model(full_script, full_result_path)
            final_result_file = sep_ext(full_result_file)
            theta_num = count_theta(full_script)
            base_theta_num = count_theta(base_script)
            thetas = read_theta(final_result_file, base_theta_num, theta_num)
            cov = read_cov(full_cov_file, base_theta_num, theta_num)
            full_obj = read_obj(final_result_file)
            full_mod_time = time.time() - full_mod_time_0
###############################################################################
            os.chdir(current_path)
            cwdir = create_folder(current_path, 'test_wam_%d' % i)
            file_list = ['nmfe74.bat', 'data_sim.ctl', 'data_sim.csv',
                         'base_model.ctl', 'base_model.ext',
                         'foce_script.ctl', 'full_script_0.ctl',
                         'full_script_0.ext', 'full_sim_data.csv',
                         'final_full_script_0.ext']
            for file_ in file_list:
                shutil.copy(os.path.join(current_path, file_), cwdir)
            os.chdir(cwdir)

            wam_time_0 = time.time()
            wam_result = wam(thetas, cov, base_theta_num,
                             foce_script, final_result_file,
                             insert_variables, nobs)
            wam_time = time.time() - wam_time_0 + full_mod_time

            wam_title = [('Rank_Sim', 'Theta Selected', 'SBC_Sim', 'LAMDA_Sim',
                          'SBC_Real', 'LAMDA_Real', 'Rank_Real')]
            wam_table = wam_title + wam_result
            with open('wald_result.csv', 'wb') as csvfile:
                writer = csv.writer(csvfile)
                [writer.writerow(r) for r in wam_table]
            with open('full_model_time.txt', 'w') as target:
                target.write(str(full_mod_time))
            with open('wam_time.txt', 'w') as target:
                target.write('The script take %.4f s.' % wam_time)

            # find the best model
            wam_result = pd.read_csv('wald_result.csv')
            best_index = wam_result[wam_result.Rank_Real == 1].loc[:, 'Rank_Sim']
            best_script = 'best_' + 'wald_script%d.ctl' % (int(best_index)-1)
            shutil.copy('wald_script%d.ctl' % (int(best_index)-1),
                        best_script)
            best_result_out = best_script.replace('ctl', 'out')

            # change the result name
            outlines = []
            with open(best_script, 'r') as target:
                for line in target:
                    if line.startswith('FILE='):
                        line = 'FILE=wald_script.csv NOAPPEND NOPRINT\n'
                    outlines.append(line)
            with open(best_script, 'w') as target:
                for line in outlines:
                    target.write(line)

            # rerun the model and get the IPRED
            run_model(best_script, best_result_out)
            data_name = 'wald_script.csv'
            result_csv = pd.read_csv(data_name, delim_whitespace=True,
                                     skiprows=1)
            # result_csv = result_csv[result_csv.ID != 'ID']
            result_csv.apply(pd.to_numeric, errors='coerce')
            result_csv = result_csv.convert_objects(convert_numeric=True)
            result_csv_obs = result_csv[result_csv.EVID < 1]
            PRED = result_csv_obs['PRED']

            # change the simulation file
            outlines = []
            with open('data_sim.ctl', 'r') as target:
                for line in target:
                    if line.startswith('$SIMULATION'):
                        line = '$SIMULATION (12345) ONLYSIMULATION SUBPROBLEM=100\n'
                    if line.startswith('$TABLE'):
                        line = line.replace('data_sim', 'data_sim_RMSE')
                    outlines.append(line)
            with open('data_sim_RMSE.ctl', 'w') as target:
                for line in outlines:
                    target.write(line)
            # run the simulation file and clean the results
            run_model('data_sim_RMSE.ctl', 'data_sim_RMSE.out')
            data_name = 'data_sim_RMSE.csv'
            result_csv = pd.read_csv(data_name, delim_whitespace=True, skiprows=1)
            result_csv = result_csv[result_csv.ID != 'ID']
            result_csv = result_csv[result_csv.ID != 'TABLE']
            result_csv.apply(pd.to_numeric, errors='coerce')
            result_csv = result_csv.convert_objects(convert_numeric=True)
            result_csv_obs = result_csv[result_csv.EVID < 1]
            DV = result_csv_obs['DV']
            DIFF = np.tile(PRED, 100) - DV
            rmse = RMSE_Value(DIFF)
            WAM_RMSE.append(rmse)
###############################################################################
            quantile_dict = {1: (0.05, 0.05),
                             2: (0.05, 0.01),
                             3: (0.01, 0.05),
                             4: (0.01, 0.01)}
            for index in quantile_dict:
                os.chdir(current_path)
                p_value_1, p_value_2 = quantile_dict[index]
                # inner elimination significance level
                quantile_1 = 1 - p_value_1
                # elimination significance level in NONMEM
                quantile_2 = 1 - p_value_2
                cwdir = create_folder(current_path, 'test_wald_%d_i' % i + str(index))
                for file_ in file_list:
                    shutil.copy(os.path.join(current_path, file_), cwdir)
                os.chdir(cwdir)

                # run wald method with backward elimintation
                bw_time_0 = time.time()
                base_script = 'base_model.ctl'
                bw_list, back_select_set = bw_wald(thetas, cov,
                                                   full_script, base_script,
                                                   nobs,
                                                   insert_variables,
                                                   cwdir,
                                                   covariates,
                                                   data,
                                                   quantile_2,
                                                   quantile_1)
                back_select_set = [("The best set is " + str(back_select_set), )]
                bw_time = time.time() - bw_time_0 + full_mod_time

                bw_title = [('Variable Deleted', 'LRT')]
                bw_table = bw_title + bw_list + back_select_set
                with open('backward_result.csv', 'wb') as csvfile:
                    writer = csv.writer(csvfile)
                    [writer.writerow(r) for r in bw_table]
                with open('wam_be_time.txt', 'w') as target:
                    target.write('The script take %.4f s.' % bw_time)
                # find the best model
                for file_ in os.listdir(cwdir):
                    if file_.startswith('best_bw_') and file_.endswith('.ctl'):
                        best_script = file_
                        best_result_out = best_script.replace('ctl', 'out')

                # change the result name
                outlines = []
                with open(best_script, 'r') as target:
                    for line in target:
                        if line.startswith('FILE='):
                            line = 'FILE=wambe.csv NOAPPEND NOPRINT\n'
                        outlines.append(line)
                with open(best_script, 'w') as target:
                    for line in outlines:
                        target.write(line)

                # rerun the model and get the IPRED
                run_model(best_script, best_result_out)
                data_name = 'wambe.csv'
                result_csv = pd.read_csv(data_name, delim_whitespace=True,
                                         skiprows=1)
                # result_csv = result_csv[result_csv.ID != 'ID']
                result_csv.apply(pd.to_numeric, errors='coerce')
                result_csv = result_csv.convert_objects(convert_numeric=True)
                result_csv_obs = result_csv[result_csv.EVID < 1]
                PRED = result_csv_obs['PRED']
                DIFF = np.tile(PRED, 100) - DV
                rmse = RMSE_Value(DIFF)
                WAMBE_RMSE.append(rmse)
###############################################################################
            quantile_dict = {0.95: '005',
                             0.99: '001'}
            for quantile in quantile_dict:
                p_value = quantile_dict[quantile]
                # create folder for stepwise method
                os.chdir(current_path)
                cwdir = create_folder(current_path, 'test_sw_%d_' % i + p_value)
                for file_ in file_list:
                    shutil.copy(os.path.join(current_path, file_), cwdir)
                os.chdir(cwdir)

                start_time = time.time()
                final_script, final_result_file = \
                    covariate_model_step(foce_script, base_script,
                                         cwdir, insert_variables,
                                         data, nobs,
                                         covariates,
                                         quantile, out_quantile=0.99)
                end_time = time.time()
                run_time = end_time - start_time
                with open('sw_set.txt', 'a') as target:
                    target.write('The script take %.4f s.' % run_time)
                # find the best model
                for file_ in os.listdir(cwdir):
                    if file_.startswith('best_sw_') and file_.endswith('.ctl'):
                        best_script = file_
                        best_result_out = best_script.replace('ctl', 'out')

                # change the result name
                outlines = []
                with open(best_script, 'r') as target:
                    for line in target:
                        if line.startswith('FILE='):
                            line = 'FILE=sw.csv NOAPPEND NOPRINT\n'
                        outlines.append(line)
                with open(best_script, 'w') as target:
                    for line in outlines:
                        target.write(line)

                # rerun the model and get the IPRED
                run_model(best_script, best_result_out)
                data_name = 'sw.csv'
                result_csv = pd.read_csv(data_name, delim_whitespace=True,
                                         skiprows=1)
                # result_csv = result_csv[result_csv.ID != 'ID']
                result_csv.apply(pd.to_numeric, errors='coerce')
                result_csv = result_csv.convert_objects(convert_numeric=True)
                result_csv_obs = result_csv[result_csv.EVID < 1]
                PRED = result_csv_obs['PRED']
                DIFF = np.tile(PRED, 100) - DV
                rmse = RMSE_Value(DIFF)
                SW_RMSE.append(rmse)
            i += 1
        except:
            continue

    WAMBE1_RMSE = [WAMBE_RMSE[index]
        for index in range(len(WAMBE_RMSE)) if index % 4 == 0]
    WAMBE2_RMSE = [WAMBE_RMSE[index]
        for index in range(len(WAMBE_RMSE)) if index % 4 == 1]
    WAMBE3_RMSE = [WAMBE_RMSE[index]
        for index in range(len(WAMBE_RMSE)) if index % 4 == 2]
    WAMBE4_RMSE = [WAMBE_RMSE[index]
        for index in range(len(WAMBE_RMSE)) if index % 4 == 3]
    SW_1_RMSE = [SW_RMSE[index]
        for index in range(len(SW_RMSE)) if index % 2 == 0]
    SW_2_RMSE = [SW_RMSE[index]
        for index in range(len(SW_RMSE)) if index % 2 == 1]

    df = pd.DataFrame(data = [SW_1_RMSE, SW_2_RMSE, WAM_RMSE,
                              WAMBE1_RMSE, WAMBE2_RMSE, WAMBE3_RMSE,
                              WAMBE4_RMSE])
    df = df.transpose()
    df.columns = ['SCM1', 'SCM2', 'WAM', 'M1', 'M2', 'M3', 'M4']
    os.chdir(current_path)
    df.to_csv(current_path + '\\rituximab_RMSE.csv', index=False)
