import numpy as N
import sys

question_max_retry = 0

def clean(list):
    # type(list) = list of strings (usually obtained with the ".readlines()" method)
    # removes "\n" and "\r" and empty lines from given list
    L = len(list)
    for i in range(L):
        list[L-1-i] = list[L-1-i].replace('\n','')
        list[L-1-i] = list[L-1-i].replace('\r','')
        if list[L-1-i].split() == []:
            list.pop(L-1-i)

def error_exit(error_message):
    print '\nError:', error_message
    print '         >>> EXITING NOW <<<'
    sys.exit()

def print_warning(warning_message):
    print '\nWarning:', warning_message,'\n--------'

def print_comment(comment_message):
    print '\nComment:', comment_message,'\n--------'

def raw_input_int(question):
    question_retry = bool(True)
    question_current_retry = N.int(0)
    while (question_retry and question_current_retry <= question_max_retry):
        input_string = raw_input(question)
        try:
            return N.int(input_string)
        except:
            print ' ... cannot convert "%s" to int ...' %input_string
            question_current_retry += 1
    if question_max_retry == 0:
        error_message = 'After one trial, could not convert input to int'
    else:
        error_message = 'After %s trials, could not convert any input to int' %(question_max_retry+1)
    error_exit(error_message)

def raw_input_float(question):
    question_retry = bool(True)
    question_current_retry = N.int(0)
    while (question_retry and question_current_retry <= question_max_retry):
        input_string = raw_input(question)
        try:
            return N.float(input_string)
        except:
            print ' ... cannot convert "%s" to float ...' %input_string
            question_current_retry += 1
    if question_max_retry == 0:
        error_message = 'After one trial, could not convert input to int'
    else:
        error_message = 'After %s trials, could not convert any input to int' %(question_max_retry+1)
    error_exit(error_message)

#def set_ecut_etotal(ecuts,etotals,etotal_tol_ha=0.001,prtvol=0,plot_figures=0):
def set_convergence_parameter_scalar(parameter_values,property_values,convergence_criterium,\
        convergence_type='standard',parameter_type='unknown',property_type='unknown',sort_parameter=0,prtvol=0,plot_figures=0):
    iconverged_parameter = N.float(-1)
    converged_parameter = N.float(-1)
    if sort_parameter == 1:
        [params,indices] = N.sort(parameter_values)
        parameter_values = params
        props = property_values[indices]
        property_values = props
    if prtvol > 0:
        print '===>>> Convergence of property "%s" with respect to parameter "%s" <<<===' %(property_type,parameter_type)
        print ' Criterium used for convergence                  : %s' %convergence_criterium
        print ' Highest (supposed converged) value of parameter : %s' %parameter_values[-1]
        print ' Second highest value of parameter               : %s' %parameter_values[-2]
    for iparameter in range(len(parameter_values)-1):
        if (iconverged_parameter >= 0):
            if (N.abs(property_values[iparameter] - property_values[-1]) > convergence_criterium):
                iconverged_parameter = -1
                converged_parameter = N.float(-1)
                if prtvol > 0:
                    print '... but then not converged for the higher parameter %s' %parameter_values[iparameter]
        else:
            if (N.abs(property_values[iparameter] - property_values[-1]) < convergence_criterium):
                iconverged_parameter = iparameter
                converged_parameter = N.float(parameter_values[iparameter])
                if prtvol > 0:
                    print 'Converged parameter found : %s' %converged_parameter
    if iconverged_parameter == -1:
        print '!!! No converged parameter found !!!'
    if plot_figures == 1:
        #print 'Plot the figures (not yet done ...)'
        import matplotlib.pyplot as plt
        plt.figure()
        plt.hold('on')
        plt.grid('on')
        plt.plot(parameter_values,property_values,'x-',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
        plt.title('Convergence of property "%s" with respect to parameter "%s"' %(property_type,parameter_type))
        plt.xlabel('%s' %parameter_type)
        plt.ylabel('%s' %property_type)
        plt.show()

def set_convergence_parameter_multiple_scalar(parameter_values,property_values,convergence_criterium,\
        convergence_type='standard',parameter_type='unknown',property_type='unknown',sort_parameter=0,prtvol=0,plot_figures=0):
    property_shape = N.shape(property_values)[1:]
    print property_shape,N.prod(property_shape)
    iconverged_parameter = N.float(-1)
    converged_parameter = N.float(-1)
    if sort_parameter == 1:
        [params,indices] = N.sort(parameter_values)
        parameter_values = params
        props = property_values[indices]
        property_values = props
    if prtvol > 0:
        print '===>>> Convergence of property "%s" with respect to parameter "%s" <<<===' %(property_type,parameter_type)
        print ' Criterium used for convergence                  : %s' %convergence_criterium
        print ' Highest (supposed converged) value of parameter : %s' %parameter_values[-1]
        print ' Second highest value of parameter               : %s' %parameter_values[-2]
    ref_flattened_property = property_values[-1].flatten()
    for iparameter in range(len(parameter_values)-1):
        if (iconverged_parameter >= 0):
            for iprop in range(len(property_values[iparameter].flatten())):
                if (N.abs(property_values[iparameter].flatten()[iprop] - ref_flattened_property[iprop]) > convergence_criterium):
                    iconverged_parameter = -1
                    converged_parameter = N.float(-1)
                    if prtvol > 0:
                        print '... but then not converged for the higher parameter %s' %parameter_values[iparameter]
                    break
        else:
            testconv = 1
            for iprop in range(len(property_values[iparameter].flatten())):
                if not (N.abs(property_values[iparameter].flatten()[iprop] - ref_flattened_property[iprop]) < convergence_criterium):
                    testconv = 0
                    break
            if testconv == 1:
                iconverged_parameter = iparameter
                converged_parameter = N.float(parameter_values[iparameter])
                if prtvol > 0:
                    print 'Converged parameter found : %s' %converged_parameter
    if plot_figures == 1:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.hold('on')
        plt.grid('on')
        plt.plot(parameter_values,property_values,'x-',linewidth=1.5,markersize=8.0,markeredgewidth=1.5)
        plt.title('Convergence of property "%s" with respect to parameter "%s"' %(property_type,parameter_type))
        plt.xlabel('%s' %parameter_type)
        plt.ylabel('%s' %property_type)
        plt.show()

#conveig = wavefunctions_list[ctrl.wfk_conv_indices[-1]].align()
#for iwfk in ctrl.wfk_conv_indices:
#    wfk = wavefunctions_list[iwfk]
#    if (iconverged_parameter >= 0):
#        if prtvol == 0:
#            return converged_parameter
#        elif prtvol == 1:
#            print converged_parameter
#            return converged_parameter
#        elif prtvol > 1:
#            print 'Final converged parameter : %s' %converged_parameter
#            print 'Difference property(converged_parameter) - property(highest_parameter) : %s' \
#                               %(property_values[iconverged_parameter]-property_values[-1])
#            return converged_parameter
#        else:
#            error_message = 'The argument prtvol=%s is not allowed' %prtvol
#            error_exit(error_message)
#    else:
#        if prtvol == 0:
#            return None
#        elif prtvol == 1:
#            print '-1'
#            return None
#        elif prtvol > 1:
#            warning_message = 'No converged parameter was found'
#            print_warning(warning_message)
#            return None
#        else:
#            error_message = 'The argument prtvol=%s is not allowed' %prtvol
#            error_exit(error_message)
