


def convert_grasp(filelist):
    """
    convert all grasp files in directory using grasp2stokes
    """
    for i,f in enumerate(filelist):
        lines=[]
        g=open('parfile'+str(i),'w')
        if 'y' in f: copol='y'
        if 'x' in f: copol='x'
        lines.append('grasp_file=' + f + '\n')
        lines.append('grasp_format=' + f[-3:]+'\n')
        lines.append('grasp_copol=' + copol + '\n')
        lines.append('grasp_norm=four_pi\n')
        lines.append('stokes_file_polar=' + f[:-3] + 'stokes\n')
        g.writelines(lines)
        g.close()