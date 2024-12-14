#!/gpfs/data/yarmarkovichlab/Frank/pan_cancer/antigen_portal/spectrum_env/bin/python3.8

import sys
import os
import numpy as np
import pandas as pd
from pyteomics import mzml
from pyteomics import mgf
import matplotlib.pyplot as plt
import matplotlib as mpl
import multiprocessing as mp
import subprocess
from tqdm import tqdm
import json
import argparse 
from dash import Dash, Input, Output, callback, dash_table, html
import re
from ast import literal_eval

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus




def each_chunk_func(chunk_raws,technology):
    if technology == 'orbitrap':
        # each is a absolute path
        for raw in chunk_raws:
            f = raw.split('/')[-1]
            d = raw.split('/')[-2]
            # copy to tmpdir
            destination = os.path.join(tmp_dir,f)
            subprocess.run(['cp',raw,destination])
            # run the program and when completed, send back to mzml_dir
            cmd = 'singularity run -B {}:/data --writable {} wine msconvert /data/{} --outdir /data --zlib'.format(tmp_dir,SANDBOX_PATH,f)
            result = subprocess.run(cmd,shell=True,check=True)
            if result.returncode == 0:
                mzml_f = f.replace('.raw','.mzML')
                source = os.path.join(tmp_dir,mzml_f)
                destination = os.path.join(mzml_dir,d)
                if not os.path.exists(destination):
                    os.makedirs(destination)
                subprocess.run(['mv',source,destination])
                subprocess.run(['rm',os.path.join(tmp_dir,f)])
            else:
                print('error out when converting {}'.format(raw))
                continue

    elif technology == 'bruker':
        if config['bruker_format'] == 'hdf':  # if mgf, then you don't need to convert
            bruker_script = fixed_config['bruker_script']
            for raw in chunk_raws:
                f = raw.split('/')[-1]
                d = raw.split('/')[-2]
                destination = os.path.join(mzml_dir,d)
                if not os.path.exists(destination):
                    os.makedirs(destination)
                cmd = '{} --inpdir {} --outdir {} --running_mode convert'.format(bruker_script,raw,destination)
                subprocess.run(cmd,shell=True)
        elif config['bruker_format'] == 'mgf':  # just need to copy mgf to the mzml, please put the mgf in the folder with the .d and same name
            for raw in chunk_raws:
                f = raw.split('/')[-1]
                d = raw.split('/')[-2]
                destination = os.path.join(mzml_dir,d)
                if not os.path.exists(destination):
                    os.makedirs(destination)  
                mgf_file_path = os.path.join(raw,f.replace('.d','.mgf'))   # mgf should be within .d folder, and name should be the same of .d folder
                subprocess.run(['cp',mgf_file_path,destination])


def split_array_to_chunks(array,cores=None):
    if not isinstance(array,list):
        raise Exception('split_array_to_chunks function works for list, not ndarray')
    array_index = np.arange(len(array))
    if cores is None:
        cores = mp.cpu_count()
    sub_indices = np.array_split(array_index,cores)
    sub_arrays = []
    for sub_index in sub_indices:
        item_in_group = []
        for i in sub_index:
            item_in_group.append(array[i])
        sub_arrays.append(item_in_group)
    return sub_arrays

def annotate_ion_type(annotation, ion_types="abym"):
    if annotation.ion_type[0] in ion_types:
        if abs(annotation.isotope) == 1:
            iso = "+i" if annotation.isotope > 0 else "-i"
        elif annotation.isotope != 0:
            iso = f"{annotation.isotope:+}i"
        else:
            iso = ""
        nl = {"-NH3": "*", "-H2O": "o"}.get(annotation.neutral_loss, "")
        return f"{annotation.ion_type}{iso}{'+' * annotation.charge}{nl}"
    else:
        return ""

def draw_spectrum(mzml_path,scan,pep,d,f,rt,mass,charge,score,ma,pass_dict=None):
    if technology == 'orbitrap':
        for i,spectrum in enumerate(mzml.read(mzml_path)):
            if i == scan - 1:
                mz = spectrum['m/z array']
                intensity = spectrum['intensity array']
                rt = spectrum['scanList']['scan'][0]['scan start time']
                precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                precursor_charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']
                sus_spectrum = sus.MsmsSpectrum(spectrum['id'], precursor_mz, precursor_charge, mz, intensity, rt)
                break

        if ma == 'ITMS':
            tol_mass, tol_mode = [0.5,'Da']
        elif ma == 'FTMS':
            tol_mass, tol_mode = [20,'ppm']
        elif ma == 'TOF':
            tol_mass, tol_mode = [25,'ppm']
        elif ma == 'Unknown':
            tol_mass, tol_mode = [20,'ppm']


        spectrum = sus_spectrum
        p_fragment_tol_mass, p_fragment_tol_mode = tol_mass, tol_mode
        f_fragment_tol_mass, f_fragment_tol_mode = tol_mass, tol_mode
        peptide = pep


    elif technology == 'bruker':

        if config['bruker_format'] == 'hdf':
            # first need to use the selected frames/scans to get the dataframes in txt format in tmpdir
            bruker_script = fixed_config['bruker_script']
            cmd = '{} --running_mode read --hdf {} --fs "{}" --pep {} --tmpdir {}'.format(bruker_script,mzml_path,str(pass_dict),pep,tmp_dir)
            subprocess.run(cmd,shell=True)
            # now read in
            df = pd.read_csv(os.path.join(tmp_dir,'{}.txt'.format(pep)),sep='\t',index_col=0)   # the tmp scans
            mz_list = []
            intensity_list = []
            df['mz_values'] = [round(item) for item in df['mz_values']]  # not perfect but it is a makeshift solution
            for mz,sub_df in df.groupby(by='mz_values'):
                mz_list.append(mz)
                intensity_list.append(sub_df['intensity_values'].values.sum())

            mz = mz_list
            intensity = intensity_list
            rt = rt
            precursor_mz = mass
            precursor_charge = charge
            spectrum = sus.MsmsSpectrum(pep, precursor_mz, precursor_charge, mz, intensity, rt)

            p_fragment_tol_mass, p_fragment_tol_mode = 0.5, 'Da'
            f_fragment_tol_mass, f_fragment_tol_mode = 0.5, 'Da'
            peptide = pep

        elif config['bruker_format'] == 'mgf':
            
            pass_dict = [item.split('TITLE=')[1] for item in pass_dict]
            combined_mz = []
            combined_intensity = []

            with mgf.read(mzml_path,use_header=True) as reader:
                for title in pass_dict:
                    for spectrum in reader:
                        if spectrum.get('params').get('title') == title :
                            mz = spectrum['m/z array']
                            intensity = spectrum['intensity array']
                            combined_mz.extend(mz)
                            combined_intensity.extend(intensity)

            rt = rt
            precursor_mz = mass
            precursor_charge = charge
            spectrum = sus.MsmsSpectrum('test', precursor_mz, precursor_charge,combined_mz, combined_intensity, rt)


            p_fragment_tol_mass, p_fragment_tol_mode = 25, 'ppm'
            f_fragment_tol_mass, f_fragment_tol_mode = 25, 'ppm'
            peptide = pep



    spectrum = (
        spectrum.set_mz_range(min_mz=0, max_mz=5000)
        .remove_precursor_peak(p_fragment_tol_mass, p_fragment_tol_mode)
        .filter_intensity(min_intensity=0.0001, max_num_peaks=50000)
        .annotate_proforma(
            peptide, f_fragment_tol_mass, f_fragment_tol_mode, ion_types="aby", neutral_losses={"NH3": -17.026549, "H2O": -18.010565}
        )
    )

    fig, ax = plt.subplots(figsize=(12, 6))
    sup.spectrum(spectrum, annot_fmt=annotate_ion_type, grid=False, ax=ax)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_title('{}\n{};{};{}\n{};{};{};{}'.format(pep,d,f,scan,rt,mass,charge,round(score,2)),fontsize=10)
    if not os.path.exists(assets_dir):
        os.mkdir(assets_dir)
    plt.savefig(os.path.join(assets_dir,'spectrum_{}.png'.format(pep)), bbox_inches='tight')
    plt.close()

def hla_formatting(pre,pre_type,post_type):
    if pre_type == 'netMHCpan_output' and post_type == 'netMHCpan_input':  # HLA-A*01:01 to HLA-A01:01
        post = [hla.replace('*','') for hla in pre]
    elif pre_type == 'netMHCpan_input' and post_type == 'netMHCpan_output':  # HLA-A01:01 to HLA-A*01:01
        post = [hla[:5] + '*' + hla[5:] for hla in pre]
    elif pre_type == 'netMHCpan_output' and post_type == 'deepimmuno':  # HLA-A*01:01 to HLA-A*0101
        post = [hla.replace(':','') for hla in pre]
    elif pre_type == 'deepimmuno' and post_type == 'netMHCpan_output': # HLA-A*0101 to HLA-A*01:01
        post = []
        for hla in pre:
            first = hla[:6] # HLA-A*
            second = hla[6:8] # 01, it is so far always 2 digit, no case exceed 2 digit
            third = hla[8:] # 01 or 101, it could be 3 digit
            now = first + second + ':' + third
            post.append(now)
    elif pre_type == 'deepimmuno_nostar' and post_type == 'netMHCpan_input':   # HLA-A0101 to HLA-A01:01
        post = []
        for hla in pre:
            first = hla[:5] # HLA-A
            second = hla[5:7] # 01, it is so far always 2 digit, no case exceed 2 digit
            third = hla[7:] # 01 or 101, it could be 3 digit
            now = first + second + ':' + third
            post.append(now)
    return post


def help_draw_differential_and_intensity(typ,source,pep):
    new_json = './config_tmp_{}.json'.format(os.getpid())
    with open(template_json,'r') as f:
        config = json.load(f)
    ori_path = config['python_script']['tunable_parameters']['outdir'] 
    config['python_interact']['nuorf']['obj'] = pep
    config['python_interact']['image_format'] = 'png'
    config['python_interact']['nuorf']['immuno_dir'] = raw_dir
    config['python_interact']['nuorf']['technology'] = technology
    config['python_interact']['nuorf']['final_path'] = final_path
    config['python_interact']['nuorf']['col'] = 'percentile'

    neoverse_script = fixed_config['neoverse_script']

    if typ == 'self_gene':
        ensg,enst,symbol = source.split('|')[:3]
        config['python_interact']['gene']['markers'] = {symbol:ensg}
        with open(new_json,'w') as f:
            json.dump(config,f)
        
        # draw diff
        cmd = '{} --config {} --running_mode interactive --sub_mode plot_gene'.format(neoverse_script,new_json)
        result = subprocess.run(cmd,shell=True,check=True)
        if result.returncode == 0:
            source = os.path.join(ori_path,'{}_{}_expr_{}.png'.format(ensg,symbol,'boxplot+boxplot'))
            destination = os.path.join(assets_dir)
            subprocess.run(['mv',source,destination])

    elif typ == 'splicing':
        coords = source.split('|')[0]
        config['python_interact']['splicing']['events'] = [coords]
        with open(new_json,'w') as f:
            json.dump(config,f)
        # draw diff
        cmd = '{} --config {} --running_mode interactive --sub_mode plot_splicing'.format(neoverse_script,new_json)
        result = subprocess.run(cmd,shell=True,check=True)
        if result.returncode == 0:
            source = os.path.join(ori_path,'{}_splicing.png'.format(coords))
            destination = os.path.join(assets_dir)
            subprocess.run(['mv',source,destination])

    elif typ == 'TE_chimeric_transcript':
        coords = source.split('|')[0]
        config['python_interact']['splicing']['events'] = [coords]
        erv = source.split('|')[4].split(',')[1]
        config['python_interact']['erv']['markers'] = [erv]
        with open(new_json,'w') as f:
            json.dump(config,f)

        # draw diff 1
        cmd = '{} --config {} --running_mode interactive --sub_mode plot_splicing'.format(neoverse_script,new_json)
        result = subprocess.run(cmd,shell=True,check=True)
        if result.returncode == 0:
            source = os.path.join(ori_path,'{}_splicing.png'.format(coords))
            destination = os.path.join(assets_dir)
            subprocess.run(['mv',source,destination])

        # draw diff 2
        cmd = '{} --config {} --running_mode interactive --sub_mode plot_erv'.format(neoverse_script,new_json)
        result = subprocess.run(cmd,shell=True,check=True)
        if result.returncode == 0:
            source = os.path.join(ori_path,'{}_expr.png'.format(erv.replace('/','_')))
            destination = os.path.join(assets_dir)
            subprocess.run(['mv',source,destination])


    elif typ == 'ERV':
        erv = source.split('|')[0].split(':')[0]
        config['python_interact']['erv']['markers'] = [erv]
        with open(new_json,'w') as f:
            json.dump(config,f)
        # draw diff
        cmd = '{} --config {} --running_mode interactive --sub_mode plot_erv'.format(neoverse_script,new_json)
        result = subprocess.run(cmd,shell=True,check=True)
        if result.returncode == 0:
            source = os.path.join(ori_path,'{}_expr.png'.format(erv.replace('/','_')))
            destination = os.path.join(assets_dir)
            subprocess.run(['mv',source,destination])

    elif typ == 'intron_retention':
        intron = '|'.join(source.split('|')[:7]).replace('|',',')
        config['python_interact']['intron_retention']['introns'] = [intron]
        with open(new_json,'w') as f:
            json.dump(config,f)
        # draw diff
        cmd = '{} --config {} --running_mode interactive --sub_mode plot_ir'.format(neoverse_script,new_json)
        result = subprocess.run(cmd,shell=True,check=True)
        if result.returncode == 0:
            source = os.path.join(ori_path,'{}_expr.png'.format(intron))
            destination = os.path.join(assets_dir)
            subprocess.run(['mv',source,destination])

    elif typ == 'nuORF':
        nuorf = source.split('|')[0]
        config['python_interact']['nuorf']['t_dic'] = {cancer:1}

        with open(new_json,'w') as f:
            json.dump(config,f)
        # draw diff
        cmd = '{} --config {} --running_mode interactive --sub_mode plot_peptide_quali'.format(neoverse_script,new_json)
        result = subprocess.run(cmd,shell=True,check=True)
        if result.returncode == 0:
            source = os.path.join(ori_path,'nuorf_{}_{}.png'.format('qualitative',pep))
            destination = os.path.join(assets_dir)
            subprocess.run(['mv',source,destination])

    elif typ == 'fusion' or typ == 'variant' or typ == 'pathogen':
        with open(new_json,'w') as f:
            json.dump(config,f)


    # draw intensity

    cmd = '{} --config {} --running_mode interactive --sub_mode plot_peptide_quanti'.format(neoverse_script,new_json)
    result = subprocess.run(cmd,shell=True,check=True)
    if result.returncode == 0:
        source = os.path.join(ori_path,'{}_{}.png'.format(pep,'percentile'))
        destination = os.path.join(assets_dir)
        subprocess.run(['mv',source,destination])

        source = os.path.join(ori_path,'{}_rank_abundance.png'.format(pep))
        destination = os.path.join(assets_dir)
        subprocess.run(['mv',source,destination])

    # remove config
    subprocess.run(['rm',new_json])



@callback(Output('output_header','children'),Output('output_text','children'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))
def click_table(active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        antigen_value = data[row]['pep']
        cell_value = data[row][active_cell['column_id']]
        return 'selected_antigen is: {}'.format(antigen_value),'selected cell is: {}'.format(cell_value)
    else:
        return 'selected_antigen is: None', 'selected cell is: None'


@callback(Output('psm','src'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))
def draw_psm(active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        cell_value = data[row]['pep']  # must be pep
        pep = cell_value
        typ = antigen2typ[pep]
        source = antigen2source[pep]
        antigen_file = os.path.join(antigen_dir,'{}_neoantigen.txt'.format(typ))
        df = pd.read_csv(antigen_file,sep='\t')
        s = df.loc[df['Sequence']==pep,:].sort_values(by='Score').iloc[-1,:]
        d = s.loc['level_0']
        f = s.loc['Raw file']
        scan_number = s.loc['Scan number']
        if technology == 'orbitrap':
            rt = s.loc['Retention time']
        elif technology == 'bruker':
            rt = s.loc['Precursor retention time']
        mass = s.loc['Mass']
        charge = s.loc['Charge']
        score = s.loc['Score']
        if technology == 'orbitrap':
            ma = s.loc['Mass analyzer']
        elif technology == 'bruker':
            ma = None
        if technology == 'orbitrap':
            mzml_path = os.path.join(mzml_dir,d,'{}.mzML'.format(f))
            pass_dict = None
        elif technology == 'bruker':
            if config['bruker_format'] == 'hdf':
                mzml_path = os.path.join(mzml_dir,d,'{}.hdf'.format(f))
                precursor_ids = [int(item) for item in s.loc['PASEF precursor IDs'].split(';')]
                pasef_file = os.path.join(raw_dir,d,'combined','txt','pasefMsmsScans.txt')
                pasef = pd.read_csv(pasef_file,sep='\t')
                pass_dict = {}
                for pid in precursor_ids:
                    pasef_pid = pasef.loc[pasef['Precursor']==pid,:]
                    scans = (pasef_pid['ScanNumBegin'].iloc[0],pasef_pid['ScanNumEnd'].iloc[0])
                    pass_dict[pid] = scans
            elif config['bruker_format'] == 'mgf':
                mzml_path = os.path.join(mzml_dir,d,'{}.mgf'.format(f))
                precursor_ids = [int(item) for item in s.loc['PASEF precursor IDs'].split(';')]
                pasef_file = os.path.join(raw_dir,d,'combined','txt','pasefMsmsScans.txt')
                pasef = pd.read_csv(pasef_file,sep='\t')
                pass_dict = []   # a list of titles
                frames = []
                for pid in precursor_ids:
                    pasef_pid = pasef.loc[pasef['Precursor']==pid,:]
                    frames.extend(pasef_pid['Frame'].values.tolist())

                frames = [int(f_) for f_ in frames]
                frames = sorted(frames)

                combine = ''
                combine += str(frames[0]) 
                for i,f_ in enumerate(frames):
                    if i < len(frames) - 1:
                        if f_ + 1 == frames[i+1]:
                            continue
                        else:
                            combine += '-{};{}'.format(str(f_),str(frames[i+1]))
                    else:
                        if f_ == frames[i-1] + 1:
                            combine += '-{}'.format(str(f_))
                        else:
                            break

                cs_list = []    # [#21015-21016,#21058]
                for cs in combine.split(';'):
                    cs_list.append('#{}'.format(cs))
                

                key_a = 'TITLE'
                m_z_1, m_z_2 = str(m_z_reported).split('.')
                allowable = [int(m_z_1)-1,int(m_z_1)+1]
                key_b_pat = re.compile(r'\+PASEF\((\d+)\.')
                for cs in cs_list:   # how to match is still a problem
                    key_c = cs
                    cmd1 = 'grep "{}" {} | grep "{}"'.format(key_a,mzml_path,key_c)
                    all_rows = subprocess.run(cmd1,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
                    candidates = []
                    for row in all_rows:
                        match = re.search(key_b_pat,row)
                        if match:
                            part = int(match.groups(1)[0])
                            if part >= allowable[0] and part <= allowable[1]:
                                candidates.append(row)
                    if len(candidates) == 1:
                        pass_dict.extend(candidates)
                    else:
                        print(cmd1)
                        raise Exception('can not find spectrum')
        if not os.path.exists(os.path.join(assets_dir,'spectrum_{}.png'.format(pep))):
            draw_spectrum(mzml_path,scan_number,pep,d,f,rt,mass,charge,score,ma,pass_dict)
        return app.get_asset_url('spectrum_{}.png'.format(pep))

def draw_psm_static(pep):
    typ = antigen2typ[pep]
    source = antigen2source[pep]
    antigen_file = os.path.join(antigen_dir,'{}_neoantigen.txt'.format(typ))
    df = pd.read_csv(antigen_file,sep='\t')
    s = df.loc[df['Sequence']==pep,:].sort_values(by='Score').iloc[-1,:]
    d = s.loc['level_0']
    f = s.loc['Raw file']
    scan_number = s.loc['Scan number']
    if technology == 'orbitrap':
        rt = s.loc['Retention time']
    elif technology == 'bruker':
        rt = s.loc['Precursor retention time']
    mass = s.loc['Mass']
    charge = s.loc['Charge']
    score = s.loc['Score']
    m_z_reported = s.loc['m/z']
    if technology == 'orbitrap':
        ma = s.loc['Mass analyzer']
    elif technology == 'bruker':
        ma = None
    if technology == 'orbitrap':
        mzml_path = os.path.join(mzml_dir,d,'{}.mzML'.format(f))
        pass_dict = None
    elif technology == 'bruker':
        if config['bruker_format'] == 'hdf':
            mzml_path = os.path.join(mzml_dir,d,'{}.hdf'.format(f))
            precursor_ids = [int(item) for item in s.loc['PASEF precursor IDs'].split(';')]
            pasef_file = os.path.join(raw_dir,d,'combined','txt','pasefMsmsScans.txt')
            pasef = pd.read_csv(pasef_file,sep='\t')
            pass_dict = {}
            for pid in precursor_ids:
                pasef_pid = pasef.loc[pasef['Precursor']==pid,:]
                scans = (pasef_pid['ScanNumBegin'].iloc[0],pasef_pid['ScanNumEnd'].iloc[0])
                pass_dict[pid] = scans
        elif config['bruker_format'] == 'mgf':
            mzml_path = os.path.join(mzml_dir,d,'{}.mgf'.format(f))
            precursor_ids = [int(item) for item in s.loc['PASEF precursor IDs'].split(';')]
            pasef_file = os.path.join(raw_dir,d,'combined','txt','pasefMsmsScans.txt')
            pasef = pd.read_csv(pasef_file,sep='\t')
            pass_dict = []   # a list of titles
            frames = []
            for pid in precursor_ids:
                pasef_pid = pasef.loc[pasef['Precursor']==pid,:]
                frames.extend(pasef_pid['Frame'].values.tolist())

            frames = [int(f_) for f_ in frames]
            frames = sorted(frames)


            combine = ''
            combine += str(frames[0]) 
            for i,f_ in enumerate(frames):
                if i < len(frames) - 1:
                    if f_ + 1 == frames[i+1]:
                        continue
                    else:
                        combine += '-{};{}'.format(str(f_),str(frames[i+1]))
                else:
                    if f_ == frames[i-1] + 1:
                        combine += '-{}'.format(str(f_))
                    else:
                        break

            cs_list = []    # [#21015-21016,#21058]
            for cs in combine.split(';'):
                first,second = cs.split('-')
                first,second = int(first),int(second)
                if first == second:
                    cs = str(first)
                cs_list.append('#{}'.format(cs))

        
            key_a = 'TITLE'
            m_z_1, m_z_2 = str(m_z_reported).split('.')
            allowable = [int(m_z_1)-1,int(m_z_1)+1]
            key_b_pat = re.compile(r'\+PASEF\((\d+)\.')
            for cs in cs_list:   
                key_c = cs
                cmd1 = 'grep "{}" {} | grep "{}"'.format(key_a,mzml_path,key_c)
                all_rows = subprocess.run(cmd1,shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
                candidates = []
                if len(all_rows) == 1:  # if unique, then don't bother the mz check
                    candidates.append(all_rows[0])
                else:    # if not, hopefully, mz check can find the one
                    for row in all_rows:
                        match = re.search(key_b_pat,row)
                        if match:
                            part = int(match.groups(1)[0])
                            if part >= allowable[0] and part <= allowable[1]:
                                candidates.append(row)
                if len(candidates) == 1:
                    pass_dict.extend(candidates)
                else:
                    print('ERROR: no spectrum match for {} with frame as {}'.format(pep,cs))
                    return 1  
    draw_spectrum(mzml_path,scan_number,pep,d,f,rt,mass,charge,score,ma,pass_dict)

@callback(Output('differential_1','src'),Output('differential_2','src'),Output('intensity_1','src'),Output('intensity_2','src'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))
def draw_diffential_and_intensity(active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        cell_value = data[row]['pep']  # must be pep
        pep = cell_value
        typ = antigen2typ[pep]
        source = antigen2source[pep]
        pat = re.compile(r'ENSG(\d+)\|ENST(\d+)\|')
        if ';' in source:
            sources = source.split(';')
            for item in sources:
                match = re.search(pat,item)
                if match:
                    source = item
                    break
        if ';' in source:  # still here
            sources = source.split(';')
            for item in sources:
                if 'nuORF' not in item:
                    source = item
                    break
        if ';' in source:
            raise Exception('non canonical can not be ambiguous')
        if typ == 'self_gene':
            ensg,enst,symbol = source.split('|')[:3]   
            if (not os.path.exists(os.path.join(assets_dir,'{}_{}_expr_{}.png'.format(ensg,symbol,'boxplot+boxplot')))) or (not os.path.exists(os.path.join(assets_dir,'{}_{}.png'.format(pep,'percentile')))) or (not os.path.exists(os.path.join(assets_dir,'{}_rank_abundance.png'.format(pep)))):
                help_draw_differential_and_intensity(typ,source,pep)
            return app.get_asset_url('{}_{}_expr_{}.png'.format(ensg,symbol,'boxplot+boxplot')),None,app.get_asset_url('{}_{}.png'.format(pep,'percentile')),app.get_asset_url('{}_rank_abundance.png'.format(pep))
        elif typ == 'nuORF':
            if (not os.path.exists(os.path.join(assets_dir,'nuorf_{}_{}.png'.format('qualitative',pep)))) or (not os.path.exists(os.path.join(assets_dir,'{}_{}.png'.format(pep,'percentile')))) or (not os.path.exists(os.path.join(assets_dir,'{}_rank_abundance.png'.format(pep)))):
                help_draw_differential_and_intensity(typ,source,pep)
            return app.get_asset_url('nuorf_{}_{}.png'.format('qualitative',pep)),None,app.get_asset_url('{}_{}.png'.format(pep,'percentile')),app.get_asset_url('{}_rank_abundance.png'.format(pep))
        elif typ == 'splicing':
            coords = source.split('|')[0]
            if (not os.path.exists(os.path.join(assets_dir,'{}_splicing.png'.format(coords)))) or (not os.path.exists(os.path.join(assets_dir,'{}_{}.png'.format(pep,'percentile')))) or (not os.path.exists(os.path.join(assets_dir,'{}_rank_abundance.png'.format(pep)))):
                help_draw_differential_and_intensity(typ,source,pep)
            return app.get_asset_url('{}_splicing.png'.format(coords)),None,app.get_asset_url('{}_{}.png'.format(pep,'percentile')),app.get_asset_url('{}_rank_abundance.png'.format(pep))
        elif typ == 'TE_chimeric_transcript':
            coords = source.split('|')[0]
            erv = source.split('|')[4].split(',')[1]
            if (not os.path.exists(os.path.join(assets_dir,'{}_splicing.png'.format(coords)))) or (not os.path.exists(os.path.join(assets_dir,'{}_expr.png'.format(erv)))) or (not os.path.exists(os.path.join(assets_dir,'{}_{}.png'.format(pep,'percentile')))) or (not os.path.exists(os.path.join(assets_dir,'{}_rank_abundance.png'.format(pep)))):
                help_draw_differential_and_intensity(typ,source,pep)
            return app.get_asset_url('{}_splicing.png'.format(coords)),app.get_asset_url('{}_expr.png'.format(erv)),app.get_asset_url('{}_{}.png'.format(pep,'percentile')),app.get_asset_url('{}_rank_abundance.png'.format(pep))

        elif typ == 'ERV':
            erv = source.split('|')[0].split(':')[0]
            if (not os.path.exists(os.path.join(assets_dir,'{}_expr.png'.format(erv)))) or (not os.path.exists(os.path.join(assets_dir,'{}_{}.png'.format(pep,'percentile')))) or (not os.path.exists(os.path.join(assets_dir,'{}_rank_abundance.png'.format(pep)))):
                help_draw_differential_and_intensity(typ,source,pep)
            return  app.get_asset_url('{}_expr.png'.format(erv)),None,app.get_asset_url('{}_{}.png'.format(pep,'percentile')),app.get_asset_url('{}_rank_abundance.png'.format(pep))
        
        elif typ == 'intron_retention':
            intron = '|'.join(source.split('|')[:7]).replace('|',',')
            if (not os.path.exists(os.path.join(assets_dir,'{}_expr.png'.format(intron)))) or (not os.path.exists(os.path.join(assets_dir,'{}_{}.png'.format(pep,'percentile')))) or (not os.path.exists(os.path.join(assets_dir,'{}_rank_abundance.png'.format(pep)))):
                help_draw_differential_and_intensity(typ,source,pep)
            return  app.get_asset_url('{}_expr.png'.format(intron)),None,app.get_asset_url('{}_{}.png'.format(pep,'percentile')),app.get_asset_url('{}_rank_abundance.png'.format(pep))

        elif typ == 'fusion' or typ == 'variant' or typ == 'pathogen':  # no diff plot
            if not os.path.exists(os.path.join(assets_dir,'{}_{}.png'.format(pep,'percentile'))) or (not os.path.exists(os.path.join(assets_dir,'{}_rank_abundance.png'.format(pep)))):
                help_draw_differential_and_intensity(typ,source,pep)
            return  None,None,app.get_asset_url('{}_{}.png'.format(pep,'percentile')),app.get_asset_url('{}_rank_abundance.png'.format(pep))

def draw_diffential_and_intensity_static(pep):
    typ = antigen2typ[pep]
    source = antigen2source[pep]
    pat = re.compile(r'ENSG(\d+)\|ENST(\d+)\|')
    if ';' in source:
        sources = source.split(';')
        for item in sources:
            match = re.search(pat,item)
            if match:
                source = item
                break
    if ';' in source:  # still here
        sources = source.split(';')
        for item in sources:
            if 'nuORF' not in item:
                source = item
                break
    if ';' in source:
        raise Exception('non canonical can not be ambiguous')
    help_draw_differential_and_intensity(typ,source,pep)



@callback(Output('hla_table','data'),Input('candidate','active_cell'),Input('candidate', 'data'),Input('candidate','page_current'),Input('candidate','page_size'))   
def display_hla_table(active_cell,data,page_current,page_size):
    if active_cell:
        row = active_cell['row'] + page_current * page_size
        lists = literal_eval(data[row]['additional_query'])
        df = pd.DataFrame.from_records(lists,columns=['hla','rank_pert','nM','id'])
        hla = [item.replace('*','') for item in df['hla'].tolist()]
        df['freq'] = [hla2freq[item] for item in hla]
        df.sort_values(by='id',inplace=True)
        return df.to_dict('records')


def draw_for_each_chunk(peps):
    for pep in tqdm(peps):
        draw_psm_static(pep)
        draw_diffential_and_intensity_static(pep)



if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='NeoVerse-py')
    parser.add_argument('--config',type=str,default='',help='path to the config json')
    parser.add_argument('--running_mode',type=str,default=None,help='set_up or launch_portal')
    args = parser.parse_args()

    config_path = args.config
    running_mode = args.running_mode

    with open(config_path,'r') as f:
        all_config = json.load(f)['antigen_portal']

    config = all_config['tunable_parameters']
    fixed_config = all_config['fixed_parameters']

    us_hla = pd.read_csv(fixed_config['us_hla_path'],sep=',',index_col=0)
    us_hla.index = hla_formatting(us_hla.index.to_list(),'deepimmuno_nostar','netMHCpan_input')
    hla2freq = us_hla['Percent US population'].to_dict()

    raw_dir = config['raw_dir']
    mzml_dir = config['mzml_dir']
    technology = config['technology']



    if running_mode == 'set_up':

        tmp_dir = fixed_config['tmp_dir']
        cores = config['cores']
        SANDBOX_PATH = fixed_config['SANDBOX_PATH']

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        if technology == 'orbitrap':
            all_raws = subprocess.run('find {} -type f -name "{}"'.format(raw_dir,'*.raw'),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
        elif technology == 'bruker':
            all_raws = subprocess.run('find {} -type d -name "{}"'.format(raw_dir,'*.d'),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]

        list_of_chunk_raws = split_array_to_chunks(all_raws,cores)
        pool = mp.Pool(processes=cores)
        r = [pool.apply_async(func=each_chunk_func,args=(chunk_raws,technology,)) for chunk_raws in list_of_chunk_raws]
        pool.close()
        pool.join()

    elif running_mode == 'set_up_rescue':
        raw_dir = config['raw_dir']
        mzml_dir = config['mzml_dir']
        technology = config['technology']

        tmp_dir = fixed_config['tmp_dir']
        SANDBOX_PATH = fixed_config['SANDBOX_PATH']

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        if technology == 'orbitrap':
            all_raws = subprocess.run('find {} -type f -name "{}"'.format(raw_dir,'*.raw'),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
        elif technology == 'bruker':
            all_raws = subprocess.run('find {} -type d -name "{}"'.format(raw_dir,'*.d'),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]

        # all you should have
        compound_id_dict = {'/'.join(item.split('/')[-2:]).split('.')[0]:item for item in all_raws}
        compound_id = list(compound_id_dict.keys())


        # now check what have been succeessfully converted
        if technology == 'orbitrap':
            all_mzmls = subprocess.run('find {} -type f -name "{}"'.format(mzml_dir,'*.mzML'),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
        elif technology == 'bruker':
            if config['bruker_format'] == 'hdf':
                all_mzmls = subprocess.run('find {} -type f -name "{}"'.format(mzml_dir,'*.hdf'),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
            elif config['bruker_format'] == 'mgf':
                all_mzmls = subprocess.run('find {} -type f -name "{}"'.format(mzml_dir,'*.mgf'),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
        
        converted_compound_id = ['/'.join(item.split('/')[-2:]) for item in all_mzmls]  # d/f.mzML
        converted_compound_id = [item.split('.')[0] for item in converted_compound_id]  # d/f

        # check and rescue sequentially
        missing = set(compound_id).difference(set(converted_compound_id))
        if len(missing) == 0:
            print('done')
        else:
            all_raws = [compound_id_dict[miss] for miss in missing]
            each_chunk_func(all_raws,technology)

    elif running_mode == 'generate_figures':

        template_json = config['template_json']
        assets_dir = config['assets_dir']
        technology = config['technology']
        tmp_dir = fixed_config['tmp_dir']

        if not os.path.exists(assets_dir):
            os.makedirs(assets_dir)

        cancer = config['cancer']
        antigen_dir = config['antigen_dir']
        final_path = os.path.join(antigen_dir,'final_enhanced.txt')

        final = pd.read_csv(final_path,sep='\t')
        cond = [False if '[]' in item else True for item in final['presented_by_each_sample_hla']]
        final = final.loc[cond,:]
        final = final.loc[final['unique']!=False,:]
        selected_columns = ['pep','typ','source','detailed_intensity','highest_score','depmap_median','hla_ligand_atlas','best_pep','n_psm','presented_by_each_sample_hla','additional_query','median_tumor','max_median_gtex']
        final = final.loc[:,selected_columns]

        antigen2typ = pd.Series(index=final['pep'].values,data=final['typ'].values).to_dict()
        antigen2source = pd.Series(index=final['pep'].values,data=final['source'].values).to_dict()

        cores = config['cores']

        all_peps = final['pep'].values.tolist()
        list_of_chunk_peps = split_array_to_chunks(all_peps,cores)
        pool = mp.Pool(processes=cores)
        r = [pool.apply_async(func=draw_for_each_chunk,args=(chunk_peps,)) for chunk_peps in list_of_chunk_peps]
        pool.close()
        pool.join()

    elif running_mode == 'generate_figures_rescue':

        template_json = config['template_json']
        assets_dir = config['assets_dir']
        technology = config['technology']
        tmp_dir = fixed_config['tmp_dir']

        if not os.path.exists(assets_dir):
            os.makedirs(assets_dir)

        cancer = config['cancer']
        antigen_dir = config['antigen_dir']
        final_path = os.path.join(antigen_dir,'final_enhanced.txt')

        final = pd.read_csv(final_path,sep='\t')
        cond = [False if '[]' in item else True for item in final['presented_by_each_sample_hla']]
        final = final.loc[cond,:]
        final = final.loc[final['unique']!=False,:]
        selected_columns = ['pep','typ','source','detailed_intensity','highest_score','depmap_median','hla_ligand_atlas','best_pep','n_psm','presented_by_each_sample_hla','additional_query','median_tumor','max_median_gtex']
        final = final.loc[:,selected_columns]

        antigen2typ = pd.Series(index=final['pep'].values,data=final['typ'].values).to_dict()
        antigen2source = pd.Series(index=final['pep'].values,data=final['source'].values).to_dict()

        cores = config['cores']

        rescue_peps = []
        all_peps = final['pep'].values.tolist()
        for pep in all_peps:
            typ = antigen2typ[pep]
            source = antigen2source[pep]

            # modify source
            pat = re.compile(r'ENSG(\d+)\|ENST(\d+)\|')
            if ';' in source:
                sources = source.split(';')
                for item in sources:
                    match = re.search(pat,item)
                    if match:
                        source = item
                        break
            if ';' in source:  # still here
                sources = source.split(';')
                for item in sources:
                    if 'nuORF' not in item:
                        source = item
                        break
            if ';' in source:
                raise Exception('non canonical can not be ambiguous')

            psm_plot = os.path.join(assets_dir,'spectrum_{}.png'.format(pep))
            intensity_plot_1 = os.path.join(assets_dir,'{}_{}.png'.format(pep,'percentile'))
            intensity_plot_2 = os.path.join(assets_dir,'{}_rank_abundance.png'.format(pep))

            if typ == 'self_gene':
                ensg,enst,symbol = source.split('|')[:3]
                diff_plot = [os.path.join(assets_dir,'{}_{}_expr_{}.png'.format(ensg,symbol,'boxplot+boxplot'))]
            elif typ == 'splicing':
                coords = source.split('|')[0]
                diff_plot = [os.path.join(assets_dir,'{}_splicing.png'.format(coords))]
            elif typ == 'TE_chimeric_transcript':
                coords = source.split('|')[0]
                erv = source.split('|')[4].split(',')[1]
                diff_plot = [os.path.join(assets_dir,'{}_splicing.png'.format(coords)),os.path.join(assets_dir,'{}_expr.png'.format(erv.replace('/','_')))]
            elif typ == 'ERV':
                erv = source.split('|')[0].split(':')[0]
                diff_plot = [os.path.join(assets_dir,'{}_expr.png'.format(erv.replace('/','_')))]
            elif typ == 'intron_retention':
                intron = '|'.join(source.split('|')[:7]).replace('|',',')
                diff_plot = [os.path.join(assets_dir,'{}_expr.png'.format(intron))]
            elif typ == 'nuORF':
                nuorf = source.split('|')[0]
                diff_plot = [os.path.join(assets_dir,'nuorf_{}_{}.png'.format('qualitative',pep))]

            required_plots = [psm_plot,intensity_plot_1,intensity_plot_2] + diff_plot
            cond = all([os.path.exists(plot) for plot in required_plots])
            if not cond:
                rescue_peps.append(pep)

        if len(rescue_peps) == 0:
            print('done')
        else:
            print(rescue_peps)
            all_peps = rescue_peps
            list_of_chunk_peps = split_array_to_chunks(all_peps,cores)
            pool = mp.Pool(processes=cores)
            r = [pool.apply_async(func=draw_for_each_chunk,args=(chunk_peps,)) for chunk_peps in list_of_chunk_peps]
            pool.close()
            pool.join()



    elif running_mode == 'launch_portal':

        template_json = config['template_json']
        assets_dir = config['assets_dir']
        typs = config['type']
        technology = config['technology']
        tmp_dir = fixed_config['tmp_dir']

        if not os.path.exists(assets_dir):
            os.makedirs(assets_dir)
        
        # build the app
        cancer = config['cancer']
        antigen_dir = config['antigen_dir']
        final_path = os.path.join(antigen_dir,'final_enhanced.txt')

        final = pd.read_csv(os.path.join(antigen_dir,'final_enhanced.txt'),sep='\t')
        cond = [False if '[]' in item else True for item in final['presented_by_each_sample_hla']]
        final = final.loc[cond,:]
        final = final.loc[final['unique']!=False,:]
        selected_columns = ['pep','typ','source','detailed_intensity','highest_score','depmap_median','hla_ligand_atlas','best_pep','n_psm','presented_by_each_sample_hla','additional_query','median_tumor','max_median_gtex']
        final = final.loc[:,selected_columns]
        final = final.loc[final['typ'].isin(typs),:]
        final = final.sort_values(by=config['sort_by'],ascending=config['ascending'])

        antigen2typ = pd.Series(index=final['pep'].values,data=final['typ'].values).to_dict()
        antigen2source = pd.Series(index=final['pep'].values,data=final['source'].values).to_dict()

        # start to build app
        app = Dash(__name__,assets_folder=assets_dir)
        app.layout = html.Div([
        html.Div(html.H1('NeoVerse Antigen Portal {}'.format(cancer)),style={'text-align':'center'}),
        html.Div(dash_table.DataTable(id='candidate',
                                      data=final.to_dict('records'),
                                      page_size=20,page_current=0,page_action='native',

                                      style_data={
                                            'width': '50px', 'minWidth': '50px', 'maxWidth': '50px',
                                            'overflow': 'hidden',
                                            'textOverflow': 'ellipsis',
                                      },

                                      )),
        html.Div([html.H2(id='output_header'),html.P(id='output_text')]),
        html.Div(html.H2('psm plot')),
        html.Img(id='psm',width='50%',height='50%',style={'border-style':'dashed'}),

        html.Div([
            html.H2('differential plots'),
            html.Img(id='differential_1',width='45%',height='50%',style={'border-style':'dashed','float':'left'}),
            html.Img(id='differential_2',width='45%',height='50%',style={'border-style':'dashed','float':'right'})
        ],style={'overflow':'hidden'}),

        html.Div([
            html.H2('intensity plots'),
            html.Img(id='intensity_1',width='45%',height='50%',style={'border-style':'dashed','float':'left'}),
            html.Img(id='intensity_2',width='45%',height='50%',style={'border-style':'dashed','float':'right'}),
        ],style={'overflow':'hidden'}),

        html.Div(html.H2('HLA table')),
        html.Div(dash_table.DataTable(id='hla_table',
                                      page_size=10,page_current=0,page_action='native',   
                                      )
                )
                
            ])





        host = subprocess.run(['hostname'],stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[0]
        port = 8050
        app.run(host=host,port=port)



