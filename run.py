"""Run the GBD 2021 glaucoma-related vision loss analysis."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=ROOT/'data/local')
    parser.add_argument('--output-dir', type=Path, default=ROOT/'results')
    parser.add_argument('--stage', choices=['all','historical','bapc','hfpg','tables','figures'], default='all')
    parser.add_argument('--model-cache-dir', type=Path, help='Reuse existing fitted RDS models instead of fitting those models again; use only models for the same inputs and configuration.')
    args = parser.parse_args()
    data = args.data_dir.expanduser().resolve()
    out = args.output_dir.expanduser().resolve()
    if data == out or data in out.parents or out in data.parents:
        parser.error('Input and output directories must be separate and must not contain each other.')
    if out == ROOT or out in ROOT.parents or out == ROOT/'src' or ROOT/'src' in out.parents or out == ROOT/'data' or ROOT/'data' in out.parents:
        parser.error('Output must not overwrite source code or the input data directory.')
    required = ['support/standard_weights_source.csv','support/SDI_support_1990_2021.csv',
                'support/population_support_global_2011_2050.csv','support/countries.csv',
                'support/world_map_support.geojson','hfpg/metadata/hfpg_download_manifest.json']
    missing = [str(data/name) for name in required if not (data/name).is_file()]
    if not list((data/'cause').glob('*.csv')):
        missing.append(str(data/'cause/*.csv'))
    if not missing:
        manifest = json.loads((data/'hfpg/metadata/hfpg_download_manifest.json').read_text())
        for entry in manifest['files']:
            source=(data/'hfpg'/entry['path']).resolve()
            if data/'hfpg' not in source.parents:
                parser.error('HFPG manifest paths must remain within the hfpg input directory.')
            if not source.is_file(): missing.append(str(source))
    if missing:
        parser.error('Missing input files:\n'+'\n'.join(missing)+'\nSee data/README.md for input formats and acquisition instructions.')
    stages = ['historical','bapc','hfpg','tables','figures'] if args.stage=='all' else [args.stage]
    for folder in ['data','tables','figures']:(out/folder).mkdir(parents=True,exist_ok=True)
    env = os.environ.copy()
    env.update(GLAUCOMA_INPUT=str(data), GLAUCOMA_OUTPUT=str(out),
               MPLCONFIGDIR=str(out/'.matplotlib'), LC_ALL='en_US.UTF-8')
    if 'historical' in stages:
        for name in ['standard_weights_source.csv','SDI_support_1990_2021.csv']:
            shutil.copy2(data/'support'/name,out/'data'/name)
    if args.model_cache_dir:
        cache=args.model_cache_dir.expanduser().resolve()
        if not cache.is_dir():parser.error('Model cache directory does not exist: '+str(cache))
        if cache != out/'data/bapc_models':
            target=out/'data/bapc_models';target.mkdir(exist_ok=True)
            for p in cache.glob('*.rds'):shutil.copy2(p,target/p.name)
    modules={'historical':'prepare_data.py','bapc':'bapc.R','hfpg':'hfpg.py','tables':'generate_tables.py','figures':'figures.py'}
    for stage in stages:
        module=ROOT/'src'/modules[stage]
        command=['Rscript',str(module)] if module.suffix=='.R' else [sys.executable,str(module)]
        print('Running '+stage,flush=True)
        subprocess.run(command,cwd=ROOT,env=env,check=True)
    print('Results: '+str(out))

if __name__=='__main__':main()
